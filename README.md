# Orangutan prediction pipeline (simplified)

This repository contains a lightweight workflow to:

- Prepare predictor tables from **already-aligned rasters** (alignment/clipping can be done in ArcGIS/QGIS).
- Fit a negative binomial abundance model on transect observations.
- Predict abundance on a raster-derived grid for one or more years.

The pipeline currently works with the best model only, but we can expand to full model comparison.
The logic is: 

- **Rasters (aligned) → extracted predictors in a buffer → long tables → scaled wide tables → model → predictions**

---

## Repository structure

- `predictor_config.R`
  - Single source of truth for:
    - which predictors are used everywhere (`MODEL_PREDICTORS`)
    - transforms applied before scaling (`apply_predictor_transforms()`)
    - candidate model term universe (`get_m_terms()`, `get_candidate_model_config()`)

- `prepare_predictors_from_aligned_rasters.R`
  - Main *data preparation* entrypoint.
  - Builds the prediction grid from a reference raster and extracts predictors for:
    - transects (buffered extraction)
    - grid cells (cell-level extraction)

- `helpers/scale.predictors.R`
  - Functions to scale predictors and cast long→wide for modelling/prediction.

- `best_model/`
  - `01_fit_best_model.R`: fit best NB model on transect observations
  - `02_cross_validation_leave_one_year_out.R`: leave-one-year-out CV on fitted RHS
  - `03_predict_best_model.R`: predict on grid for multiple years
  - `04_bootstrap_best_model.R`: coefficient uncertainty via parametric bootstrap

---

## Key data concepts

### 1) Long predictor table

Most scripts expect predictors in *long format*:

- `id` (unique transect id or grid cell id)
- `year` (integer)
- `predictor` (string, must match names in `MODEL_PREDICTORS`)
- `value` (numeric)

Example row:

| id | year | predictor | value |
|---|---:|---|---:|
| ground_12 | 2010 | temp_mean | 26.3 |

### 2) Geography table

The pipeline uses **geography tables** to store the geometry (coordinates) of the objects you extract/predict on.

There are two kinds of geography tables:

- `geography_observation.csv`
  - Geometry for the **observation units** used to fit the model.
  - In the current workflow these are typically **transects**.
- `geography_<year>.csv`
  - Geometry for the **prediction grid cells** for a given year.
  - In the current workflow these are **grid cells** derived from the reference raster.

All geography tables are in the model CRS (typically AEA meters) and contain:

- `id`
- `x_start`, `y_start`, `x_end`, `y_end`

The model scripts compute center coordinates from these:

- `x_center = mean(x_start, x_end)`
- `y_center = mean(y_start, y_end)`

### 3) Transect table

The abundance model uses (at minimum):

- `id`
- `year`
- `length_km`
- `nr_nests`

`nest_decay` can be provided as a column or set/hard-coded inside modelling scripts.

---

## Step-by-step workflow

### Step 0 — Align / clip rasters externally (ArcGIS/QGIS)

This pipeline assumes you already produced aligned rasters:

- same CRS
- same resolution
- same extent / alignment
- same pixel origin

You also need one raster to act as the **reference raster** for grid construction.

#### Raster naming convention

For each predictor name `p` in `MODEL_PREDICTORS`:

- Static predictor:
  - `p.tif`
- Year-specific predictor:
  - `p_<year>.tif` (preferred when available)

At runtime the preparation script tries `p_<year>.tif` first, then falls back to `p.tif`.

---

### Step 1 — Configure predictors and transforms

Edit `predictor_config.R`:

- **`MODEL_PREDICTORS`**
  - The predictor names you want to use across the entire pipeline.
  - These names must match raster filenames (see naming convention) and will become column names later.

- **`apply_predictor_transforms()`**
  - Transform rules applied to the *long predictor table* before scaling.
  - Typical transform choices (examples):
    - `log(x + 1)` for strictly-positive, right-skewed predictors (e.g. `human_pop_dens`, distance-to-feature layers)
    - `sqrt(x)` for non-negative predictors where variance increases with the mean (often percentages/areas, some distances)
  - Important: these are **heuristics**, not rules. You should validate transform choices using:
    - distribution plots (histograms/densities)
    - residual diagnostics
    - out-of-sample performance
    - domain knowledge about how the predictor should relate to the response
  - Currently implemented transforms in this repo:
    - `log(x + 1)` for `human_pop_dens`, `plantation_distance`, `pulp_distance`, `palm_distance`
    - `sqrt(x)` for `distance_PA` and `deforestation_gaveau`

If you add/remove predictors, update `MODEL_PREDICTORS` and (optionally) update transforms.

---

### Step 2 — Build grid and extract predictors (transects + grid)

Optional but recommended: if you have a study-area boundary (AOI polygon), you can provide it to the extraction script so that it only builds/extracts grid cells inside your AOI (reduces runtime and avoids ocean/outside-area predictions).

Open `prepare_predictors_from_aligned_rasters.R`, edit the **USER SETTINGS (EDIT!)** block at the top (paths, years, buffer), then click **Source** in RStudio.

Optional: the same script can also be run from the command line.

Command line example:

```bash
Rscript prepare_predictors_from_aligned_rasters.R \
  --raster-dir "/path/to/aligned_rasters" \
  --reference-raster "/path/to/reference.tif" \
  --transects-csv "/path/to/transects.csv" \
  --out-dir "/path/to/pipeline_inputs" \
  --years "1999,2000,2001,2015" \
  --buffer-km 20 \
  --mask-shp "/optional/aoi_polygon.shp"
```

#### What this produces

In `--out-dir`:

- `geography_observation.csv`
- `transects.csv`
- `predictors_observation_<buffer_km>.csv`  (long)

And for each year `Y` in `--years`:

- `geography_Y.csv`  (grid geometry)
- `predictors_abundance_Y.csv` (long)

#### Optional: AOI polygon mask (`--mask-shp`)

- If supplied, the script keeps only grid cells whose **cell centers** are inside the polygon.
- This is useful to avoid extracting/predicting over ocean or outside your study region.

---

### Step 3 — Fit the model

Edit `best_model/01_fit_best_model.R`:

- Set:
  - `INPUT_DIR` to the folder created in Step 2
  - `OUTPUT_DIR` to a model output folder

Then run the script. It will write (among others):

- `best_glm_nb_model.rds`
- `predictors_observation_scaled.csv`
- `best_model_rhs.csv`

---

### Step 4 — Predict on the grid

Edit `best_model/03_predict_best_model.R`:

- Set:
  - `INPUT_DIR` to the folder created in Step 2
  - `MODEL_DIR` to the folder created in Step 3
  - `OUTPUT_DIR` to a prediction output folder

Optionally edit:

- `YEARS_TO_PREDICT`

Then run. It writes per year:

- `predictors_grid_scaled_<year>.csv`
- `abundance_pred_per_cell_<year>.csv`

---

## How buffered extraction works (continuous vs “categorical % coverage”)

The preparation script buffers transect lines and extracts raster values inside the buffer.

### Continuous predictors

For continuous rasters (e.g. `temp_mean`):

- extract all pixel values within the buffer
- use `mean(values, na.rm=TRUE)` as the transect-level predictor

### “Categorical % coverage” predictors

There are two common ways to do categorical predictors:

#### Recommended for this simplified pipeline: **binary 0/1 rasters per class**

Example:

- `lowland_forest.tif` is a raster where each pixel is:
  - `1` if the pixel is lowland forest
  - `0` otherwise

Then:

- `mean(lowland_forest values in buffer)` = proportion of buffer area covered by lowland forest

This yields a number in `[0, 1]` (multiply by 100 if you want percent).

This is why `lowland_forest`, `lower_montane_forest`, etc. work well as predictors in this simplified pipeline.

#### Optional future improvement: single coded categorical raster

If you have a single raster with integer codes (1..N), then you must:

- count pixels per class within the buffer
- create one column per class (expansion)
- divide each class count by total pixels (proportions)

This requires extra custom code (a small post-processing step) and is intentionally avoided in the current simplified workflow. If you later decide to support coded categorical rasters, you can add a step that expands a coded raster into multiple proportion predictors.

---

## What you must adapt (common edit points)

### 1) Paths (placeholders)

In `best_model/*.R` scripts, set:

- `INPUT_DIR`, `OUTPUT_DIR`, `MODEL_DIR`

They are intentionally placeholders like `"/path/to/your/input"`.

### 2) Predictor list

Edit only:

- `predictor_config.R` → `MODEL_PREDICTORS`

Do not hard-code predictor lists inside `best_model/*.R`. The scripts now source the config.

### 3) Predictor transforms

Edit only:

- `predictor_config.R` → `apply_predictor_transforms()`

This keeps transformations consistent between training and prediction.

### 4) `nest_decay`

You can:

- include a `nest_decay` column in `transects.csv`, **or**
- hard-code it inside `best_model/01_fit_best_model.R` before computing `offset_term`.

---

## Notes / assumptions

- All coordinates are assumed to already be in a planar CRS (e.g. AEA meters) consistent with the aligned rasters.
- Grid ids are generated as `grid_1`, `grid_2`, ... and remain consistent across years.
- The pipeline currently models `nr_nests` with a NB GLM; if you change response definitions you’ll need to update `best_model/01_fit_best_model.R`.
