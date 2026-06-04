# Orangutan prediction pipeline (simplified)

This repository contains a lightweight workflow to:

- Prepare predictor tables from **already-aligned rasters** (alignment/clipping can be done in ArcGIS/QGIS).
- Fit a **negative-binomial GAM** abundance model with a **2D spatial smooth** on transect observations.
- Predict (spatially interpolate) abundance on a raster-derived grid for one or more years.

The pipeline works with the best model only, but we can expand to full model comparison.
The logic is: 

- **Rasters (aligned) → extracted predictors in a buffer → long tables → scaled wide tables → NB GAM (env terms + `s(x, y)`) → predictions**

---

## Two analysis tracks

The pipeline supports **two distinct questions** that must not share a model. See [`README_two_tracks.md`](README_two_tracks.md) for full design rationale.

| | **Track A — spatial interpolation** | **Track B — degradation / re-occupation** |
|---|---|---|
| Question | What is orangutan abundance across the landscape in a given year? | Does logging drive a temporary abundance change, and do orangutans return? |
| Unit | grid cell | transect **segment** × year |
| Model | NB GAM with `s(x_km, y_km)` spatial smooth | NB GAM, recovery smooth `s(years_since_last_logging, by = ever_logged)`, transect + year random effects |
| Role of `s(x,y)` | **central** — interpolates across unsurveyed space | **absent** (or nuisance mop-up if residuals cluster spatially) |
| Output | wall-to-wall density surface | degradation contrast + recovery trajectory in occupied habitat |
| Scripts | `prepare_predictors_*`, `best_model/` | `prepare_segments_track_b.R`, `fit_track_b.R` |

> **Why they must stay separate:** Track A's `s(x_km, y_km)` smooth would absorb the spatially structured logging signal that Track B is designed to measure. Track B's transect random effect and recovery smooth answer a change-over-time question that Track A's per-year surface cannot.

### Model form

The model is a negative-binomial GAM:

```
nr_nests ~ <environmental predictors> + s(x_km, y_km) + offset(offset_term)
```

- Environmental predictors enter as **parametric (linear) terms** — these are the interpretable coefficients.
- `s(x_km, y_km)` is an isotropic thin-plate spline over geographic coordinates (in km): this is the **spatial-interpolation** component.

> ## ⚠️ Key trade-off: spatial confounding
>
> The spatial smooth `s(x_km, y_km)` and any **spatially structured environmental predictor** compete for the same variance. Because almost every environmental layer (deforestation, distances, human population, etc.) is spatially patterned, the smooth will **absorb part of their signal**.
>
> Consequences you must report:
> - Environmental coefficients are typically **smaller** (shrunk toward zero) and have **wider standard errors** than in a non-spatial model.
> - The same predictor can look important without the smooth and unimportant with it — neither is "wrong"; they answer different questions.
> - The smooth itself is **not** a single coefficient: interpret it via its EDF and a plot of the fitted surface, not a number in the coefficients table.
>
> **Bottom line:** the GAM is excellent for *spatial prediction/interpolation*, but environmental coefficients must be interpreted *conditional on the spatial smooth*. If clean coefficient inference is the primary goal, a spatial GLMM (Matérn random field) would be the more appropriate tool.

---

## Repository structure

- `predictor_config.R`
  - Single source of truth for:
    - which predictors are used everywhere (`MODEL_PREDICTORS`)
    - transforms applied before scaling (`apply_predictor_transforms()`)
    - candidate model term universe (`get_m_terms()`, `get_candidate_model_config()`)
    - the spatial smooth term (`SPATIAL_SMOOTH_TERM`, `SPATIAL_COORDS`)

- `prepare_predictors_from_aligned_rasters.R`
  - Main *data preparation* entrypoint for **Track A**.
  - Builds the prediction grid from a reference raster and extracts predictors for:
    - transects (buffered extraction)
    - grid cells (cell-level extraction)

- `prepare_segments_track_b.R`
  - **Track B** data preparation.
  - Splits each transect into fixed-length segments (≈250 m), assigns nest GPS points to segments, extracts per-segment buffered covariates including degradation timing from a `year_last_logged` raster, and builds the offset term.
  - Output: `segment_year_table.csv` — one row per segment × year, ready for `fit_track_b.R`.

- `fit_track_b.R`
  - **Track B** model fitting + diagnostics.
  - Fits the NB GAM with the by-factor recovery smooth, runs `gam.check`, tests residual spatial autocorrelation (Moran's I), and performs two sensitivity checks (regression-to-the-mean; `distance_road` as mediator vs. confounder).

- `helpers/scale.predictors.R`
  - Functions to scale predictors and cast long→wide for modelling/prediction.
- `helpers/rogers_model_functions.R`
  - `built.all.models()`: enumerates the all-subsets candidate model matrix (used by `00`).

- `best_model/`
  - `00_select_best_model.R`: **(one-off)** all-subsets search (each candidate = env subset **+ fixed `s(x_km, y_km)`**) by AIC; prints best RHS and writes `all_models_aic.csv` + `best_model_rhs.csv`
  - `01_fit_best_model.R`: fit the pre-specified best NB GAM (set `BEST_MODEL_RHS` at the top); writes `best_gam_nb_model.rds`, parametric `best_model_coefficients.csv`, and `best_model_smooth_terms.csv`
  - `02_cross_validation_leave_one_year_out.R`: leave-one-year-out CV — **not suitable for the spatial use case** (see caveat in the file); for temporal diagnostics only
  - `02b_spatial_cross_validation.R`: **spatial block CV** — groups transects into spatial blocks/folds and measures interpolation skill (R², RMSE, Spearman). Use this to validate spatial prediction.
  - `03_predict_best_model.R`: spatially interpolate on grid for multiple years, with a configurable prediction footprint
  - `04_bootstrap_best_model.R`: coefficient uncertainty via GAM posterior simulation (parametric terms)

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
| ground_12 | 2010 | slope | xx |
| ground_12 | 2016 | temp_mean | 26.3 |
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
- gps location of nest
- perpendicular distance to nest
- think about how to code 0 nest transects (to do!)

- think how we aggregate, whether we deal with individual nests as a response or transects segment

`nest_decay` can be provided as a column or set/hard-coded inside modelling scripts.
- 602 days

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

- add script for predictor interrogation for distribution and correlation analysis


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

> **Tip — keep `MODEL_PREDICTORS` small before running model selection.**
> All-subsets search in `00_select_best_model.R` fits **2^N candidate models**, where N is the number of predictors.
> With 16 predictors that is ~65,000 models; with 10 it is ~1,000.
> Start with a biologically motivated shortlist, run selection, then expand only if needed.

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

### Step 3a — Select the best model (one-off)

If you do not yet know your best model formula, run `best_model/00_select_best_model.R`:

- Set `INPUT_DIR` and `OUTPUT_DIR` at the top.
- The script fits all candidate NB **GAMs** (all-subsets over `MODEL_PREDICTORS`, each with the **fixed spatial smooth `s(x_km, y_km)`**) and selects the lowest-AIC one.
- Only the **environmental** terms are searched; the spatial smooth is always present.
- The best RHS (environmental terms only) is printed to the console and written to `best_model_rhs.csv`.
- A full AIC table is written to `all_models_aic.csv`.

**This step is slow** (up to 2^16 candidate GAM fits). Run it once, then record the result.

### Step 3b — Fit the best model

Edit `best_model/01_fit_best_model.R`:

- Set:
  - `INPUT_DIR` to the folder created in Step 2
  - `OUTPUT_DIR` to a model output folder
  - `BEST_MODEL_RHS` to the **environmental** formula RHS from Step 3a (or from prior knowledge). The spatial smooth and offset are appended automatically — do **not** include them here.

Then run the script. It will write (among others):

- `best_gam_nb_model.rds` (the fitted NB GAM)
- `predictors_observation_scaled.csv`
- `best_model_rhs.csv`
- `best_model_coefficients.csv` (parametric/environmental terms with SEs)
- `best_model_smooth_terms.csv` (EDF/significance of the spatial smooth)

---

### Step 4 — Predict on the grid

Edit `best_model/03_predict_best_model.R`:

- Set:
  - `INPUT_DIR` to the folder created in Step 2
  - `MODEL_DIR` to the folder created in Step 3
  - `OUTPUT_DIR` to a prediction output folder

Optionally edit:

- `YEARS_TO_PREDICT`

#### Prediction footprint (extrapolation handling)

Because the spatial smooth extrapolates poorly outside the survey area, you choose a **footprint** that defines where predictions are trustworthy:

- **`FOOTPRINT_MODE`**
  - `"none"` — predict everywhere
  - `"buffer"` — keep cells within `FOOTPRINT_BUFFER_KM` of any transect center
  - `"convex_hull"` — convex hull of transects, optionally expanded by `FOOTPRINT_BUFFER_KM`
  - `"shapefile"` — a polygon you provide via `FOOTPRINT_SHP` (must be in the **same planar CRS** as the coordinates)
- **`FOOTPRINT_ACTION`**
  - `"flag"` — keep all cells but add an `in_footprint` column to the output
  - `"clip"` — drop out-of-footprint cells entirely

Then run. It writes per year:

- `predictors_grid_scaled_<year>.csv`
- `abundance_pred_per_cell_<year>.csv` — columns: `density_ou_per_km2`, `cell_area_km2`, `abundance_pred` (absolute orangutans in the cell), `in_footprint`
- `abundance_total_<year>.csv` — landscape population total summed over in-footprint cells

**Units:** with `offset = 0` the model predicts orangutan **density** (individuals/km²). The script multiplies by each cell's area (from the cell bounds, assuming coordinates in meters) to get **absolute abundance per cell**, which sums to the landscape total.

### Step 5 — Validate spatially

Run `best_model/02b_spatial_cross_validation.R` to assess interpolation skill:

- Set `MODEL_DIR` / `OUTPUT_DIR`, and optionally `BLOCK_SIZE_KM` and `N_FOLDS`.
- Transects are grouped into spatial blocks, blocks assigned to folds; each fold is held out and predicted.
- Writes `cv_spatial_blocks.csv` (per-fold R², RMSE, Spearman) and `cv_spatial_summary.csv`.

> Prefer this over `02_cross_validation_leave_one_year_out.R` for the spatial use case. Year-based CV is only a temporal diagnostic and is degenerate when you have few survey years.

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

### 4) `BEST_MODEL_RHS`

Set this in `best_model/01_fit_best_model.R` to the **environmental** terms of your chosen model, e.g.:

```r
BEST_MODEL_RHS <- "dem + slope + lowland_forest + human_pop_dens + distance_PA"
```

Do not include `offset(offset_term)` **or** the spatial smooth `s(x_km, y_km)` here — both are appended automatically (the smooth is configured via `SPATIAL_SMOOTH_TERM` in `predictor_config.R`).

### 5) `nest_decay`

You can:

- include a `nest_decay` column in `transects.csv`, **or**
- hard-code it inside `best_model/01_fit_best_model.R` before computing `offset_term`.

---

## Track B workflow

### Step B1 — Prepare segment-year table

Open `prepare_segments_track_b.R` and edit the **USER SETTINGS (EDIT!)** block:

- `transects_shp` — LINESTRING shapefile, one row per transect × year, must carry a transect id column and a year column.
- `nests_shp` — POINT shapefile, one row per nest, must carry a year column and a transect id column.
- `raster_dir` — same aligned rasters as Track A.
- `logging_year_rast` — raster where each pixel value is the most recent year that pixel was logged (`NA` = never logged).
- `out_dir` — output folder.
- `env_covariates` — named vector of covariate rasters to extract (e.g. `c("dem", "slope", "distance_road")`). For each name `p` the script tries `p_<year>.tif` first then falls back to `p.tif`.
- `segment_length_m` — recommended 200–250 m.
- `buffer_m` — **keep ≤ `segment_length_m`**; a wider buffer causes adjacent segment buffers to overlap, manufacturing spatial autocorrelation and eroding the within-transect contrast.
- `nest_decay_days` — mean nest decay time in days. **Scrutinise this parameter first** — if decay differs between degraded and intact forest, a constant decay term in the offset can manufacture the very trend Track B tests for.

Then **Source** the script (or run via `Rscript prepare_segments_track_b.R --help` for CLI options).

#### What this produces

In `out_dir`:

- `segment_year_table.csv` — one row per segment × year with columns: `seg_id`, `transect_id`, `year`, `seg_index`, `length_km`, `nr_nests`, `prop_degraded_1yr_prior`, `years_since_last_logging`, `ever_logged`, environmental covariates, `offset_term`, `x_km`, `y_km`, `baseline_nests_c`.
- `segments_geometry.gpkg` — segment linestrings for mapping/QC.

Diagnostics printed to console: zero-nest fraction, buffer/segment ratio warning, logged vs. never-logged counts.

### Step B2 — Fit the recovery model

Edit `fit_track_b.R`:

- `INPUT_CSV` — path to `segment_year_table.csv` from Step B1.
- `OUTPUT_DIR` — model output folder.
- `RECOVERY_K` — basis dimension for the recovery smooth (default 4; raise only if `gam.check` flags `k` too low).
- `USE_DISTANCE_ROAD` — include `distance_road` as a covariate (default `TRUE`; sensitivity check refits without it).

The script auto-detects whether there are enough logged segment-years for a spline (≥ 30) or falls back to a quadratic recovery term.

#### What this produces

In `OUTPUT_DIR`:

- `track_b_gam.rds` — fitted NB GAM.
- `model_summary.txt`, `parametric_coefficients.csv`, `smooth_terms.csv`.
- `recovery_curve.png` — **headline plot**: deviation from never-logged baseline vs. years since last logging. A dip-then-return pattern = re-occupation signal.
- `gam_check.png` / `gam_check.txt` — residual diagnostics.
- `residuals_in_space.png` + `spatial_autocorr.txt` — Moran's I on between-transect deviance residuals. If strongly positive, refit adding `s(x_km, y_km)` as a nuisance term.
- `sensitivity_baseline.csv` — key terms with vs. without `baseline_nests_c`.
- `sensitivity_road.csv` — key terms with vs. without `distance_road`.

### Track B interpretation reminders

- **Defensible result**: the degradation/recovery *contrast* within the revisited set (revisit selection was independent of disturbance).
- **Caveat — absolute trajectory**: transects were selected on high baseline occupancy, so regression to the mean can mimic a real decline. The `baseline_nests_c` sensitivity check quantifies this; report the with-baseline version.
- **Caveat — recolonisation**: transects empty at baseline were dropped. Recovery is censored at the baseline-occupied level — we can see a segment dip and return toward its own baseline, but not orangutans appearing in previously empty forest.

---

## Notes / assumptions

- All coordinates are assumed to already be in a planar CRS (e.g. AEA meters) consistent with the aligned rasters.
- Grid ids are generated as `grid_1`, `grid_2`, ... and remain consistent across years.
- The pipeline models `nr_nests` with a **negative-binomial GAM** (`mgcv::gam`, `family = nb()`) including an isotropic spatial smooth `s(x_km, y_km)`; if you change response definitions you’ll need to update `best_model/01_fit_best_model.R`.
- The spatial smooth uses coordinates in **km** (built from the unscaled cell/transect centers) so x and y share a common isotropic scale. Predictions for grid cells far outside the transect footprint are **extrapolations** and unreliable — use the prediction footprint (Step 4) to flag/clip them.
- **Validate spatially**, not temporally: use `02b_spatial_cross_validation.R` (spatial block CV). The leave-one-year-out script is for temporal diagnostics only.
- Required packages: `mgcv` (model), and `sp` + `rgeos` + `rgdal` (prediction footprint, same stack as the prepare script).
- Per-cell output is **absolute abundance**: the density prediction (individuals/km², from `offset = 0`) is multiplied by each cell's area (km², derived from the cell bounds). Cells sum to a landscape total (`abundance_total_<year>.csv`). This assumes coordinates are in **meters**.
