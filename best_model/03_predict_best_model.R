suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(sp))

INPUT_DIR <- "/path/to/your/input"
MODEL_DIR <- "/path/to/your/output"
OUTPUT_DIR <- "/path/to/your/output"

# --- Prediction footprint (extrapolation handling) ---
# The spatial smooth s(x_km, y_km) extrapolates poorly far from transects.
# Restrict or flag grid cells to a footprint where interpolation is trustworthy.
# FOOTPRINT_MODE:
#   "none"        -> predict everywhere (no footprint)
#   "buffer"      -> region within FOOTPRINT_BUFFER_KM of any transect center
#   "convex_hull" -> convex hull of transects, optionally expanded by FOOTPRINT_BUFFER_KM
#   "shapefile"   -> a polygon you provide (must be in the SAME planar CRS as the coordinates)
FOOTPRINT_MODE      <- "buffer"
FOOTPRINT_BUFFER_KM <- 20
FOOTPRINT_SHP       <- ""
# FOOTPRINT_ACTION: "flag" adds an in_footprint column; "clip" drops cells outside the footprint.
FOOTPRINT_ACTION    <- "flag"

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
script_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
repo_dir <- dirname(script_dir)

source(file.path(repo_dir, "predictor_config.R"))
source(file.path(repo_dir, "helpers", "scale.predictors.R"))

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

predictors_obs <- read.csv(file.path(MODEL_DIR, "predictors_observation_scaled.csv"), stringsAsFactors = FALSE)
best_fit <- readRDS(file.path(MODEL_DIR, "best_gam_nb_model.rds"))

predictor_names_for_scaling <- PREDICTORS_FOR_SCALING
predictor_names_add <- PREDICTORS_ADD

# Build the prediction footprint once (uses transect centers from training data).
build_prediction_footprint <- function() {
  if (FOOTPRINT_MODE == "none") return(NULL)

  if (FOOTPRINT_MODE == "shapefile") {
    if (!nzchar(FOOTPRINT_SHP)) stop("FOOTPRINT_SHP must be set when FOOTPRINT_MODE = 'shapefile'")
    poly <- rgdal::readOGR(dsn = dirname(FOOTPRINT_SHP),
                           layer = tools::file_path_sans_ext(basename(FOOTPRINT_SHP)),
                           verbose = FALSE)
    return(as(poly, "SpatialPolygons"))
  }

  tx <- unique(predictors_obs[, c("unscaled_x_center", "unscaled_y_center")])
  tx <- tx[stats::complete.cases(tx), , drop = FALSE]
  pts <- SpatialPoints(as.matrix(tx))

  if (FOOTPRINT_MODE == "buffer") {
    return(rgeos::gBuffer(pts, width = FOOTPRINT_BUFFER_KM * 1000, byid = FALSE))
  }
  if (FOOTPRINT_MODE == "convex_hull") {
    hull <- rgeos::gConvexHull(pts)
    if (FOOTPRINT_BUFFER_KM > 0) hull <- rgeos::gBuffer(hull, width = FOOTPRINT_BUFFER_KM * 1000)
    return(hull)
  }
  stop("Unknown FOOTPRINT_MODE: ", FOOTPRINT_MODE)
}

footprint <- build_prediction_footprint()

# TRUE/FALSE: are these (unscaled) coordinates inside the footprint?
in_footprint <- function(x, y) {
  if (is.null(footprint)) return(rep(TRUE, length(x)))
  grid_pts <- SpatialPoints(cbind(x, y))
  suppressWarnings(proj4string(grid_pts) <- proj4string(footprint))
  !is.na(over(grid_pts, footprint))
}

# Set this to the years you want to predict
YEARS_TO_PREDICT <- sort(unique(predictors_obs$unscaled_year))

for (year_to_predict in YEARS_TO_PREDICT) {
  geography_path <- file.path(INPUT_DIR, paste0("geography_", year_to_predict, ".csv"))
  predictors_path <- file.path(INPUT_DIR, paste0("predictors_abundance_", year_to_predict, ".csv"))

  geography <- read.csv(geography_path, stringsAsFactors = FALSE)
  predictors <- read.csv(predictors_path, stringsAsFactors = FALSE)

  geography$unscaled_x_center <- rowMeans(cbind(geography$x_start, geography$x_end), na.rm = TRUE)
  geography$unscaled_y_center <- rowMeans(cbind(geography$y_start, geography$y_end), na.rm = TRUE)

  # Cell area in km^2 from the cell bounds (coordinates assumed in meters).
  geography$cell_area_km2 <- abs((geography$x_end - geography$x_start) / 1000 *
                                 (geography$y_end - geography$y_start) / 1000)

  predictors <- dplyr::rename(predictors, unscaled_value = value)
  predictors <- apply_predictor_transforms(predictors)

  predictors_grid <- scale.predictors.grid(
    predictor_names_for_scaling = predictor_names_for_scaling,
    predictor_names_add = predictor_names_add,
    predictors = predictors,
    predictors_obs = predictors_obs,
    geography = geography
  )

  # Coordinates in km for the isotropic spatial smooth (must match training)
  predictors_grid$x_km <- predictors_grid$unscaled_x_center / 1000
  predictors_grid$y_km <- predictors_grid$unscaled_y_center / 1000

  # Footprint membership (extrapolation handling)
  predictors_grid$in_footprint <- in_footprint(predictors_grid$unscaled_x_center,
                                                predictors_grid$unscaled_y_center)

  # "clip" drops out-of-footprint cells entirely; "flag" keeps them with a flag.
  if (FOOTPRINT_ACTION == "clip") {
    predictors_grid <- predictors_grid[predictors_grid$in_footprint, , drop = FALSE]
  }

  # Prediction: offset is not used for grid prediction
  predictors_grid_pred <- predictors_grid
  predictors_grid_pred$offset_term <- 0

  # With offset = 0 the prediction is orangutan DENSITY (individuals per km^2).
  density_pred <- as.numeric(predict(best_fit, newdata = predictors_grid_pred, type = "response"))

  # Absolute abundance per cell = density (per km^2) * cell area (km^2).
  cell_area_km2 <- geography$cell_area_km2[match(predictors_grid$id, geography$id)]
  abundance_pred <- density_pred * cell_area_km2

  pred_per_cell <- data.frame(
    id = predictors_grid$id,
    year = year_to_predict,
    density_ou_per_km2 = density_pred,
    cell_area_km2 = cell_area_km2,
    abundance_pred = abundance_pred,
    in_footprint = predictors_grid$in_footprint
  )

  # Landscape total over in-footprint cells (the summable population estimate).
  total_abundance <- sum(pred_per_cell$abundance_pred[pred_per_cell$in_footprint], na.rm = TRUE)
  abundance_total <- data.frame(
    year = year_to_predict,
    total_abundance_in_footprint = total_abundance,
    n_cells_in_footprint = sum(pred_per_cell$in_footprint, na.rm = TRUE)
  )

  write.csv(predictors_grid, file.path(OUTPUT_DIR, paste0("predictors_grid_scaled_", year_to_predict, ".csv")), row.names = FALSE)
  write.csv(pred_per_cell, file.path(OUTPUT_DIR, paste0("abundance_pred_per_cell_", year_to_predict, ".csv")), row.names = FALSE)
  write.csv(abundance_total, file.path(OUTPUT_DIR, paste0("abundance_total_", year_to_predict, ".csv")), row.names = FALSE)
}
