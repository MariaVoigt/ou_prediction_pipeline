suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))

INPUT_DIR <- "/path/to/your/input"
MODEL_DIR <- "/path/to/your/output"
OUTPUT_DIR <- "/path/to/your/output"

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
script_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
repo_dir <- dirname(script_dir)

source(file.path(repo_dir, "predictor_config.R"))
source(file.path(repo_dir, "helpers", "scale.predictors.R"))

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

predictors_obs <- read.csv(file.path(MODEL_DIR, "predictors_observation_scaled.csv"), stringsAsFactors = FALSE)
best_fit <- readRDS(file.path(MODEL_DIR, "best_glm_nb_model.rds"))

predictor_names_for_scaling <- PREDICTORS_FOR_SCALING
predictor_names_add <- PREDICTORS_ADD

# Set this to the years you want to predict
YEARS_TO_PREDICT <- sort(unique(predictors_obs$unscaled_year))

for (year_to_predict in YEARS_TO_PREDICT) {
  geography_path <- file.path(INPUT_DIR, paste0("geography_", year_to_predict, ".csv"))
  predictors_path <- file.path(INPUT_DIR, paste0("predictors_abundance_", year_to_predict, ".csv"))

  geography <- read.csv(geography_path, stringsAsFactors = FALSE)
  predictors <- read.csv(predictors_path, stringsAsFactors = FALSE)

  geography$unscaled_x_center <- rowMeans(cbind(geography$x_start, geography$x_end), na.rm = TRUE)
  geography$unscaled_y_center <- rowMeans(cbind(geography$y_start, geography$y_end), na.rm = TRUE)

  predictors <- dplyr::rename(predictors, unscaled_value = value)
  predictors <- apply_predictor_transforms(predictors)

  predictors_grid <- scale.predictors.grid(
    predictor_names_for_scaling = predictor_names_for_scaling,
    predictor_names_add = predictor_names_add,
    predictors = predictors,
    predictors_obs = predictors_obs,
    geography = geography
  )

  # Prediction: offset is not used for grid prediction
  predictors_grid_pred <- predictors_grid
  predictors_grid_pred$offset_term <- 0

  pred_nests <- predict(best_fit, newdata = predictors_grid_pred, type = "response")

  pred_per_cell <- data.frame(
    id = predictors_grid$id,
    year = year_to_predict,
    abundance_pred = as.numeric(pred_nests)
  )

  write.csv(predictors_grid, file.path(OUTPUT_DIR, paste0("predictors_grid_scaled_", year_to_predict, ".csv")), row.names = FALSE)
  write.csv(pred_per_cell, file.path(OUTPUT_DIR, paste0("abundance_pred_per_cell_", year_to_predict, ".csv")), row.names = FALSE)
}
