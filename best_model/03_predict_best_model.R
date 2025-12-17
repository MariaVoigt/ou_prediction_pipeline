suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))

INPUT_DIR <- "/path/to/your/input"
MODEL_DIR <- "/path/to/your/output"
OUTPUT_DIR <- "/path/to/your/output"

FUN_DIR <- file.path(getwd(), "src", "functions")
source(file.path(FUN_DIR, "project_functions", "scale.predictors.R"))

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

predictors_obs <- read.csv(file.path(MODEL_DIR, "predictors_observation_scaled.csv"), stringsAsFactors = FALSE)
best_fit <- readRDS(file.path(MODEL_DIR, "best_glm_nb_model.rds"))

predictor_names_for_scaling <- c(
  "dem",
  "slope",
  "temp_mean",
  "rain_dry",
  "rain_var",
  "ou_killings",
  "ou_killing_prediction",
  "human_pop_dens",
  "perc_muslim",
  "peatswamp",
  "lowland_forest",
  "lower_montane_forest",
  "road_dens",
  "distance_PA",
  "fire_dens",
  "deforestation_hansen",
  "deforestation_gaveau",
  "plantation_distance",
  "pulp_distance",
  "palm_distance",
  "dom_T_OC",
  "dom_T_PH"
)

predictor_names_add <- c("year", "x_center", "y_center")

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
