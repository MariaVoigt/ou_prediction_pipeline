suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(reshape2))

INPUT_DIR <- "/path/to/your/input"
OUTPUT_DIR <- "/path/to/your/output"

FUN_DIR <- file.path(getwd(), "src", "functions")
source(file.path(FUN_DIR, "project_functions", "scale.predictors.R"))
source(file.path(FUN_DIR, "roger_functions", "rogers_model_functions.R"))

ESW <- 0.01595
NCS <- 1.12
PNB <- 0.88

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

geography <- read.csv(file.path(INPUT_DIR, "geography_observation.csv"), stringsAsFactors = FALSE)
transects <- read.csv(file.path(INPUT_DIR, "transects.csv"), stringsAsFactors = FALSE)
predictors <- read.csv(file.path(INPUT_DIR, "predictors_observation_20.csv"), stringsAsFactors = FALSE)

geography$unscaled_x_center <- rowMeans(cbind(geography$x_start, geography$x_end), na.rm = TRUE)
geography$unscaled_y_center <- rowMeans(cbind(geography$y_start, geography$y_end), na.rm = TRUE)

predictors[predictors$predictor == "distance_PA", "value"] <- sqrt(predictors[predictors$predictor == "distance_PA", "value"])
predictors[predictors$predictor == "human_pop_dens", "value"] <- log(predictors[predictors$predictor == "human_pop_dens", "value"] + 1)
predictors[predictors$predictor == "deforestation_gaveau", "value"] <- sqrt(predictors[predictors$predictor == "deforestation_gaveau", "value"])
predictors[predictors$predictor == "plantation_distance", "value"] <- log(predictors[predictors$predictor == "plantation_distance", "value"] + 1)
predictors[predictors$predictor == "pulp_distance", "value"] <- log(predictors[predictors$predictor == "pulp_distance", "value"] + 1)
predictors[predictors$predictor == "palm_distance", "value"] <- log(predictors[predictors$predictor == "palm_distance", "value"] + 1)

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

predictors <- dplyr::select(predictors, id, predictor, unscaled_year = year, unscaled_value = value) %>%
  inner_join(transects, by = "id")

predictors <- predictors[!is.na(predictors$unscaled_value), ]

predictors_obs <- scale.predictors.observation(
  predictor_names_for_scaling = predictor_names_for_scaling,
  predictor_names_add = predictor_names_add,
  predictors = predictors,
  geography = geography
)

predictors_obs <- geography %>%
  dplyr::select(-c(year, unscaled_x_center, unscaled_y_center)) %>%
  inner_join(transects, by = "id") %>%
  inner_join(predictors_obs, by = "id") %>%
  arrange(id) %>%
  as.data.frame()

predictors_obs$ou_dens <- (predictors_obs$nr_nests / (predictors_obs$length_km * ESW * 2)) *
  (1 / (predictors_obs$nest_decay * NCS * PNB))

predictors_obs$offset_term <- log(predictors_obs$length_km * ESW * 2 * predictors_obs$nest_decay * NCS * PNB)

m_terms <- c(
  "1",
  "year",
  "temp_mean",
  "rain_var",
  "rain_dry",
  "dom_T_OC",
  "peatswamp",
  "lowland_forest",
  "lower_montane_forest",
  "deforestation_hansen",
  "human_pop_dens",
  "ou_killing_prediction",
  "perc_muslim",
  "I(rain_dry^2)"
)

all_model_terms <- built.all.models(
  env.cov.names = c(
    "year",
    "temp_mean",
    "rain_var",
    "rain_dry",
    "dom_T_OC",
    "peatswamp",
    "lowland_forest",
    "lower_montane_forest",
    "deforestation_hansen",
    "human_pop_dens",
    "ou_killing_prediction",
    "perc_muslim"
  ),
  env.cov.int = list(),
  env.cov.2 = c("rain_dry")
)

results <- vector("list", nrow(all_model_terms))

for (i in seq_len(nrow(all_model_terms))) {
  rhs <- paste(m_terms[all_model_terms[i, ] == 1], collapse = "+")
  fml <- as.formula(paste0("nr_nests ~ ", rhs, " + offset(offset_term)"))
  fit <- glm.nb(fml, data = predictors_obs, control = glm.control(maxit = 250))
  results[[i]] <- list(
    model_rhs = rhs,
    aic = extractAIC(fit)[2],
    fit = fit
  )
}

aics <- vapply(results, function(x) x$aic, numeric(1))
best_i <- which.min(aics)
best_fit <- results[[best_i]]$fit
best_rhs <- results[[best_i]]$model_rhs

write.csv(predictors_obs, file = file.path(OUTPUT_DIR, "predictors_observation_scaled.csv"), row.names = FALSE)

best_rhs_tbl <- data.frame(best_model_rhs = best_rhs)
write.csv(best_rhs_tbl, file = file.path(OUTPUT_DIR, "best_model_rhs.csv"), row.names = FALSE)

best_coef_tbl <- data.frame(term = names(coef(best_fit)), estimate = as.numeric(coef(best_fit)))
write.csv(best_coef_tbl, file = file.path(OUTPUT_DIR, "best_model_coefficients.csv"), row.names = FALSE)

saveRDS(best_fit, file = file.path(OUTPUT_DIR, "best_glm_nb_model.rds"))

best_summary <- data.frame(
  model_rhs = best_rhs,
  AIC = extractAIC(best_fit)[2]
)
write.csv(best_summary, file = file.path(OUTPUT_DIR, "best_model_summary.csv"), row.names = FALSE)
