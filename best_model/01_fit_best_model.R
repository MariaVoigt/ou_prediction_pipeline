suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(reshape2))

INPUT_DIR <- "/path/to/your/input"
OUTPUT_DIR <- "/path/to/your/output"

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
script_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
repo_dir <- dirname(script_dir)

source(file.path(repo_dir, "predictor_config.R"))
source(file.path(repo_dir, "helpers", "scale.predictors.R"))
source(file.path(repo_dir, "helpers", "rogers_model_functions.R"))

ESW <- 0.01595
NCS <- 1.12
PNB <- 0.88

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

geography <- read.csv(file.path(INPUT_DIR, "geography_observation.csv"), stringsAsFactors = FALSE)
transects <- read.csv(file.path(INPUT_DIR, "transects.csv"), stringsAsFactors = FALSE)
predictors <- read.csv(file.path(INPUT_DIR, "predictors_observation_20.csv"), stringsAsFactors = FALSE)

geography$unscaled_x_center <- rowMeans(cbind(geography$x_start, geography$x_end), na.rm = TRUE)
geography$unscaled_y_center <- rowMeans(cbind(geography$y_start, geography$y_end), na.rm = TRUE)

predictors <- apply_predictor_transforms(predictors)

predictor_names_for_scaling <- PREDICTORS_FOR_SCALING
predictor_names_add <- PREDICTORS_ADD

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

m_terms <- get_m_terms()

all_model_terms <- built.all.models(
  env.cov.names = get_candidate_model_config()$env.cov.names,
  env.cov.int = get_candidate_model_config()$env.cov.int,
  env.cov.2 = get_candidate_model_config()$env.cov.2
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
