suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(reshape2))

# One-off model selection script.
# Run this to determine the best model RHS by AIC (all-subsets search),
# then copy the result into BEST_MODEL_RHS in 01_fit_best_model.R.

INPUT_DIR  <- "/path/to/your/input"
OUTPUT_DIR <- "/path/to/your/output"

# Buffer (km) used when preparing predictors; must match prepare step --buffer-km
BUFFER_KM <- 20

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
script_dir  <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
repo_dir    <- dirname(script_dir)

source(file.path(repo_dir, "predictor_config.R"))
source(file.path(repo_dir, "helpers", "scale.predictors.R"))
source(file.path(repo_dir, "helpers", "rogers_model_functions.R"))

ESW <- 0.01595
NCS <- 1.12
PNB <- 0.88

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

geography  <- read.csv(file.path(INPUT_DIR, "geography_observation.csv"), stringsAsFactors = FALSE)
transects  <- read.csv(file.path(INPUT_DIR, "transects.csv"), stringsAsFactors = FALSE)
predictors <- read.csv(file.path(INPUT_DIR, paste0("predictors_observation_", BUFFER_KM, ".csv")), stringsAsFactors = FALSE)

geography$unscaled_x_center <- rowMeans(cbind(geography$x_start, geography$x_end), na.rm = TRUE)
geography$unscaled_y_center <- rowMeans(cbind(geography$y_start, geography$y_end), na.rm = TRUE)

predictors <- apply_predictor_transforms(predictors)

predictors <- dplyr::select(predictors, id, predictor, unscaled_year = year, unscaled_value = value) %>%
  inner_join(transects, by = "id")

predictors <- predictors[!is.na(predictors$unscaled_value), ]

predictors_obs <- scale.predictors.observation(
  predictor_names_for_scaling = PREDICTORS_FOR_SCALING,
  predictor_names_add         = PREDICTORS_ADD,
  predictors                  = predictors,
  geography                   = geography
)

predictors_obs <- geography %>%
  dplyr::select(-c(unscaled_x_center, unscaled_y_center)) %>%
  inner_join(transects, by = "id") %>%
  inner_join(predictors_obs, by = "id") %>%
  arrange(id) %>%
  as.data.frame()

predictors_obs$ou_dens <- (predictors_obs$nr_nests / (predictors_obs$length_km * ESW * 2)) *
  (1 / (predictors_obs$nest_decay * NCS * PNB))

predictors_obs$offset_term <- log(predictors_obs$length_km * ESW * 2 * predictors_obs$nest_decay * NCS * PNB)

# Coordinates in km for the isotropic spatial smooth (see predictor_config.R)
predictors_obs$x_km <- predictors_obs$unscaled_x_center / 1000
predictors_obs$y_km <- predictors_obs$unscaled_y_center / 1000

m_terms <- get_m_terms()

all_model_terms <- built.all.models(
  env.cov.names = get_candidate_model_config()$env.cov.names,
  env.cov.int   = get_candidate_model_config()$env.cov.int,
  env.cov.2     = get_candidate_model_config()$env.cov.2
)

message("Fitting ", nrow(all_model_terms), " candidate models ...")

results <- vector("list", nrow(all_model_terms))

for (i in seq_len(nrow(all_model_terms))) {
  rhs <- paste(m_terms[all_model_terms[i, ] == 1], collapse = "+")
  # The spatial smooth is a fixed term in every candidate (not part of the search).
  fml <- as.formula(paste0("nr_nests ~ ", rhs, " + ", SPATIAL_SMOOTH_TERM, " + offset(offset_term)"))
  fit <- try(gam(fml, data = predictors_obs, family = nb(), method = "ML"), silent = TRUE)
  results[[i]] <- list(
    model_rhs = rhs,
    aic       = if (inherits(fit, "try-error")) Inf else AIC(fit)
  )
}

aics    <- vapply(results, function(x) x$aic, numeric(1))
best_i  <- which.min(aics)
best_rhs <- results[[best_i]]$model_rhs

message("Best model RHS (AIC = ", round(aics[best_i], 2), "):")
message("  ", best_rhs)
message("\nCopy this into BEST_MODEL_RHS in 01_fit_best_model.R")

aic_tbl <- data.frame(
  model_rhs = vapply(results, function(x) x$model_rhs, character(1)),
  AIC       = aics
) %>% dplyr::arrange(AIC)

write.csv(aic_tbl,
          file      = file.path(OUTPUT_DIR, "all_models_aic.csv"),
          row.names = FALSE)

best_summary <- data.frame(best_model_rhs = best_rhs, AIC = aics[best_i])
write.csv(best_summary,
          file      = file.path(OUTPUT_DIR, "best_model_rhs.csv"),
          row.names = FALSE)
