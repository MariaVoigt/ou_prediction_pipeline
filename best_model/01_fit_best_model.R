suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(reshape2))

INPUT_DIR <- "/path/to/your/input"
OUTPUT_DIR <- "/path/to/your/output"

# Buffer (km) used when preparing predictors; must match prepare step --buffer-km
BUFFER_KM <- 20

# Right-hand side of the best model formula (edit this!)
# Environmental (parametric) terms only. The spatial smooth and offset are
# appended automatically.
BEST_MODEL_RHS <- "PREDICTOR1 + PREDICTOR2 + ..."

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
script_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
repo_dir <- dirname(script_dir)

source(file.path(repo_dir, "predictor_config.R"))
source(file.path(repo_dir, "helpers", "scale.predictors.R"))

ESW <- 0.01595
NCS <- 1.12
PNB <- 0.88

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

geography <- read.csv(file.path(INPUT_DIR, "geography_observation.csv"), stringsAsFactors = FALSE)
transects <- read.csv(file.path(INPUT_DIR, "transects.csv"), stringsAsFactors = FALSE)
predictors <- read.csv(file.path(INPUT_DIR, paste0("predictors_observation_", BUFFER_KM, ".csv")), stringsAsFactors = FALSE)

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

best_rhs <- BEST_MODEL_RHS
fml <- as.formula(paste0("nr_nests ~ ", best_rhs, " + ", SPATIAL_SMOOTH_TERM, " + offset(offset_term)"))
best_fit <- gam(fml, data = predictors_obs, family = nb(), method = "REML")

write.csv(predictors_obs, file = file.path(OUTPUT_DIR, "predictors_observation_scaled.csv"), row.names = FALSE)

best_rhs_tbl <- data.frame(best_model_rhs = best_rhs)
write.csv(best_rhs_tbl, file = file.path(OUTPUT_DIR, "best_model_rhs.csv"), row.names = FALSE)

smry <- summary(best_fit)

# Parametric (environmental) coefficients -- these are the interpretable terms.
# NOTE: estimates are conditional on the spatial smooth; spatially structured
# predictors can be partly absorbed by s(x_km, y_km) (spatial confounding).
p_tbl <- as.data.frame(smry$p.table)
names(p_tbl) <- c("estimate", "std_error", "statistic", "p_value")
p_tbl$term <- rownames(p_tbl)
p_tbl <- p_tbl[, c("term", "estimate", "std_error", "statistic", "p_value")]
write.csv(p_tbl, file = file.path(OUTPUT_DIR, "best_model_coefficients.csv"), row.names = FALSE)

# Smooth-term summary (EDF / significance of the spatial surface)
s_tbl <- as.data.frame(smry$s.table)
names(s_tbl) <- c("edf", "ref_df", "statistic", "p_value")
s_tbl$smooth <- rownames(s_tbl)
s_tbl <- s_tbl[, c("smooth", "edf", "ref_df", "statistic", "p_value")]
write.csv(s_tbl, file = file.path(OUTPUT_DIR, "best_model_smooth_terms.csv"), row.names = FALSE)

saveRDS(best_fit, file = file.path(OUTPUT_DIR, "best_gam_nb_model.rds"))

best_summary <- data.frame(
  model_rhs = best_rhs,
  AIC = AIC(best_fit),
  deviance_explained = smry$dev.expl
)
write.csv(best_summary, file = file.path(OUTPUT_DIR, "best_model_summary.csv"), row.names = FALSE)
