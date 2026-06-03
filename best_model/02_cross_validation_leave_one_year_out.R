suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(mgcv))

# CAVEAT: leave-one-year-out CV is NOT appropriate for the spatial-interpolation
# use case (especially with only a few survey years -- e.g. 2 years would train
# on a single year). To evaluate spatial prediction skill, replace this with
# spatial block CV (hold out geographic blocks/transect clusters). This script
# is kept only for temporal diagnostics when several years are available.

MODEL_DIR <- "/path/to/your/output"
OUTPUT_DIR <- "/path/to/your/output"

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
script_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
repo_dir <- dirname(script_dir)
source(file.path(repo_dir, "predictor_config.R"))

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

predictors_obs <- read.csv(file.path(MODEL_DIR, "predictors_observation_scaled.csv"), stringsAsFactors = FALSE)
best_rhs_tbl <- read.csv(file.path(MODEL_DIR, "best_model_rhs.csv"), stringsAsFactors = FALSE)
best_rhs <- best_rhs_tbl$best_model_rhs[1]

fml <- as.formula(paste0("nr_nests ~ ", best_rhs, " + ", SPATIAL_SMOOTH_TERM, " + offset(offset_term)"))

years <- sort(unique(predictors_obs$unscaled_year))

cv_res <- lapply(years, function(y) {
  train <- predictors_obs[predictors_obs$unscaled_year != y, , drop = FALSE]
  test <- predictors_obs[predictors_obs$unscaled_year == y, , drop = FALSE]

  fit <- gam(fml, data = train, family = nb(), method = "REML")

  test_pred_data <- test
  test_pred_data$offset_term <- 0

  pred_nests <- predict(fit, newdata = test_pred_data, type = "response")

  cross_lm <- lm(log(test$ou_dens + 1) ~ log(pred_nests + 1))

  data.frame(
    excluded_year = y,
    n_test = nrow(test),
    R2_cross = summary(cross_lm)$r.squared
  )
})

cv_res <- bind_rows(cv_res)

write.csv(cv_res, file = file.path(OUTPUT_DIR, "cv_leave_one_year_out.csv"), row.names = FALSE)
