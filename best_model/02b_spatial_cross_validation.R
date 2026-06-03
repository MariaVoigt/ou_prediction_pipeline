suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(mgcv))

# Spatial block cross-validation for the NB GAM.
# Transects are grouped into square spatial blocks (BLOCK_SIZE_KM); blocks are
# assigned to N_FOLDS folds. Each fold is held out in turn, the model is refit
# on the remaining folds, and prediction skill is measured on the held-out
# block(s). This honestly assesses *spatial* prediction (interpolation) skill,
# unlike leave-one-year-out CV.

MODEL_DIR <- "/path/to/your/output"
OUTPUT_DIR <- "/path/to/your/output"

BLOCK_SIZE_KM <- 50   # size of square spatial blocks
N_FOLDS       <- 5    # number of spatial folds
SET_SEED      <- 1

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
script_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
repo_dir <- dirname(script_dir)
source(file.path(repo_dir, "predictor_config.R"))

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
set.seed(SET_SEED)

predictors_obs <- read.csv(file.path(MODEL_DIR, "predictors_observation_scaled.csv"), stringsAsFactors = FALSE)
best_rhs_tbl <- read.csv(file.path(MODEL_DIR, "best_model_rhs.csv"), stringsAsFactors = FALSE)
best_rhs <- best_rhs_tbl$best_model_rhs[1]

fml <- as.formula(paste0("nr_nests ~ ", best_rhs, " + ", SPATIAL_SMOOTH_TERM, " + offset(offset_term)"))

# Assign each observation to a spatial block, then blocks to folds.
predictors_obs$block <- paste(floor(predictors_obs$x_km / BLOCK_SIZE_KM),
                              floor(predictors_obs$y_km / BLOCK_SIZE_KM),
                              sep = "_")

blocks <- unique(predictors_obs$block)
n_folds <- min(N_FOLDS, length(blocks))
fold_of_block <- setNames(sample(rep(seq_len(n_folds), length.out = length(blocks))), blocks)
predictors_obs$fold <- fold_of_block[predictors_obs$block]

cv_res <- lapply(seq_len(n_folds), function(k) {
  train <- predictors_obs[predictors_obs$fold != k, , drop = FALSE]
  test  <- predictors_obs[predictors_obs$fold == k, , drop = FALSE]

  fit <- tryCatch(gam(fml, data = train, family = nb(), method = "REML"),
                  error = function(e) NULL)
  if (is.null(fit)) {
    return(data.frame(fold = k, n_test = nrow(test), n_blocks = length(unique(test$block)),
                      R2_cross = NA_real_, rmse = NA_real_, spearman = NA_real_))
  }

  test_pred_data <- test
  test_pred_data$offset_term <- 0
  pred_nests <- as.numeric(predict(fit, newdata = test_pred_data, type = "response"))

  cross_lm <- lm(log(test$ou_dens + 1) ~ log(pred_nests + 1))

  data.frame(
    fold     = k,
    n_test   = nrow(test),
    n_blocks = length(unique(test$block)),
    R2_cross = summary(cross_lm)$r.squared,
    rmse     = sqrt(mean((test$ou_dens - pred_nests)^2, na.rm = TRUE)),
    spearman = suppressWarnings(cor(test$ou_dens, pred_nests, method = "spearman", use = "complete.obs"))
  )
})

cv_res <- bind_rows(cv_res)

cv_overall <- data.frame(
  block_size_km = BLOCK_SIZE_KM,
  n_folds       = n_folds,
  n_blocks      = length(blocks),
  mean_R2_cross = mean(cv_res$R2_cross, na.rm = TRUE),
  mean_rmse     = mean(cv_res$rmse, na.rm = TRUE),
  mean_spearman = mean(cv_res$spearman, na.rm = TRUE)
)

write.csv(cv_res, file = file.path(OUTPUT_DIR, "cv_spatial_blocks.csv"), row.names = FALSE)
write.csv(cv_overall, file = file.path(OUTPUT_DIR, "cv_spatial_summary.csv"), row.names = FALSE)

message("Spatial block CV (", n_folds, " folds, ", BLOCK_SIZE_KM, " km blocks):")
message("  mean R2 = ", round(cv_overall$mean_R2_cross, 3),
        " | mean RMSE = ", round(cv_overall$mean_rmse, 3),
        " | mean Spearman = ", round(cv_overall$mean_spearman, 3))
