suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(MASS))

MODEL_DIR <- "/path/to/your/output"
OUTPUT_DIR <- "/path/to/your/output"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

predictors_obs <- read.csv(file.path(MODEL_DIR, "predictors_observation_scaled.csv"), stringsAsFactors = FALSE)
best_rhs_tbl <- read.csv(file.path(MODEL_DIR, "best_model_rhs.csv"), stringsAsFactors = FALSE)
best_rhs <- best_rhs_tbl$best_model_rhs[1]

fml <- as.formula(paste0("nr_nests ~ ", best_rhs, " + offset(offset_term)"))

years <- sort(unique(predictors_obs$unscaled_year))

cv_res <- lapply(years, function(y) {
  train <- predictors_obs[predictors_obs$unscaled_year != y, , drop = FALSE]
  test <- predictors_obs[predictors_obs$unscaled_year == y, , drop = FALSE]

  fit <- glm.nb(fml, data = train, control = glm.control(maxit = 500))

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
