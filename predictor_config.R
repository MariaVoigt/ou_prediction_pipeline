# Shared configuration for lightweight model-averaging pipeline

# Predictor set to be used throughout (single source of truth)
MODEL_PREDICTORS <- c(
  "year",
  "temp_mean",
  "rain_var",
  "rain_dry",
  "dom_T_OC",
  "peatswamp",
  "lowland_forest",
  "lower_montane_forest",
  "deforestation_hansen",
  "human_pop_dens"
)

# Derived terms (still based only on predictors above)
INCLUDE_RAIN_DRY_QUADRATIC <- TRUE

# Predictors that are scaled from the long predictor table
# (exclude year/x/y since those come from geography)
PREDICTORS_FOR_SCALING <- setdiff(MODEL_PREDICTORS, c("year", "x_center", "y_center"))

# Additional predictors computed from geography and scaled
PREDICTORS_ADD <- c("year", "x_center", "y_center")

# Transformation rules applied before scaling.
# Implemented as a function to keep behavior explicit.
apply_predictor_transforms <- function(predictors_df) {
  if (!all(c("predictor", "value") %in% names(predictors_df))) {
    stop("predictors_df must contain columns: predictor, value")
  }

  # Only transform predictors we actually use.
  if ("human_pop_dens" %in% unique(predictors_df$predictor)) {
    predictors_df[predictors_df$predictor == "human_pop_dens", "value"] <-
      log(predictors_df[predictors_df$predictor == "human_pop_dens", "value"] + 1)
  }

  return(predictors_df)
}

# Model term universe (m_terms) used for candidate model generation
get_m_terms <- function() {
  base <- MODEL_PREDICTORS
  # built.all.models expects "1" for intercept inclusion
  out <- c("1", base)
  if (isTRUE(INCLUDE_RAIN_DRY_QUADRATIC)) {
    out <- c(out, "I(rain_dry^2)")
  }
  return(out)
}

# Candidate model config (for built.all.models)
get_candidate_model_config <- function() {
  env_cov_names <- MODEL_PREDICTORS
  env_cov_2 <- character(0)
  if (isTRUE(INCLUDE_RAIN_DRY_QUADRATIC)) {
    env_cov_2 <- c("rain_dry")
  }
  list(
    env.cov.names = env_cov_names,
    env.cov.int = list(),
    env.cov.2 = env_cov_2
  )
}
