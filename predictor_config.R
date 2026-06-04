# Shared configuration for lightweight model-averaging pipeline

# Predictor set to be used throughout (single source of truth)
MODEL_PREDICTORS <- c(
  "dem",
  "slope",
  "road_distance",
  "proportion_degraded_forest_1yrs"
  "years_since_newest_change" # 10 years - value ten, last year it would be 1, two years ago 2
)



# newly degraded forest int he last two years (to do- create layer)

# Predictors that are scaled from the long predictor table
# (exclude year/x/y since those come from geography)
PREDICTORS_FOR_SCALING <- setdiff(MODEL_PREDICTORS, c("year", "x_center", "y_center"))

# Additional predictors computed from geography and scaled
PREDICTORS_ADD <- c("year", "x_center", "y_center")

# --- Spatial model (GAM 2D smooth) ---
# The pipeline fits a negative-binomial GAM with an isotropic thin-plate
# spline over geographic coordinates: this is the spatial-interpolation term.
# Coordinates are used in KILOMETERS (built in the scripts from the unscaled
# center coordinates) so that x and y share a common, isotropic scale.
# NOTE: do NOT z-score x/y separately for the smooth -- that distorts distances.
SPATIAL_COORDS <- c("x_km", "y_km")
SPATIAL_SMOOTH_TERM <- 's(x_km, y_km, bs = "tp")'

# Transformation rules applied before scaling.
# Implemented as a function to keep behavior explicit.
apply_predictor_transforms <- function(predictors_df) {
  if (!all(c("predictor", "value") %in% names(predictors_df))) {
    stop("predictors_df must contain columns: predictor, value")
  }

  if ("human_pop_dens" %in% unique(predictors_df$predictor)) {
    predictors_df[predictors_df$predictor == "human_pop_dens", "value"] <-
      log(predictors_df[predictors_df$predictor == "human_pop_dens", "value"] + 1)
  }

  if ("distance_PA" %in% unique(predictors_df$predictor)) {
    predictors_df[predictors_df$predictor == "distance_PA", "value"] <-
      sqrt(predictors_df[predictors_df$predictor == "distance_PA", "value"])
  }

  if ("deforestation_gaveau" %in% unique(predictors_df$predictor)) {
    predictors_df[predictors_df$predictor == "deforestation_gaveau", "value"] <-
      sqrt(predictors_df[predictors_df$predictor == "deforestation_gaveau", "value"])
  }

  for (nm in c("plantation_distance", "pulp_distance", "palm_distance")) {
    if (nm %in% unique(predictors_df$predictor)) {
      predictors_df[predictors_df$predictor == nm, "value"] <-
        log(predictors_df[predictors_df$predictor == nm, "value"] + 1)
    }
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
