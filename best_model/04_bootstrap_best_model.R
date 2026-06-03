suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(mgcv))

MODEL_DIR <- "/path/to/your/output"
OUTPUT_DIR <- "/path/to/your/output"

N_BOOT <- 1000
SET_SEED <- 1

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
set.seed(SET_SEED)

best_fit <- readRDS(file.path(MODEL_DIR, "best_gam_nb_model.rds"))

# Coefficient uncertainty for a GAM via posterior (Bayesian) simulation:
# draw coefficients from N(beta_hat, Vp). This covers the parametric
# (environmental) terms AND the spatial-smooth basis coefficients.
beta_hat <- coef(best_fit)
Vp <- vcov(best_fit)

beta_draws <- mgcv::rmvn(N_BOOT, beta_hat, Vp)
colnames(beta_draws) <- names(beta_hat)

# Report only the parametric (interpretable) terms; smooth-basis coefficients
# are not individually meaningful.
param_terms <- rownames(summary(best_fit)$p.table)

beta_draws_param <- as.data.frame(beta_draws[, param_terms, drop = FALSE])

write.csv(beta_draws_param,
          file = file.path(OUTPUT_DIR, "best_model_posterior_coefficients.csv"),
          row.names = FALSE)

coef_summary <- lapply(param_terms, function(nm) {
  vals <- beta_draws_param[[nm]]
  data.frame(
    term     = nm,
    estimate = beta_hat[[nm]],
    q025     = as.numeric(quantile(vals, probs = 0.025, na.rm = TRUE)),
    q500     = as.numeric(quantile(vals, probs = 0.5,   na.rm = TRUE)),
    q975     = as.numeric(quantile(vals, probs = 0.975, na.rm = TRUE))
  )
}) %>% bind_rows()

write.csv(coef_summary,
          file = file.path(OUTPUT_DIR, "best_model_posterior_coef_summary.csv"),
          row.names = FALSE)

# Negative-binomial theta (dispersion) point estimate from the fitted GAM.
theta_hat <- best_fit$family$getTheta(TRUE)
theta_summary <- data.frame(theta_hat = theta_hat, n_boot = N_BOOT)
write.csv(theta_summary,
          file = file.path(OUTPUT_DIR, "best_model_theta_summary.csv"),
          row.names = FALSE)
