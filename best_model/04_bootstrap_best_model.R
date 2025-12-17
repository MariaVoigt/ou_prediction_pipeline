suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(MASS))

MODEL_DIR <- "/path/to/your/output"
OUTPUT_DIR <- "/path/to/your/output"

N_BOOT <- 200
SET_SEED <- 1

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
set.seed(SET_SEED)

predictors_obs <- read.csv(file.path(MODEL_DIR, "predictors_observation_scaled.csv"), stringsAsFactors = FALSE)
best_rhs_tbl <- read.csv(file.path(MODEL_DIR, "best_model_rhs.csv"), stringsAsFactors = FALSE)
best_rhs <- best_rhs_tbl$best_model_rhs[1]

fml <- as.formula(paste0("nr_nests ~ ", best_rhs, " + offset(offset_term)"))

best_fit <- readRDS(file.path(MODEL_DIR, "best_glm_nb_model.rds"))

# Design matrix used for simulation
X <- model.matrix(fml, data = predictors_obs)
beta_hat <- coef(best_fit)
V <- vcov(best_fit)

# Theta uncertainty (may be NA if not available)
theta_hat <- best_fit$theta
se_theta <- tryCatch(summary(best_fit)$SE.theta, error = function(e) NA_real_)

boot_coefs <- matrix(NA_real_, nrow = N_BOOT, ncol = length(beta_hat))
colnames(boot_coefs) <- names(beta_hat)

boot_theta <- rep(NA_real_, N_BOOT)

attempts <- 0
kept <- 0

while (kept < N_BOOT) {
  attempts <- attempts + 1

  beta_draw <- MASS::mvrnorm(n = 1, mu = beta_hat, Sigma = V)
  names(beta_draw) <- names(beta_hat)

  theta_draw <- theta_hat
  if (!is.na(se_theta)) {
    theta_draw <- rnorm(1, mean = theta_hat, sd = se_theta)
  }
  if (!is.finite(theta_draw) || theta_draw <= 0) next

  eta <- as.numeric(X %*% beta_draw) + predictors_obs$offset_term
  mu <- exp(eta)

  y_sim <- rnbinom(n = nrow(predictors_obs), mu = mu, size = theta_draw)

  sim_data <- predictors_obs
  sim_data$nr_nests <- y_sim

  fit_sim <- try(glm.nb(fml, data = sim_data, control = glm.control(maxit = 500)), silent = TRUE)
  if (inherits(fit_sim, "try-error")) next

  kept <- kept + 1
  boot_coefs[kept, ] <- coef(fit_sim)
  boot_theta[kept] <- fit_sim$theta
}

boot_coefs <- as.data.frame(boot_coefs)

write.csv(boot_coefs, file = file.path(OUTPUT_DIR, "best_model_bootstrap_coefficients.csv"), row.names = FALSE)

coef_summary <- lapply(names(beta_hat), function(nm) {
  vals <- boot_coefs[[nm]]
  data.frame(
    term = nm,
    estimate = beta_hat[[nm]],
    q025 = as.numeric(quantile(vals, probs = 0.025, na.rm = TRUE)),
    q500 = as.numeric(quantile(vals, probs = 0.5, na.rm = TRUE)),
    q975 = as.numeric(quantile(vals, probs = 0.975, na.rm = TRUE))
  )
}) %>% bind_rows()

write.csv(coef_summary, file = file.path(OUTPUT_DIR, "best_model_bootstrap_coef_summary.csv"), row.names = FALSE)

theta_summary <- data.frame(
  theta_hat = theta_hat,
  theta_q025 = as.numeric(quantile(boot_theta, probs = 0.025, na.rm = TRUE)),
  theta_q500 = as.numeric(quantile(boot_theta, probs = 0.5, na.rm = TRUE)),
  theta_q975 = as.numeric(quantile(boot_theta, probs = 0.975, na.rm = TRUE)),
  n_boot = N_BOOT,
  attempts = attempts
)

write.csv(theta_summary, file = file.path(OUTPUT_DIR, "best_model_bootstrap_theta_summary.csv"), row.names = FALSE)
