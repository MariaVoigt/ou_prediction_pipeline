# =============================================================================
# fit_track_b.R
#
# Track B (spatio-temporal degradation / re-occupation) MODEL FITTING + DIAGNOSTICS.
#
# Consumes segment_year_table.csv from prepare_segments_track_b.R and:
#   1. fits the NB GAM with the by-factor recovery smooth
#   2. runs gam.check (incl. whether k is high enough for the recovery smooth)
#   3. tests residual spatial autocorrelation across segments
#   4. SENSITIVITY: refits without baseline_nests_c   (regression-to-the-mean)
#   5. SENSITIVITY: refits without distance_road       (mediator vs confounder)
#   6. plots + writes the recovery curve and coefficient tables
#
# The headline result is the RECOVERY SMOOTH: plot of s(years_since_last_logging)
# for logged segments. Your hypothesis = dip below 0 just after logging, climbing
# back toward 0 (the never-logged baseline) as years pass.
#
# Honest framing this script enforces in its output:
#   - the degradation CONTRAST is the defensible test (revisit was disturbance-
#     independent), so it gets the headline.
#   - the absolute TRAJECTORY level is biased by selection on baseline occupancy;
#     the baseline_nests_c sensitivity check quantifies how much.
#   - RECOLONISATION of never-occupied habitat is unobservable here (those
#     transects were dropped) -> recovery is censored at the baseline level.
# =============================================================================

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
  library(readr)
})

# ----------------------------------------------------------------------------
# USER SETTINGS (EDIT!)
# ----------------------------------------------------------------------------
INPUT_CSV  <- "/path/to/track_b_inputs/segment_year_table.csv"
OUTPUT_DIR <- "/path/to/track_b_model"
RECOVERY_K <- 4          # basis dim for the recovery smooth; lower if few logged years
USE_DISTANCE_ROAD <- TRUE # included by default; sensitivity check refits without it

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# Load + prepare types
# ----------------------------------------------------------------------------
message("Loading segment-year table ...")
d <- read_csv(INPUT_CSV, show_col_types = FALSE)

req <- c("nr_nests", "years_since_last_logging", "ever_logged",
         "prop_degraded_1yr_prior", "baseline_nests_c", "offset_term",
         "transect_id", "year", "x_km", "y_km")
missing <- setdiff(req, names(d))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

d$transect_id <- factor(d$transect_id)
d$year_f      <- factor(d$year)
# ever_logged must be an ORDERED factor with never-logged (0) as the reference
# level, so the by-smooth is zero on never-logged segments. (Written ordered by
# the prep script, but re-coerce defensively after CSV round-trip.)
d$ever_logged <- ordered(d$ever_logged, levels = c(0, 1))

n_logged <- sum(as.integer(as.character(d$ever_logged)) == 1)
message(sprintf("rows: %d | logged: %d | never-logged: %d | zero-nest: %d (%.0f%%)",
                nrow(d), n_logged, nrow(d) - n_logged,
                sum(d$nr_nests == 0), 100 * mean(d$nr_nests == 0)))
if (n_logged < 30) {
  message("NOTE: few logged segment-years. Falling back to a quadratic recovery")
  message("      term instead of a spline (set USE_QUADRATIC_RECOVERY below).")
}
USE_QUADRATIC_RECOVERY <- n_logged < 30

# ----------------------------------------------------------------------------
# Build the formula
# ----------------------------------------------------------------------------
recovery_term <- if (USE_QUADRATIC_RECOVERY) {
  "poly(years_since_last_logging, 2):ever_logged"
} else {
  sprintf("s(years_since_last_logging, by = ever_logged, k = %d)", RECOVERY_K)
}

env_terms <- c("dem", "slope")
if (USE_DISTANCE_ROAD) env_terms <- c(env_terms, "distance_road")
env_terms <- env_terms[env_terms %in% names(d)]  # only those actually present

build_formula <- function(include_baseline = TRUE, include_road = TRUE) {
  et <- c("dem", "slope")
  if (include_road && "distance_road" %in% names(d)) et <- c(et, "distance_road")
  et <- et[et %in% names(d)]
  rhs <- c(
    "ever_logged",
    recovery_term,
    "prop_degraded_1yr_prior",
    if (include_baseline) "baseline_nests_c" else NULL,
    et,
    "s(transect_id, bs = 're')",
    "s(year_f, bs = 're')",
    "offset(offset_term)"
  )
  as.formula(paste("nr_nests ~", paste(rhs, collapse = " + ")))
}

# ----------------------------------------------------------------------------
# 1. Fit the main model
# ----------------------------------------------------------------------------
message("\nFitting main model ...")
f_main <- build_formula(include_baseline = TRUE, include_road = USE_DISTANCE_ROAD)
message("  formula: ", deparse1(f_main))
m <- gam(f_main, family = nb(), method = "REML", data = d)

sink(file.path(OUTPUT_DIR, "model_summary.txt"))
cat("=== TRACK B MAIN MODEL ===\n\n"); print(summary(m))
cat("\n=== AIC ===\n"); print(AIC(m))
sink()
saveRDS(m, file.path(OUTPUT_DIR, "track_b_gam.rds"))

# parametric coefficient table
ptab <- as.data.frame(summary(m)$p.table)
ptab$term <- rownames(ptab)
write_csv(ptab, file.path(OUTPUT_DIR, "parametric_coefficients.csv"))

# smooth-term table (incl. the recovery smooth EDF/p if spline used)
if (!USE_QUADRATIC_RECOVERY) {
  stab <- as.data.frame(summary(m)$s.table)
  stab$term <- rownames(stab)
  write_csv(stab, file.path(OUTPUT_DIR, "smooth_terms.csv"))
}

# ----------------------------------------------------------------------------
# 2. gam.check: is k high enough? residual behaviour?
# ----------------------------------------------------------------------------
message("Running gam.check ...")
png(file.path(OUTPUT_DIR, "gam_check.png"), width = 1000, height = 1000)
par(mfrow = c(2, 2)); gam.check(m); dev.off()
sink(file.path(OUTPUT_DIR, "gam_check.txt")); gam.check(m); sink()

# ----------------------------------------------------------------------------
# 3. THE headline plot: recovery curve (logged segments)
# ----------------------------------------------------------------------------
message("Plotting recovery curve ...")
png(file.path(OUTPUT_DIR, "recovery_curve.png"), width = 900, height = 650)
if (!USE_QUADRATIC_RECOVERY) {
  # find the index of the by-ever_logged smooth over years_since_last_logging
  smooth_labels <- sapply(m$smooth, function(s) s$label)
  rec_idx <- which(grepl("years_since_last_logging", smooth_labels))
  plot(m, select = rec_idx,
       shade = TRUE, seWithMean = TRUE,
       main = "Re-occupation: deviation from never-logged baseline",
       xlab = "Years since last logging", ylab = "Effect on log nest density")
  abline(h = 0, lty = 2, col = "grey50")
} else {
  # quadratic: predict over a grid for logged segments
  rng <- range(d$years_since_last_logging[as.integer(as.character(d$ever_logged)) == 1])
  nd <- data.frame(years_since_last_logging = seq(rng[1], rng[2], length.out = 50))
  nd$ever_logged <- ordered(1, levels = c(0, 1))
  for (v in c("prop_degraded_1yr_prior", "baseline_nests_c", env_terms))
    nd[[v]] <- mean(d[[v]], na.rm = TRUE)
  nd$offset_term <- 0
  nd$transect_id <- d$transect_id[1]; nd$year_f <- d$year_f[1]
  pr <- predict(m, nd, type = "link", se.fit = TRUE,
                exclude = c("s(transect_id)", "s(year_f)"))
  plot(nd$years_since_last_logging, pr$fit, type = "l", lwd = 2,
       main = "Re-occupation (quadratic) — logged segments",
       xlab = "Years since last logging", ylab = "log nest density (centred)")
  lines(nd$years_since_last_logging, pr$fit + 2 * pr$se.fit, lty = 2)
  lines(nd$years_since_last_logging, pr$fit - 2 * pr$se.fit, lty = 2)
}
dev.off()

# ----------------------------------------------------------------------------
# 4. Residual spatial autocorrelation
#    Segments within a transect are handled by the RE; this checks whether
#    residuals are still spatially clustered ACROSS transects (would justify
#    adding s(x_km, y_km) as a NUISANCE term — not the Track A interpolation use).
# ----------------------------------------------------------------------------
message("Testing residual spatial autocorrelation ...")
d$resid_dev <- residuals(m, type = "deviance")
# collapse to transect centroids to test the between-transect pattern
tr <- d %>% group_by(transect_id) %>%
  summarise(x = mean(x_km), y = mean(y_km), r = mean(resid_dev), .groups = "drop")

moran_result <- tryCatch({
  # simple distance-band Moran's I without extra packages
  coords <- as.matrix(tr[, c("x", "y")])
  D <- as.matrix(dist(coords))
  W <- 1 / D; diag(W) <- 0; W[!is.finite(W)] <- 0
  W <- W / rowSums(W); W[!is.finite(W)] <- 0
  z <- tr$r - mean(tr$r)
  I <- (length(z) / sum(W)) * (t(z) %*% (W %*% z)) / sum(z^2)
  as.numeric(I)
}, error = function(e) NA_real_)

sink(file.path(OUTPUT_DIR, "spatial_autocorr.txt"))
cat("Between-transect Moran's I on deviance residuals:", round(moran_result, 4), "\n")
cat("(near 0 = no residual spatial structure; strongly positive => consider\n")
cat(" adding s(x_km, y_km) as a NUISANCE term and refitting.)\n")
sink()
message(sprintf("  Moran's I (between-transect): %.4f", moran_result))

# residual-vs-space plot
png(file.path(OUTPUT_DIR, "residuals_in_space.png"), width = 800, height = 700)
with(tr, {
  cols <- ifelse(r >= 0, "red", "blue")
  plot(x, y, cex = sqrt(abs(r)) + 0.3, col = cols, pch = 19,
       main = "Mean deviance residual per transect (red +, blue -)",
       xlab = "x (km)", ylab = "y (km)", asp = 1)
})
dev.off()

# ----------------------------------------------------------------------------
# 5. SENSITIVITY: regression-to-the-mean (with vs without baseline_nests_c)
# ----------------------------------------------------------------------------
message("Sensitivity: regression-to-the-mean (drop baseline_nests_c) ...")
m_nobase <- gam(build_formula(include_baseline = FALSE, include_road = USE_DISTANCE_ROAD),
                family = nb(), method = "REML", data = d)

extract_key <- function(model, label) {
  pt <- summary(model)$p.table
  key <- intersect(c("prop_degraded_1yr_prior", "ever_logged.L", "dem", "slope",
                     "distance_road"), rownames(pt))
  data.frame(model = label, term = key,
             estimate = pt[key, 1], se = pt[key, 2], p = pt[key, 4],
             row.names = NULL)
}
sens_base <- bind_rows(extract_key(m, "with_baseline"),
                       extract_key(m_nobase, "without_baseline"))
write_csv(sens_base, file.path(OUTPUT_DIR, "sensitivity_baseline.csv"))

# ----------------------------------------------------------------------------
# 6. SENSITIVITY: distance_road as mediator vs confounder (with vs without)
# ----------------------------------------------------------------------------
sens_road <- NULL
if (USE_DISTANCE_ROAD && "distance_road" %in% names(d)) {
  message("Sensitivity: distance_road mediator check (drop distance_road) ...")
  m_noroad <- gam(build_formula(include_baseline = TRUE, include_road = FALSE),
                  family = nb(), method = "REML", data = d)
  sens_road <- bind_rows(extract_key(m, "with_road"),
                         extract_key(m_noroad, "without_road"))
  write_csv(sens_road, file.path(OUTPUT_DIR, "sensitivity_road.csv"))
}

# ----------------------------------------------------------------------------
# Console summary
# ----------------------------------------------------------------------------
message("\n================ TRACK B SUMMARY ================")
message("Main model written to: ", file.path(OUTPUT_DIR, "track_b_gam.rds"))
message("Headline plot:         ", file.path(OUTPUT_DIR, "recovery_curve.png"))
message("\nKey terms (main model):")
print(extract_key(m, "main"), row.names = FALSE)
message("\nRegression-to-the-mean check (does dropping baseline_nests_c move terms?):")
print(sens_base, row.names = FALSE)
if (!is.null(sens_road)) {
  message("\ndistance_road check (does dropping it move the degradation terms? -> mediator):")
  print(sens_road, row.names = FALSE)
}
message("\nInterpretation reminders:")
message(" - Recovery curve dipping then returning toward 0 = re-occupation signal.")
message(" - If degradation terms move a LOT when baseline_nests_c is dropped, your")
message("   absolute trajectory is selection-sensitive; report the with-baseline version.")
message(" - If degradation terms shrink when distance_road is ADDED, roads may be a")
message("   mediator of the degradation effect, not a confounder; report both.")
message(" - Moran's I far above 0 -> refit adding s(x_km, y_km) as a nuisance term.")
message("=================================================\n")
