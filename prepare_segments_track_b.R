# =============================================================================
# prepare_segments_track_b.R
#
# Track B (spatio-temporal degradation / re-occupation) data preparation.
#
# Builds the SEGMENT-YEAR table that the longitudinal NB GLMM consumes:
#   - splits each transect line into fixed-length segments
#   - assigns nest GPS points to the segment they fall on
#   - extracts per-segment buffered covariates (incl. degradation timing)
#   - constructs the effort + nest-parameter offset
#
# This is the piece the original pipeline does NOT have. It deliberately works
# at SEGMENT scale with a TIGHT buffer, so within-transect degradation
# contrasts are preserved instead of being averaged away.
#
# Output: one row per segment-year, ready for glmmTMB / mgcv.
#
# CRS assumption: everything is in a planar CRS in METERS (e.g. AEA), the same
# CRS as the aligned rasters — identical to the rest of the pipeline.
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(optparse)
})

# ----------------------------------------------------------------------------
# USER SETTINGS (EDIT!) — or override on the command line (see bottom).
# ----------------------------------------------------------------------------
DEFAULTS <- list(
  transects_shp   = "/path/to/transects_lines.shp",  # LINESTRING per transect, must carry a transect id + year
  nests_shp       = "/path/to/nests_points.shp",      # POINT per nest, must carry a year (+ ideally a transect id)
  raster_dir      = "/path/to/aligned_rasters",       # same aligned rasters as Track A
  logging_year_rast = "/path/to/year_last_logged.tif",# raster: most recent year a pixel was logged (NA = never)
  out_dir         = "/path/to/track_b_inputs",

  # Environmental covariates to extract per segment (named explicitly, not glob).
  # For each name p the script tries p_<year>.tif (time-varying) then p.tif (static).
  #   - dem, slope          : static -> dem.tif, slope.tif
  #   - distance_road       : if roads change, supply distance_road_<year>.tif;
  #                           else a single distance_road.tif (but see caveat:
  #                           logging roads ARE part of the degradation signal)
  #   - add forest-type binary layers (e.g. lowland_forest.tif) the same way
  env_covariates  = c("dem", "slope", "distance_road"),

  id_col          = "transect_id",   # column naming each transect uniquely (on transects)
  nest_id_col     = "transect_id",   # column on the NESTS layer giving each nest's transect
  year_col        = "year",          # survey year column on BOTH transects and nests

  segment_length_m = 250,            # 200-250 m recommended (see README discussion)
  buffer_m         = 500,            # TIGHT buffer around each segment (NOT 3 km!)

  # nest -> density offset parameters (orangutan nest method).
  # density = N / (L * 2*w * p * r * t)  -> all the constants go in a log-offset.
  strip_half_width_m = 25,           # w: effective half-strip width (one side), meters
  prop_builders      = 0.9,          # p: proportion of population that builds nests
  nest_prod_rate     = 1.0,          # r: nests built per individual per day
  nest_decay_days    = 200           # t: mean nest decay time (days). CHECK per-year/per-condition!
)

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

# Split a single LINESTRING into fixed-length segments (last segment may be shorter).
segment_one_line <- function(line_geom, seg_len_m, tid, yr) {
  total_len <- as.numeric(st_length(line_geom))
  if (total_len == 0) return(NULL)
  n_seg <- max(1L, ceiling(total_len / seg_len_m))
  # densify then cut by cumulative distance
  breaks <- seq(0, total_len, length.out = n_seg + 1L)
  coords <- st_coordinates(line_geom)[, c("X", "Y"), drop = FALSE]
  # cumulative distance along the vertices
  d <- c(0, cumsum(sqrt(rowSums(diff(coords)^2))))
  interp <- function(dist) {
    j <- findInterval(dist, d, all.inside = TRUE)
    f <- (dist - d[j]) / (d[j + 1] - d[j])
    f[!is.finite(f)] <- 0
    cbind(coords[j, 1] + f * (coords[j + 1, 1] - coords[j, 1]),
          coords[j, 2] + f * (coords[j + 1, 2] - coords[j, 2]))
  }
  segs <- lapply(seq_len(n_seg), function(k) {
    pts <- interp(c(breaks[k], breaks[k + 1]))
    g <- st_linestring(pts)
    list(
      seg_id = sprintf("%s__y%s__seg%02d", tid, yr, k),
      transect_id = tid,
      year = yr,
      seg_index = k,
      length_km = as.numeric(st_length(st_sfc(g))) / 1000,
      geometry = g
    )
  })
  segs
}

build_segments <- function(transects, id_col, year_col, seg_len_m) {
  rows <- list()
  for (i in seq_len(nrow(transects))) {
    tid <- as.character(transects[[id_col]][i])
    yr  <- transects[[year_col]][i]
    geom <- st_geometry(transects)[[i]]
    # handle MULTILINESTRING by merging to a single line where possible
    if (inherits(geom, "MULTILINESTRING")) geom <- st_line_merge(st_sfc(geom, crs = st_crs(transects)))[[1]]
    segs <- segment_one_line(st_sfc(geom, crs = st_crs(transects))[[1]], seg_len_m, tid, yr)
    rows <- c(rows, segs)
  }
  geoms <- st_sfc(lapply(rows, `[[`, "geometry"), crs = st_crs(transects))
  st_sf(
    seg_id      = vapply(rows, `[[`, "", "seg_id"),
    transect_id = vapply(rows, `[[`, "", "transect_id"),
    year        = vapply(rows, function(r) r$year, numeric(1)),
    seg_index   = vapply(rows, function(r) r$seg_index, numeric(1)),
    length_km   = vapply(rows, function(r) r$length_km, numeric(1)),
    geometry    = geoms
  )
}

# Assign each nest to a segment. Nests carry a transect id, so we match WITHIN
# the correct transect+year first (robust), then snap to the nearest segment of
# that transect. This avoids a nest being grabbed by a neighbouring transect's
# segment, which plain nearest-feature snapping can do where lines run close.
assign_nests_to_segments <- function(nests, segments, id_col, year_col,
                                      nest_id_col, max_snap_m = 100) {
  counts <- segments %>% st_drop_geometry() %>% transmute(seg_id, nr_nests = 0L)
  for (yr in sort(unique(segments$year))) {
    for (tid in unique(segments$transect_id[segments$year == yr])) {
      seg_ty <- segments[segments$year == yr & segments$transect_id == tid, ]
      nest_ty <- nests[nests[[year_col]] == yr &
                         as.character(nests[[nest_id_col]]) == tid, ]
      if (nrow(nest_ty) == 0 || nrow(seg_ty) == 0) next
      nrst <- st_nearest_feature(nest_ty, seg_ty)
      dist <- as.numeric(st_distance(nest_ty, seg_ty[nrst, ], by_element = TRUE))
      keep <- dist <= max_snap_m
      tab  <- table(seg_ty$seg_id[nrst[keep]])
      idx  <- match(names(tab), counts$seg_id)
      counts$nr_nests[idx] <- counts$nr_nests[idx] + as.integer(tab)
      if (any(!keep)) {
        message(sprintf("  transect %s year %s: %d nests > %dm from segment, dropped",
                        tid, yr, sum(!keep), max_snap_m))
      }
    }
  }
  counts
}

# Per-segment buffered mean of a continuous/binary raster (binary -> proportion cover).
extract_buffered_mean <- function(segments, rast, buffer_m, colname) {
  buf <- st_buffer(segments, dist = buffer_m)
  v   <- terra::extract(rast, terra::vect(buf), fun = mean, na.rm = TRUE, ID = FALSE)
  setNames(data.frame(seg_id = segments$seg_id, val = v[[1]]), c("seg_id", colname))
}

# Degradation TIMING covariates from a "year last logged" raster.
#   prop_degraded_1yr_prior : share of buffer logged in window (survey_year-1 .. survey_year)
#   years_since_last_logging: survey_year - (most recent logged year in buffer)
#   ever_logged             : 1 if any pixel in buffer was logged on/before survey year, else 0
#
# Never-logged segments get ever_logged = 0 and years_since_last_logging set to a
# finite REFERENCE value (0) so the row is NOT dropped. The model must include
# ever_logged and fit the recovery smooth with `by = ever_logged`, so the smooth
# is estimated ONLY over logged segments and never-logged ones act as the flat
# reference baseline (see model block at the end / README).
extract_degradation_timing <- function(segments, logging_year_rast, buffer_m,
                                        never_logged_ref = 0) {
  buf <- st_buffer(segments, dist = buffer_m)
  res <- lapply(seq_len(nrow(segments)), function(i) {
    yr  <- segments$year[i]
    cell_years <- terra::extract(logging_year_rast, terra::vect(buf[i, ]), ID = FALSE)[[1]]
    logged_years <- cell_years[is.finite(cell_years) & cell_years <= yr]
    if (length(logged_years) == 0) {
      return(data.frame(seg_id = segments$seg_id[i],
                        prop_degraded_1yr_prior = 0,
                        years_since_last_logging = never_logged_ref,
                        ever_logged = 0L))
    }
    prop_recent <- mean(cell_years[is.finite(cell_years)] %in% c(yr - 1, yr))
    last_logged <- max(logged_years)
    data.frame(seg_id = segments$seg_id[i],
               prop_degraded_1yr_prior = prop_recent,
               years_since_last_logging = yr - last_logged,
               ever_logged = 1L)
  })
  bind_rows(res)
}

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
run <- function(cfg) {
  dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

  message("Reading transects + nests ...")
  transects <- st_read(cfg$transects_shp, quiet = TRUE)
  nests     <- st_read(cfg$nests_shp, quiet = TRUE)
  stopifnot(cfg$id_col %in% names(transects), cfg$year_col %in% names(transects))
  stopifnot(cfg$year_col %in% names(nests), cfg$nest_id_col %in% names(nests))
  # force common planar CRS
  if (st_crs(nests) != st_crs(transects)) nests <- st_transform(nests, st_crs(transects))

  message(sprintf("Segmenting transects into %dm segments ...", cfg$segment_length_m))
  segments <- build_segments(transects, cfg$id_col, cfg$year_col, cfg$segment_length_m)
  message(sprintf("  -> %d segment-year rows", nrow(segments)))

  message("Assigning nests to segments ...")
  counts <- assign_nests_to_segments(nests, segments, cfg$id_col, cfg$year_col,
                                     cfg$nest_id_col)

  message("Extracting per-segment degradation timing ...")
  logyr <- terra::rast(cfg$logging_year_rast)
  deg   <- extract_degradation_timing(segments, logyr, cfg$buffer_m)

  # Environmental covariates, named explicitly in ENV_COVARIATES (config) so the
  # set is auditable rather than "whatever .tif is in the folder". For each name
  # `p` we look for p_<year>.tif first (time-varying, e.g. distance_road if roads
  # change) then fall back to p.tif (static, e.g. dem, slope). Buffered mean works
  # for continuous layers AND binary class-cover (0/1 -> proportion cover).
  env_tabs <- list()
  for (nm in cfg$env_covariates) {
    # build per-segment values, picking the right (possibly year-specific) raster
    vals <- lapply(sort(unique(segments$year)), function(yr) {
      f_yr <- file.path(cfg$raster_dir, sprintf("%s_%d.tif", nm, yr))
      f_st <- file.path(cfg$raster_dir, sprintf("%s.tif", nm))
      f <- if (file.exists(f_yr)) f_yr else if (file.exists(f_st)) f_st else NA
      if (is.na(f)) { message(sprintf("  MISSING raster for %s (year %d) - skipped", nm, yr)); return(NULL) }
      seg_y <- segments[segments$year == yr, ]
      extract_buffered_mean(seg_y, terra::rast(f), cfg$buffer_m, nm)
    })
    env_tabs[[nm]] <- bind_rows(vals)
    message(sprintf("  covariate: %s", nm))
  }

  message("Assembling segment-year table + offset ...")
  seg_tab <- segments %>% st_drop_geometry() %>%
    left_join(counts, by = "seg_id") %>%
    mutate(nr_nests = ifelse(is.na(nr_nests), 0L, nr_nests)) %>%
    left_join(deg, by = "seg_id")
  for (t in env_tabs) seg_tab <- left_join(seg_tab, t, by = "seg_id")

  # ---- offset: log of the density-conversion denominator ----
  # area surveyed (km^2) = length_km * (2 * strip_half_width_m / 1000)
  # detectability denom (per km^2 nest->ind) = p * r * t   (decay in days)
  area_km2 <- seg_tab$length_km * (2 * cfg$strip_half_width_m / 1000)
  denom    <- cfg$prop_builders * cfg$nest_prod_rate * cfg$nest_decay_days
  seg_tab$offset_term <- log(area_km2 * denom)   # so exp(linpred) = density (ind/km^2)

  # segment center coords (for an optional residual spatial term in Track B)
  ctr <- st_coordinates(st_centroid(st_geometry(segments)))
  seg_tab$x_km <- ctr[, 1] / 1000
  seg_tab$y_km <- ctr[, 2] / 1000

  # baseline nest count per transect (partial regression-to-the-mean control)
  baseline_year <- min(seg_tab$year)
  base <- seg_tab %>% filter(year == baseline_year) %>%
    group_by(transect_id) %>% summarise(baseline_nests = sum(nr_nests), .groups = "drop")
  seg_tab <- seg_tab %>% left_join(base, by = "transect_id") %>%
    mutate(baseline_nests_c = baseline_nests - mean(baseline_nests, na.rm = TRUE))

  # ever_logged must be an ORDERED factor for mgcv's `by =` (so the recovery
  # smooth is zero on never-logged segments and only takes shape on logged ones).
  seg_tab$ever_logged <- ordered(seg_tab$ever_logged, levels = c(0, 1))

  out_csv <- file.path(cfg$out_dir, "segment_year_table.csv")
  write.csv(seg_tab, out_csv, row.names = FALSE)
  st_write(segments, file.path(cfg$out_dir, "segments_geometry.gpkg"),
           delete_dsn = TRUE, quiet = TRUE)

  # ---- quick diagnostics you actually need to look at ----
  message("\n================ DIAGNOSTICS ================")
  message(sprintf("segment-year rows: %d | zero-nest rows: %d (%.0f%%)",
                  nrow(seg_tab), sum(seg_tab$nr_nests == 0),
                  100 * mean(seg_tab$nr_nests == 0)))
  message(sprintf("buffer %dm on %dm segments -> buffer/seg length ratio %.1f",
                  cfg$buffer_m, cfg$segment_length_m, cfg$buffer_m / cfg$segment_length_m))
  if (cfg$buffer_m > cfg$segment_length_m) {
    message("  WARNING: buffer wider than segment length -> adjacent segments' buffers")
    message("           overlap, eroding the within-transect degradation contrast.")
    message("           Consider buffer_m <= segment_length_m.")
  }
  message(sprintf("years present: %s",
                  paste(sort(unique(seg_tab$year)), collapse = ", ")))
  n_logged <- sum(as.integer(as.character(seg_tab$ever_logged)) == 1)
  message(sprintf("logged segment-years: %d | never-logged: %d",
                  n_logged, nrow(seg_tab) - n_logged))
  if (n_logged < 30) {
    message("  WARNING: few logged segment-years -> the recovery smooth will be")
    message("           poorly estimated. Consider a quadratic instead of a spline.")
  }
  message(sprintf("wrote: %s", out_csv))
  message("=============================================\n")

  message("---- suggested Track B model (paste into your fitting script) ----")
  cat('
library(mgcv)
seg_tab$transect_id <- factor(seg_tab$transect_id)
seg_tab$year_f      <- factor(seg_tab$year)
# ever_logged is already an ordered factor in the table.

m <- gam(
  nr_nests ~ ever_logged                                  # level shift: logged vs never-logged baseline
           + s(years_since_last_logging, by = ever_logged, k = 4)  # recovery smooth, LOGGED segments only
           + prop_degraded_1yr_prior                       # acute/standing degradation (expect weak)
           + baseline_nests_c                              # partial regression-to-the-mean control
           + dem + slope + distance_road                   # static/other environmental covariates
           + s(transect_id, bs = "re")                     # repeated visits to same transect
           + s(year_f,      bs = "re")                     # survey-wide year variation (6 levels)
           + offset(offset_term),
  family = nb(),
  method = "REML",
  data   = seg_tab
)
summary(m)
# Read your hypothesis off the smooth: plot(m, select = <the by-ever_logged term>).
# Re-occupation = curve dips below 0 at small years_since_last_logging then climbs
# back toward 0 (= toward the never-logged baseline) as time passes.
#
# k = 4 keeps the smooth from over-wiggling with few logged years; raise only if
# gam.check() flags k too low. With very few logged segment-years, drop the smooth
# and use a plain quadratic: poly(years_since_last_logging, 2):ever_logged.
', sep = "")
}

# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
if (sys.nframe() == 0) {
  opt_list <- list(
    make_option("--transects-shp", default = DEFAULTS$transects_shp),
    make_option("--nests-shp",     default = DEFAULTS$nests_shp),
    make_option("--raster-dir",    default = DEFAULTS$raster_dir),
    make_option("--logging-year-rast", default = DEFAULTS$logging_year_rast),
    make_option("--out-dir",       default = DEFAULTS$out_dir),
    make_option("--segment-length-m", default = DEFAULTS$segment_length_m, type = "integer"),
    make_option("--buffer-m",      default = DEFAULTS$buffer_m, type = "integer")
  )
  args <- parse_args(OptionParser(option_list = opt_list))
  cfg <- modifyList(DEFAULTS, list(
    transects_shp = args$`transects-shp`,
    nests_shp     = args$`nests-shp`,
    raster_dir    = args$`raster-dir`,
    logging_year_rast = args$`logging-year-rast`,
    out_dir       = args$`out-dir`,
    segment_length_m = args$`segment-length-m`,
    buffer_m      = args$`buffer-m`
  ))
  run(cfg)
}
