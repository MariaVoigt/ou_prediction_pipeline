suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(raster))
suppressPackageStartupMessages(library(sp))
suppressPackageStartupMessages(library(rgeos))
suppressPackageStartupMessages(library(rgdal))

# ======================
# USER SETTINGS (EDIT!)
# ======================
# If you run this script by clicking "Source" in RStudio, the values below are used.
# If you run it via `Rscript prepare_predictors_from_aligned_rasters.R --help`, the CLI
# arguments are used instead.

RASTER_DIR <- "/path/to/aligned_rasters"
REFERENCE_RASTER <- "/path/to/reference.tif"
TRANSECTS_CSV <- "/path/to/transects.csv"
OUT_DIR <- "/path/to/pipeline_inputs"

# Comma-separated years, e.g. "1999,2000,2015"
YEARS <- "1999,2000,2001,2015"

# Buffer around transects for extraction (km)
BUFFER_KM <- 20

# Optional: shapefile polygon used to restrict the grid (leave as "" to disable)
MASK_SHP <- ""

# Prefix for generated grid ids
GRID_ID_PREFIX <- "grid_"

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !is.na(x)) x else y

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
script_dir <- if (!is.na(script_path) && nzchar(script_path)) dirname(script_path) else getwd()

source(file.path(script_dir, "predictor_config.R"))

short_name <- function(path) tools::file_path_sans_ext(basename(path))

read_mask_polygon <- function(mask_shp) {
  if (is.null(mask_shp) || mask_shp == "") return(NULL)
  dsn <- dirname(mask_shp)
  layer <- short_name(mask_shp)
  suppressWarnings(rgdal::readOGR(dsn = dsn, layer = layer, verbose = FALSE))
}

build_grid_from_reference <- function(reference_raster, mask_poly = NULL, id_prefix = "grid_") {
  r_ref <- raster(reference_raster)

  cells <- seq_len(ncell(r_ref))

  if (!is.null(mask_poly)) {
    xy <- xyFromCell(r_ref, cells)
    pts <- SpatialPoints(xy, proj4string = crs(r_ref))
    inside <- !is.na(over(pts, as(mask_poly, "SpatialPolygons")))
    cells <- cells[inside]
  }

  ext <- extent(r_ref)
  res_xy <- res(r_ref)
  ncol_r <- ncol(r_ref)

  rc <- rowColFromCell(r_ref, cells)
  row <- rc[, 1]
  col <- rc[, 2]

  x_min <- ext@xmin + (col - 1) * res_xy[1]
  x_max <- x_min + res_xy[1]
  y_max <- ext@ymax - (row - 1) * res_xy[2]
  y_min <- y_max - res_xy[2]

  grid <- data.frame(
    id = paste0(id_prefix, seq_along(cells)),
    cell = cells,
    x_start = x_min,
    y_start = y_min,
    x_end = x_max,
    y_end = y_max,
    stringsAsFactors = FALSE
  )

  grid
}

as_long_predictors <- function(df_wide, id_col = "id", year_col = "year") {
  out <- reshape2::melt(df_wide, id.vars = c(id_col, year_col), variable.name = "predictor", value.name = "value")
  out <- out %>% dplyr::select(id = !!id_col, year = !!year_col, predictor, value)
  out
}

resolve_predictor_rasters_for_year <- function(raster_dir, predictors, year) {
  out <- list()

  for (p in predictors) {
    p_year <- file.path(raster_dir, paste0(p, "_", year, ".tif"))
    p_static <- file.path(raster_dir, paste0(p, ".tif"))

    if (file.exists(p_year)) {
      out[[p]] <- p_year
    } else if (file.exists(p_static)) {
      out[[p]] <- p_static
    } else {
      stop(paste0("Missing raster for predictor '", p, "'. Looked for: ", p_year, " or ", p_static))
    }
  }

  out
}

stack_predictors <- function(predictor_raster_paths) {
  s <- raster::stack(unname(unlist(predictor_raster_paths)))
  names(s) <- names(predictor_raster_paths)
  s
}

extract_transect_predictors <- function(p_stack, transects_df, buffer_km) {
  req_cols <- c("id", "x_start", "y_start", "x_end", "y_end")
  miss <- setdiff(req_cols, names(transects_df))
  if (length(miss) > 0) stop(paste0("transects_df missing columns: ", paste(miss, collapse = ", ")))

  transects_df <- transects_df %>% mutate(.row_id = seq_len(n()))

  lines_list <- vector("list", nrow(transects_df))
  for (i in seq_len(nrow(transects_df))) {
    coords <- rbind(
      c(transects_df$x_start[i], transects_df$y_start[i]),
      c(transects_df$x_end[i], transects_df$y_end[i])
    )
    lines_list[[i]] <- Lines(list(Line(coords)), ID = as.character(transects_df$.row_id[i]))
  }

  sl <- SpatialLines(lines_list, proj4string = crs(p_stack))
  sldf <- SpatialLinesDataFrame(sl, data = transects_df, match.ID = TRUE)

  buf_m <- buffer_km * 1000
  sldf_buf <- rgeos::gBuffer(sldf, width = buf_m, byid = TRUE, capStyle = "ROUND", quadsegs = 10)

  raw <- raster::extract(p_stack, sldf_buf, df = TRUE)

  raw <- raw %>%
    dplyr::rename(.row_id = ID) %>%
    dplyr::mutate(.row_id = as.integer(.row_id)) %>%
    dplyr::inner_join(transects_df %>% dplyr::select(id, year, .row_id), by = ".row_id")

  agg <- raw %>%
    dplyr::group_by(id, year) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(names(p_stack)), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  as_long_predictors(agg, id_col = "id", year_col = "year")
}

extract_grid_predictors <- function(p_stack, grid_df, year) {
  vals <- raster::extract(p_stack, grid_df$cell)
  vals <- as.data.frame(vals)
  vals$id <- grid_df$id
  vals$year <- year
  vals <- vals %>% dplyr::select(id, year, dplyr::all_of(names(p_stack)))

  as_long_predictors(vals, id_col = "id", year_col = "year")
}

option_list <- list(
  make_option(c("--raster-dir"), type = "character", help = "Folder with aligned predictor .tif files"),
  make_option(c("--reference-raster"), type = "character", help = "Reference raster that defines grid (CRS/res/extent)"),
  make_option(c("--transects-csv"), type = "character", help = "CSV with transects (id, x_start,y_start,x_end,y_end, year, length_km, nr_nests, nest_decay, ...)"),
  make_option(c("--out-dir"), type = "character", help = "Output directory"),
  make_option(c("--years"), type = "character", help = "Comma-separated years to prepare grid predictors for (e.g. 1999,2000,2015)"),
  make_option(c("--buffer-km"), type = "double", default = 20, help = "Buffer (km) around transects for extraction"),
  make_option(c("--mask-shp"), type = "character", default = "", help = "Optional polygon shapefile to limit grid cells"),
  make_option(c("--grid-id-prefix"), type = "character", default = "grid_", help = "Prefix for grid ids")
)

opts <- NULL
if (interactive()) {
  opts <- list(
    raster_dir = RASTER_DIR,
    reference_raster = REFERENCE_RASTER,
    transects_csv = TRANSECTS_CSV,
    out_dir = OUT_DIR,
    years = YEARS,
    buffer_km = BUFFER_KM,
    mask_shp = MASK_SHP,
    grid_id_prefix = GRID_ID_PREFIX
  )
} else {
  opts <- parse_args(OptionParser(option_list = option_list))
}

if (is.null(opts$raster_dir) || is.null(opts$reference_raster) || is.null(opts$transects_csv) || is.null(opts$out_dir) || is.null(opts$years)) {
  stop("Missing required args. Use --help")
}

dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

transects <- read.csv(opts$transects_csv, stringsAsFactors = FALSE)
if (!all(c("id", "year") %in% names(transects))) stop("transects_csv must include at least: id, year")
transects$year <- as.integer(transects$year)

geography_obs <- transects %>%
  dplyr::select(id, x_start, y_start, x_end, y_end) %>%
  distinct()

write.csv(geography_obs, file.path(opts$out_dir, "geography_observation.csv"), row.names = FALSE)
write.csv(transects, file.path(opts$out_dir, "transects.csv"), row.names = FALSE)

predictors_obs <- NULL
transect_years <- sort(unique(transects$year))

for (y in transect_years) {
  p_paths <- resolve_predictor_rasters_for_year(opts$raster_dir, MODEL_PREDICTORS, y)
  p_stack <- stack_predictors(p_paths)

  transects_y <- transects[transects$year == y, , drop = FALSE]
  long_y <- extract_transect_predictors(p_stack, transects_y, buffer_km = opts$buffer_km)
  predictors_obs <- dplyr::bind_rows(predictors_obs, long_y)
}

predictors_obs <- apply_predictor_transforms(predictors_obs)

pred_obs_out <- file.path(opts$out_dir, paste0("predictors_observation_", as.integer(opts$buffer_km), ".csv"))
write.csv(predictors_obs, pred_obs_out, row.names = FALSE)

mask_poly <- read_mask_polygon(opts$mask_shp)

grid <- build_grid_from_reference(opts$reference_raster, mask_poly = mask_poly, id_prefix = opts$grid_id_prefix)

years_to_prepare <- as.integer(strsplit(opts$years, ",")[[1]])

for (y in years_to_prepare) {
  geography_y <- grid %>% dplyr::select(id, x_start, y_start, x_end, y_end)
  write.csv(geography_y, file.path(opts$out_dir, paste0("geography_", y, ".csv")), row.names = FALSE)

  p_paths <- resolve_predictor_rasters_for_year(opts$raster_dir, MODEL_PREDICTORS, y)
  p_stack <- stack_predictors(p_paths)

  long_grid <- extract_grid_predictors(p_stack, grid, year = y)
  long_grid <- apply_predictor_transforms(long_grid)

  write.csv(long_grid, file.path(opts$out_dir, paste0("predictors_abundance_", y, ".csv")), row.names = FALSE)
}
