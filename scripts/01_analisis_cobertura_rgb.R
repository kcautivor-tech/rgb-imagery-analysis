# =============================================================================
# PROJECT: PanGaia - 12 Kanpa IMUs, RGB vegetation-cover screening
# AUTHOR: Karla Cautivo
# REVISION: 2026-11-10; VARI + histogram Otsu, provisional raw denominator D > 5
# Run from the RStudio project root: source('scripts/01_analisis_cobertura_rgb.R')
# Install dependencies once: install.packages(c('terra', 'sf'))
#
# This estimates HIGH-VARI cover, subject to ecological validation. It does not
# remove all illumination effects or establish true/live vegetation automatically.
# Existing QGIS ExG results are not inputs to threshold fitting.
# =============================================================================

pangaia_config <- function() list(
  input_dir = 'data/raw',
  boundary_file = 'data/raw/KA_IMU_boundaries.gpkg',
  boundary_layer = NULL,  # Set explicitly if the GeoPackage has multiple layers.
  manifest_file = 'data/raw/KA_image_manifest.csv',
  output_root = 'data/processed/VARI_Otsu_D5',
  plots_root = 'outputs/VARI_Otsu_D5',
  expected_imus = 12L, epsg = 7851L,
  min_denominator_raw = 5, # STRICT D > 5 in original 0-255 RGB units.
  explicit_nodata = -9999, # Optional extra sentinel. Embedded NoData is respected.
  histogram_bins = 512L, block_rows = 100L,
  min_valid_pixels = 100L, min_valid_polygon_percent = 95,
  min_separability = 0.50, extreme_cover_percent = c(1, 99)
)

require_spatial_packages <- function() {
  for (p in c('terra', 'sf')) if (!requireNamespace(p, quietly = TRUE))
    stop("Install dependencies first: install.packages(c('terra', 'sf'))")
}

# Pure R Otsu. Histogram bins are right closed: (lower, upper], with the first
# including its lower endpoint. Class 1 therefore uses VARI > threshold.
otsu_histogram <- function(counts, breaks) {
  stopifnot(length(breaks) == length(counts) + 1L, all(diff(breaks) > 0),
            all(is.finite(counts)), all(counts >= 0))
  if (sum(counts) == 0 || sum(counts > 0) < 2L) return(NULL)
  p <- counts / sum(counts)
  centres <- (head(breaks, -1) + tail(breaks, -1)) / 2
  w <- cumsum(p); mu <- cumsum(p * centres); total <- sum(p * centres)
  score <- rep(NA_real_, length(p)); eligible <- w > 0 & w < 1
  score[eligible] <- (total * w[eligible] - mu[eligible])^2 /
    (w[eligible] * (1 - w[eligible]))
  if (!any(is.finite(score))) return(NULL)
  best <- max(score, na.rm = TRUE)
  ties <- which(is.finite(score) & abs(score - best) <= max(1e-15, best * 1e-12))
  k <- ties[ceiling(length(ties) / 2)]
  list(threshold = breaks[k + 1L],
       separability = best / sum(p * (centres - total)^2))
}

# Work directly with original integer channels, without first dividing by 255.
# The common scale cancels algebraically. Direct arithmetic makes D > 5 and
# thresholds reproducible without normalization/temporary-raster rounding.
rgb_values <- function(rgb, cfg) {
  if (length(cfg$explicit_nodata)) rgb[rgb %in% cfg$explicit_nodata] <- NA_real_
  complete <- rowSums(is.finite(rgb)) == 3L
  in_range <- rowSums(is.finite(rgb) & rgb >= 0 & rgb <= 255) == 3L
  ok <- complete & in_range
  if (any(rgb[ok, , drop = FALSE] != round(rgb[ok, , drop = FALSE])))
    stop('This calibrated workflow requires integer 0-255 RGB, not normalized or fractional channels.')
  d <- rgb[, 2] + rgb[, 1] - rgb[, 3]
  n <- rgb[, 2] - rgb[, 1]
  list(valid = ok, D = d, N = n,
       valid_D = ok & is.finite(d) & d > cfg$min_denominator_raw)
}

parse_imu_number <- function(x) {
  z <- toupper(trimws(tools::file_path_sans_ext(basename(as.character(x)))))
  pattern <- '^(KA[_ -]*)?IMU[_ -]*0*([1-9][0-9]*)([^0-9].*)?$'
  ans <- rep(NA_integer_, length(z)); hit <- grepl(pattern, z)
  ans[hit] <- as.integer(sub(pattern, '\\2', z[hit]))
  numeric_id <- grepl('^[0-9]+$', z)
  ans[numeric_id] <- as.integer(z[numeric_id])
  ans
}

read_boundaries <- function(cfg) {
  if (!file.exists(cfg$boundary_file)) stop('Missing boundary file: ', cfg$boundary_file)
  layers <- sf::st_layers(cfg$boundary_file)$name
  layer <- cfg$boundary_layer
  if (is.null(layer)) {
    if (length(layers) != 1L) stop('Set boundary_layer to the IMU polygon layer.')
    layer <- layers[1]
  }
  b <- sf::st_read(cfg$boundary_file, layer = layer, quiet = TRUE)
  if (!'imu_id' %in% names(b)) stop('Boundary layer must contain imu_id.')
  b$imu_id <- as.character(b$imu_id)
  if (anyNA(b$imu_id) || any(!nzchar(trimws(b$imu_id)))) stop('Missing boundary imu_id.')
  if (!isTRUE(sf::st_crs(b) == sf::st_crs(cfg$epsg))) stop('Boundary CRS must be EPSG:7851.')
  if (any(!sf::st_is_valid(b)) || any(sf::st_is_empty(b)) ||
      any(!as.character(sf::st_geometry_type(b)) %in% c('POLYGON', 'MULTIPOLYGON')))
    stop('Boundary geometries must be valid, nonempty polygons.')
  ids <- unique(b$imu_id)
  if (length(ids) != cfg$expected_imus) stop('Expected exactly 12 distinct Kanpa IMUs.')
  geoms <- lapply(ids, function(id) sf::st_union(sf::st_geometry(b[b$imu_id == id, ]))[[1]])
  sf::st_sf(imu_id = ids, geometry = sf::st_sfc(geoms, crs = sf::st_crs(b)))
}

image_manifest <- function(b, cfg) {
  if (file.exists(cfg$manifest_file)) {
    m <- utils::read.csv(cfg$manifest_file, colClasses = 'character', check.names = FALSE)
    if (!all(c('imu_id', 'image_file', 'rgb_order_confirmed') %in% names(m)))
      stop('Manifest needs imu_id,image_file,rgb_order_confirmed.')
    if (anyNA(m$rgb_order_confirmed) || any(!toupper(trimws(m$rgb_order_confirmed)) %in% c('TRUE', 'FALSE')))
      stop('rgb_order_confirmed must be TRUE or FALSE.')
    m$rgb_order_confirmed <- toupper(trimws(m$rgb_order_confirmed)) == 'TRUE'
  } else {
    files <- list.files(cfg$input_dir, pattern = '\\.tiff?$', ignore.case = TRUE)
    bn <- parse_imu_number(b$imu_id); fn <- parse_imu_number(files)
    if (length(files) != cfg$expected_imus || anyNA(bn) || anyNA(fn) ||
        anyDuplicated(bn) || anyDuplicated(fn) || !setequal(bn, fn))
      stop('Cannot safely map images. Supply an explicit KA_image_manifest.csv.')
    m <- data.frame(imu_id = b$imu_id, image_file = files[match(bn, fn)], rgb_order_confirmed = FALSE)
  }
  if (nrow(m) != cfg$expected_imus || anyNA(m$imu_id) || anyNA(m$image_file) ||
      any(!nzchar(m$image_file)) || anyDuplicated(m$imu_id) || anyDuplicated(m$image_file) ||
      !setequal(m$imu_id, b$imu_id)) stop('Manifest must map all 12 IDs to distinct files.')
  number <- parse_imu_number(m$imu_id)
  if (!anyNA(number)) m <- m[order(number), ]
  m$output_key <- sprintf('KA_%02d', seq_len(nrow(m)))
  m
}

prepare_image <- function(entry, boundary, cfg, audit_dir) {
  f <- entry$image_file
  if (!grepl('^([A-Za-z]:[/\\\\]|/|\\\\\\\\)', f)) f <- file.path(cfg$input_dir, f)
  if (!file.exists(f)) stop('Missing image: ', f)
  img <- terra::rast(f)
  if (terra::nlyr(img) != 3L) stop('Image must contain exactly three RGB bands.')
  if (!isTRUE(sf::st_crs(terra::crs(img)) == sf::st_crs(cfg$epsg))) stop('Image CRS must be EPSG:7851.')
  info <- terra::describe(f, print = FALSE)
  writeLines(info, file.path(audit_dir, paste0(entry$output_key, '_source_metadata.txt')))
  colors <- sub('.*ColorInterp=([A-Za-z]+).*', '\\1',
                grep('^Band [0-9]+ .*ColorInterp=', info, value = TRUE))
  flags <- 'PROVISIONAL_D_CUTOFF_REVIEW_REQUIRED'
  if (length(colors) == 3L && setequal(colors, c('Red', 'Green', 'Blue'))) {
    img <- img[[match(c('Red', 'Green', 'Blue'), colors)]]
  } else if (isTRUE(entry$rgb_order_confirmed)) {
    flags <- c(flags, 'RGB_ORDER_USER_CONFIRMED')
  } else stop('RGB order not tagged. Verify bands 1=R,2=G,3=B and confirm in the manifest.')
  # Preserve the source mask/NoData. NEVER assign NAflag(img) <- 0.
  if (any(grepl('NoData Value=0($|[[:space:]])', info))) flags <- c(flags, 'SOURCE_NODATA_ZERO')
  names(img) <- c('R', 'G', 'B')
  poly <- terra::vect(boundary)
  img <- terra::crop(img, poly, snap = 'out')
  img <- terra::extend(img, terra::ext(poly), snap = 'out')
  inside <- terra::rasterize(poly, img[[1]], field = 1, background = NA, touches = FALSE)
  list(rgb = img, inside = inside, flags = flags)
}

# Every valid integer RGB pixel contributes to a compact exact joint histogram:
# denominator D = -255..510 and numerator N = -255..255. No image sampling.
count_rgb_distribution <- function(img, inside, cfg) {
  joint <- numeric(766 * 511); polygon_n <- 0
  terra::readStart(img); terra::readStart(inside)
  on.exit({terra::readStop(img); terra::readStop(inside)}, add = TRUE)
  for (row in seq.int(1L, terra::nrow(img), by = cfg$block_rows)) {
    nr <- min(cfg$block_rows, terra::nrow(img) - row + 1L)
    ins <- !is.na(terra::readValues(inside, row = row, nrows = nr))
    rgb <- terra::readValues(img, row = row, nrows = nr, mat = TRUE)
    v <- rgb_values(rgb[ins, , drop = FALSE], cfg)
    polygon_n <- polygon_n + sum(ins)
    ok <- v$valid
    index <- (v$D[ok] + 255) * 511 + (v$N[ok] + 255) + 1
    joint <- joint + tabulate(index, nbins = length(joint))
  }
  ix <- which(joint > 0) - 1L
  df <- data.frame(D = ix %/% 511 - 255, N = ix %% 511 - 255, count = joint[ix + 1L])
  df$VARI <- df$N / ifelse(df$D > 0, df$D, NA_real_)
  list(distribution = df, polygon_pixels = polygon_n)
}

fit_distribution <- function(df, cfg) {
  z <- df[df$D > cfg$min_denominator_raw, ]
  if (sum(z$count) < cfg$min_valid_pixels) stop('Too few usable VARI pixels.')
  limits <- range(z$VARI)
  if (!all(is.finite(limits)) || limits[2] <= limits[1]) stop('Constant/empty VARI: Otsu is undefined.')
  breaks <- seq(limits[1], limits[2], length.out = cfg$histogram_bins + 1L)
  bins <- cut(z$VARI, breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)
  totals <- rowsum(z$count, bins, reorder = FALSE)
  counts <- numeric(cfg$histogram_bins); counts[as.integer(rownames(totals))] <- totals[, 1]
  stopifnot(sum(counts) == sum(z$count))
  fit <- otsu_histogram(counts, breaks)
  if (is.null(fit)) stop('Otsu requires at least two occupied bins.')
  fit$histogram <- data.frame(lower = head(breaks, -1), upper = tail(breaks, -1), count = counts)
  fit$limits <- limits
  fit
}

# Classify in blocks using the SAME raw arithmetic as histogram fitting.
# Invalid/outside/excluded pixels are NA, never converted to class zero.
write_binary_mask <- function(img, inside, threshold, path, cfg) {
  result <- terra::rast(img[[1]]); names(result) <- 'high_VARI_1_low_VARI_0'
  terra::readStart(img); terra::readStart(inside)
  on.exit({terra::readStop(img); terra::readStop(inside)}, add = TRUE)
  terra::writeStart(result, path, overwrite = FALSE,
    wopt = list(datatype = 'INT1U', NAflag = 255, gdal = c('COMPRESS=LZW')))
  open <- TRUE
  on.exit(if (open) try(terra::writeStop(result), silent = TRUE), add = TRUE)
  counts <- c(high = 0, low = 0)
  for (row in seq.int(1L, terra::nrow(img), by = cfg$block_rows)) {
    nr <- min(cfg$block_rows, terra::nrow(img) - row + 1L)
    rgb <- terra::readValues(img, row = row, nrows = nr, mat = TRUE)
    ins <- !is.na(terra::readValues(inside, row = row, nrows = nr))
    rgb[!ins, ] <- NA_real_
    v <- rgb_values(rgb, cfg); ok <- ins & v$valid_D
    binary <- rep(NA_real_, nrow(rgb))
    binary[ok] <- as.integer(v$N[ok] / v$D[ok] > threshold)
    counts <- counts + c(sum(binary == 1, na.rm = TRUE), sum(binary == 0, na.rm = TRUE))
    terra::writeValues(result, binary, row, nr)
  }
  result <- terra::writeStop(result); open <- FALSE
  list(raster = result, counts = counts)
}

# Preview samples cells on a regular grid; coverage statistics use every pixel.
make_preview <- function(img, inside, binary, fit, entry, metrics, file, cfg) {
  step <- max(1L, ceiling(max(terra::nrow(img), terra::ncol(img)) / 600))
  rr <- seq.int(1L, terra::nrow(img), by = step); cc <- seq.int(1L, terra::ncol(img), by = step)
  cells <- as.vector(t(outer(rr, cc, function(r, c) (r - 1) * terra::ncol(img) + c)))
  rgb <- as.matrix(terra::extract(img, cells))
  ins <- !is.na(terra::extract(inside, cells)[, 1])
  rgb[!ins, ] <- NA_real_; v <- rgb_values(rgb, cfg)
  pred <- terra::extract(binary, cells)[, 1]
  draw <- function(colors, title) {
    graphics::plot.new(); graphics::plot.window(xlim = c(0, length(cc)/length(rr)), ylim = c(0, 1), asp = 1)
    graphics::rasterImage(as.raster(matrix(colors, nrow = length(rr), byrow = TRUE)),
      0, 0, length(cc)/length(rr), 1, interpolate = FALSE)
    graphics::title(main = title, cex.main = 0.9)
  }
  type <- if (capabilities('cairo')) 'cairo-png' else getOption('bitmapType')
  grDevices::png(file, width = 1900, height = 1500, res = 160, type = type)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(5, 5.5, 3, 1))
  colors <- rep('white', length(cells)); ok <- ins & v$valid
  colors[ok] <- grDevices::rgb(rgb[ok, 1], rgb[ok, 2], rgb[ok, 3], maxColorValue = 255)
  draw(colors, paste(entry$imu_id, 'RGB - regular sampled preview'))
  colors <- rep('white', length(cells)); colors[ins] <- '#D9D9D9'
  colors[!is.na(pred) & pred == 0] <- '#BE8A58'; colors[!is.na(pred) & pred == 1] <- '#228B22'
  draw(colors, sprintf('High VARI: %.2f%% of usable pixels', metrics$coverage_percent))
  h <- fit$histogram
  graphics::plot((h$lower + h$upper)/2, h$count, type = 'h', xlab = 'VARI',
    ylab = 'Pixel count', main = 'Full-pixel histogram and Otsu threshold')
  graphics::abline(v = fit$threshold, col = 'red3', lwd = 2)
  graphics::plot.new()
  graphics::legend('topleft', legend = c('1: high VARI', '0: low VARI', 'Excluded / NoData'),
    fill = c('#228B22', '#BE8A58', '#D9D9D9'), bty = 'n')
  graphics::text(.5, .4, sprintf('D > %g (raw RGB units)\nUsable polygon: %.2f%%\nThreshold: %.5f\nProvisional classification; review required',
    cfg$min_denominator_raw, metrics$valid_polygon_percent, fit$threshold), cex = .95)
}

process_imu <- function(entry, boundary, cfg, run_dir, plot_dir) {
  audit <- file.path(run_dir, 'audit')
  result <- data.frame(imu_id = entry$imu_id, image_file = entry$image_file,
    method = 'VARI_Otsu_raw_D_cutoff', status = 'ERROR', denominator_cutoff_raw = cfg$min_denominator_raw,
    threshold = NA_real_, otsu_separability = NA_real_,
    polygon_area_m2 = as.numeric(sf::st_area(boundary)), pixel_area_m2 = NA_real_,
    polygon_pixels = NA_real_, valid_rgb_pixels = NA_real_, invalid_rgb_pixels = NA_real_,
    nonpositive_denominator_pixels = NA_real_, small_positive_denominator_pixels = NA_real_,
    excluded_denominator_pixels = NA_real_, valid_vari_pixels = NA_real_,
    high_vari_pixels = NA_real_, low_vari_pixels = NA_real_, valid_polygon_percent = NA_real_,
    coverage_percent = NA_real_, high_vari_polygon_percent = NA_real_, high_vari_area_m2 = NA_real_,
    vari_min = NA_real_, vari_max = NA_real_, binary_geotiff = '', png = '', qc_flags = '', error = '')
  flags <- 'PROVISIONAL_D_CUTOFF_REVIEW_REQUIRED'
  tryCatch({
    message('  Checking RGB / CRS and masking boundary')
    x <- prepare_image(entry, boundary, cfg, audit); flags <- x$flags
    message('  Counting all valid pixels in blocks')
    counted <- count_rgb_distribution(x$rgb, x$inside, cfg); df <- counted$distribution
    result$polygon_pixels <- counted$polygon_pixels
    result$valid_rgb_pixels <- sum(df$count)
    result$invalid_rgb_pixels <- result$polygon_pixels - result$valid_rgb_pixels
    result$nonpositive_denominator_pixels <- sum(df$count[df$D <= 0])
    result$small_positive_denominator_pixels <- sum(df$count[df$D > 0 & df$D <= cfg$min_denominator_raw])
    result$excluded_denominator_pixels <- result$nonpositive_denominator_pixels + result$small_positive_denominator_pixels
    result$valid_vari_pixels <- sum(df$count[df$D > cfg$min_denominator_raw])
    if (result$polygon_pixels == 0) stop('No pixel centres inside the polygon.')
    result$valid_polygon_percent <- 100 * result$valid_vari_pixels / result$polygon_pixels
    result$pixel_area_m2 <- prod(terra::res(x$rgb))
    if (result$invalid_rgb_pixels > 0) flags <- c(flags, 'MISSING_OR_OUT_OF_RANGE_RGB')
    if (result$excluded_denominator_pixels > 0) flags <- c(flags, 'DENOMINATOR_PIXELS_EXCLUDED')
    if (result$valid_polygon_percent < cfg$min_valid_polygon_percent) flags <- c(flags, 'LOW_VALID_POLYGON_COVERAGE')
    fit <- fit_distribution(df, cfg)
    result$threshold <- fit$threshold; result$otsu_separability <- fit$separability
    result$vari_min <- fit$limits[1]; result$vari_max <- fit$limits[2]
    if (result$vari_min < -1 || result$vari_max > 1) flags <- c(flags, 'VARI_OUTSIDE_MINUS1_PLUS1')
    if (fit$separability < cfg$min_separability) flags <- c(flags, 'LOW_OTSU_SEPARABILITY')
    high <- df$D > cfg$min_denominator_raw & is.finite(df$VARI) & df$VARI > fit$threshold
    result$high_vari_pixels <- sum(df$count[high])
    result$low_vari_pixels <- result$valid_vari_pixels - result$high_vari_pixels
    result$coverage_percent <- 100 * result$high_vari_pixels / result$valid_vari_pixels
    result$high_vari_polygon_percent <- 100 * result$high_vari_pixels / result$polygon_pixels
    result$high_vari_area_m2 <- result$high_vari_pixels * result$pixel_area_m2
    if (result$coverage_percent <= cfg$extreme_cover_percent[1] || result$coverage_percent >= cfg$extreme_cover_percent[2])
      flags <- c(flags, 'EXTREME_CLASS_COVERAGE')
    utils::write.csv(fit$histogram, file.path(audit, paste0(entry$output_key, '_histogram.csv')), row.names = FALSE)
    saveRDS(df, file.path(audit, paste0(entry$output_key, '_joint_counts.rds')))
    message('  Writing binary GeoTIFF and checking class counts')
    tif <- file.path(run_dir, paste0(entry$output_key, '_VARI_Otsu.tif'))
    written <- write_binary_mask(x$rgb, x$inside, fit$threshold, tif, cfg)
    stopifnot(written$counts['high'] == result$high_vari_pixels,
              written$counts['low'] == result$low_vari_pixels)
    result$binary_geotiff <- tif
    png <- file.path(plot_dir, paste0(entry$output_key, '_QC.png'))
    make_preview(x$rgb, x$inside, written$raster, fit, entry, result, png, cfg)
    result$png <- png; result$status <- 'OK_REVIEW_REQUIRED'
  }, error = function(e) {
    result$error <<- conditionMessage(e); flags <<- c(flags, 'PROCESSING_FAILED')
  })
  result$qc_flags <- paste(unique(flags), collapse = ';')
  if (result$status == 'ERROR') {
    message('  ERROR: ', result$error)
    # Clear paths to incomplete artifacts; never present a partial mask as success.
    result$binary_geotiff <- ''
    result$png <- file.path(plot_dir, paste0(entry$output_key, '_ERROR.png'))
    grDevices::png(result$png, width = 1400, height = 800)
    tryCatch({graphics::plot.new(); graphics::title(main = paste(entry$imu_id, 'PROCESSING FAILED'))
      graphics::text(.5, .5, paste(strwrap(result$error, 80), collapse = '\n'))
    }, finally = grDevices::dev.off())
  }
  result
}

run_analysis <- function(cfg = pangaia_config()) {
  require_spatial_packages()
  stopifnot(cfg$expected_imus == 12L, is.finite(cfg$min_denominator_raw), cfg$min_denominator_raw >= 0,
            cfg$histogram_bins >= 2L, cfg$block_rows >= 1L)
  b <- read_boundaries(cfg); m <- image_manifest(b, cfg)
  run <- paste0(format(Sys.time(), '%Y%m%d_%H%M%S'), '_', Sys.getpid())
  run_dir <- file.path(cfg$output_root, run); plot_dir <- file.path(cfg$plots_root, run)
  if (dir.exists(run_dir) || dir.exists(plot_dir)) stop('Run directory already exists; retry in a new second.')
  for (d in c(run_dir, file.path(run_dir, 'audit'), plot_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(m, file.path(run_dir, 'audit', 'image_manifest.csv'), row.names = FALSE)
  dput(cfg, file = file.path(run_dir, 'audit', 'configuration.R'))
  writeLines(capture.output(sessionInfo()), file.path(run_dir, 'audit', 'sessionInfo.txt'))
  results <- list()
  for (i in seq_len(nrow(m))) {
    message('Processing ', i, '/12: ', m$imu_id[i])
    results[[i]] <- process_imu(m[i, ], b[b$imu_id == m$imu_id[i], ], cfg, run_dir, plot_dir)
    utils::write.csv(do.call(rbind, results), file.path(run_dir, 'KA_VARI_Otsu_results.csv'), row.names = FALSE, na = '')
  }
  results <- do.call(rbind, results)
  message('Finished: ', sum(results$status != 'ERROR'), '/12 succeeded; ', sum(results$status == 'ERROR'), ' failed.')
  message('Results: ', file.path(run_dir, 'KA_VARI_Otsu_results.csv'))
  invisible(list(results = results, run_dir = run_dir, plot_dir = plot_dir, config = cfg))
}

# Tests and companion scripts can load functions without processing images.
if (!isTRUE(getOption('pangaia.skip_main'))) run_analysis()
