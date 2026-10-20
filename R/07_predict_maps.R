# =============================================================================
# 07 - Pixel-wise AGB and uncertainty maps
# =============================================================================
# Applies a trained GPR to multispectral stacks clipped to a field boundary.
# Each output is a two-band GeoTIFF: band 1 = DM_mean, band 2 = DM_SE
# (kg DM/ha), LZW-compressed, in the CRS of the input stack.
#
# Input stacks must contain the model's bands in the order listed in
# PREDICTORS (B2 ... B12 for Sentinel-2, PLB2 ... PLB8 for PlanetScope) and the
# same reflectance scaling as the training data (surface reflectance x 10000).
#
# Usage: Rscript R/07_predict_maps.R S2 path/to/stacks path/to/field.shp outputs/maps

source("R/00_utils.R")
library(terra)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript R/07_predict_maps.R <S2|PS> <stack_dir> <field.shp> <out_dir>")
sensor <- args[1]; stack_dir <- args[2]; field_path <- args[3]; out_dir <- args[4]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obj   <- load_model(sensor)
bands <- obj$predictors
field <- vect(field_path)
field_name <- tools::file_path_sans_ext(basename(field_path))

for (f in list.files(stack_dir, pattern = "\\.tiff?$", full.names = TRUE)) {
  r <- rast(f)
  if (nlyr(r) != length(bands)) {
    message("Skipping ", basename(f), ": expected ", length(bands), " bands, found ", nlyr(r))
    next
  }
  names(r) <- bands
  r <- mask(crop(r, project(field, crs(r))), project(field, crs(r)))

  vals  <- as.data.frame(values(r))
  ok    <- which(complete.cases(vals))
  if (length(ok) == 0) { message("No valid pixels in ", basename(f)); next }

  X <- predict(obj$preproc, vals[ok, bands, drop = FALSE])
  p <- predict_gpr(obj$model, X)

  mu <- se <- rep(NA_real_, nrow(vals))
  mu[ok] <- p$mu; se[ok] <- p$se
  out <- c(setValues(r[[1]], mu), setValues(r[[1]], se))
  names(out) <- c("DM_mean", "DM_SE")

  date <- regmatches(basename(f), regexpr("\\d{8}", basename(f)))
  out_file <- file.path(out_dir, sprintf("DM_pred_%s_%s_%s.tif", tolower(sensor), field_name,
                                         if (length(date)) date else tools::file_path_sans_ext(basename(f))))
  writeRaster(out, out_file, overwrite = TRUE, gdal = c("COMPRESS=LZW", "PREDICTOR=2"), NAflag = -9999)
  message("Saved ", out_file, " (", length(ok), " pixels)")
}
