# =============================================================================
# 06 - Residual diagnostics and spatial autocorrelation
# =============================================================================
# * Histogram, residuals vs predicted, normal QQ plot
# * Shapiro-Wilk (normality) and Breusch-Pagan (heteroscedasticity) tests
# * Moran's I on standardised residuals. Repeated surveys of the same point are
#   first collapsed to one value per point with a precision-weighted mean
#   (w = 1 / se^2). kNN weights (k = 8), 9,999 permutations.
#
# Requires: outputs/test_predictions_<sensor>.csv (from 03_train_gpr.R)
# Output  : outputs/residuals_<sensor>.png, outputs/residual_tests_<sensor>.txt
# Usage   : Rscript R/06_residual_diagnostics.R S2   (or PS)

source("R/00_utils.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(lmtest)
  library(sf)
  library(spdep)
})

sensor <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(sensor)) sensor <- "S2"
tag <- tolower(sensor)
col <- SENSOR_COLOURS[sensor]

res <- read.csv(file.path(CONFIG$out_dir, paste0("test_predictions_", tag, ".csv"))) |>
  mutate(resid = observed - predicted, resid_std = resid / se)

# --- Plots ---------------------------------------------------------------------
p_hist <- ggplot(res, aes(resid)) + geom_histogram(bins = 30, fill = col, alpha = 0.8) +
  labs(x = "Residual (kg DM/ha)", y = "Count") + theme_classic()
p_fit  <- ggplot(res, aes(predicted, resid)) + geom_point(colour = col, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Predicted AGB (kg DM/ha)", y = "Residual (kg DM/ha)") + theme_classic()
p_qq   <- ggplot(res, aes(sample = resid)) + stat_qq(colour = col, alpha = 0.7) + stat_qq_line(linetype = "dashed") +
  labs(x = "Theoretical quantiles", y = "Sample quantiles") + theme_classic()

png(file.path(CONFIG$out_dir, paste0("residuals_", tag, ".png")), width = 15, height = 4.5, units = "in", res = 300)
gridExtra_available <- requireNamespace("gridExtra", quietly = TRUE)
if (gridExtra_available) gridExtra::grid.arrange(p_hist, p_fit, p_qq, ncol = 3) else print(p_fit)
dev.off()

# --- Tests ---------------------------------------------------------------------
sink(file.path(CONFIG$out_dir, paste0("residual_tests_", tag, ".txt")), split = TRUE)
cat("Residual diagnostics -", SENSOR_LABELS[sensor], "\n\n")
print(shapiro.test(res$resid))
print(bptest(lm(resid ~ predicted, data = res)))

# --- Moran's I -----------------------------------------------------------------
site_res <- res |>
  mutate(w = 1 / se^2) |>
  group_by(site, point_id) |>
  summarise(lon = mean(lon), lat = mean(lat),
            r_w = weighted.mean(resid_std, w), .groups = "drop")

pts <- st_as_sf(site_res, coords = c("lon", "lat"), crs = 4326) |> st_transform(32632)  # UTM 32N
k   <- min(8L, nrow(site_res) - 1L)
lw  <- nb2listw(knn2nb(knearneigh(st_coordinates(pts), k = k)), style = "W")
set.seed(CONFIG$seed)
moran <- moran.mc(pts$r_w, listw = lw, nsim = 9999)
cat(sprintf("\nMoran's I (point-level, k = %d, n = %d): I = %.4f, permutation p = %.4f\n",
            k, nrow(site_res), moran$statistic, moran$p.value))
sink()
