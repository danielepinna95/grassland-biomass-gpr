# =============================================================================
# 04 - Robustness, calibrated uncertainty and field-scale accuracy
# =============================================================================
# 1. Robustness: refit the GPR on 20 repeated stratified 80/20 splits.
# 2. Calibration: cross-conformal scaling of the GPR posterior SE. The 95th
#    quantile of |y - mu| / se over the out-of-fold residuals of step 1
#    rescales the intervals: PI95 = mu +/- q95 * se.
# 3. Field scale: aggregate test-set predictions per site and survey.
# 4. Accuracy by AGB class.
#
# Requires: models/gpr_<sensor>.rds (from 03_train_gpr.R)
# Output  : outputs/robustness_<sensor>.csv, outputs/uncertainty_<sensor>.csv,
#           outputs/field_level_<sensor>.csv, outputs/metrics_by_bin_<sensor>.csv
# Usage   : Rscript R/04_robustness_uncertainty.R S2   (or PS)

source("R/00_utils.R")

sensor <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(sensor)) sensor <- "S2"
preds <- PREDICTORS[[sensor]]
tag   <- tolower(sensor)

df <- load_agb_data(sensor)

# --- 1. Repeated splits ------------------------------------------------------
set.seed(CONFIG$seed)
reps <- lapply(seq_len(CONFIG$n_repeats), function(r) {
  cat("Repeat", r, "/", CONFIG$n_repeats, "\n")
  s  <- scale_split(stratified_split(df), preds)
  gp <- fit_gpr(s$X_train, s$y_train)
  p  <- predict_gpr(gp, s$X_test)
  list(metrics = data.frame(rep = r, t(metrics(s$y_test, p$mu))),
       oof     = data.frame(rep = r, y = s$y_test, mu = p$mu, se = p$se))
})
rep_metrics <- bind_rows(lapply(reps, `[[`, "metrics"))
robust <- rep_metrics |>
  summarise(across(c(R2, RMSE, MAE), list(mean = mean, sd = sd,
                                         q05 = ~ quantile(.x, 0.05), q95 = ~ quantile(.x, 0.95))))
print(round(robust, 3))
write.csv(rep_metrics, file.path(CONFIG$out_dir, paste0("robustness_", tag, ".csv")), row.names = FALSE)

# --- 2. Conformal calibration -------------------------------------------------
oof   <- bind_rows(lapply(reps, `[[`, "oof"))
score <- sort(abs(oof$y - oof$mu) / oof$se)
k     <- ceiling((length(score) + 1) * 0.95)
q95   <- score[min(k, length(score))]
cat(sprintf("Conformal scale q0.95 = %.3f\n", q95))

test <- read.csv(file.path(CONFIG$out_dir, paste0("test_predictions_", tag, ".csv")))

interval_summary <- function(y, mu, se, q) {
  lo <- mu - q * se; hi <- mu + q * se
  data.frame(coverage95      = mean(y >= lo & y <= hi),
             median_width    = median(hi - lo),
             median_rel_width = median((hi - lo) / pmax(abs(mu), 1e-6)),
             median_cv_pct   = median(100 * se / pmax(abs(mu), 1e-6)))
}

pixel_unc <- cbind(scale = "pixel", q95 = q95,
                   interval_summary(test$observed, test$predicted, test$se, q95))

# --- 3. Field scale (site x survey) -------------------------------------------
field <- test |>
  group_by(site, survey) |>
  summarise(n_points  = n(),
            observed  = mean(observed),
            predicted = mean(predicted),
            se        = sqrt(sum(se^2)) / n(),   # SE of the mean, assuming independence
            .groups   = "drop")
write.csv(field, file.path(CONFIG$out_dir, paste0("field_level_", tag, ".csv")), row.names = FALSE)

field_unc <- cbind(scale = "field", q95 = q95,
                   interval_summary(field$observed, field$predicted, field$se, q95))
unc <- rbind(pixel_unc, field_unc)
print(unc, digits = 3)
write.csv(unc, file.path(CONFIG$out_dir, paste0("uncertainty_", tag, ".csv")), row.names = FALSE)

cat("\nField-level accuracy:\n")
print(bootstrap_metrics(field$observed, field$predicted), digits = 3)

# --- 4. Accuracy by AGB class ------------------------------------------------
bins <- test |>
  mutate(bin = cut(observed, c(0, 500, 1000, 1500, 2000, Inf),
                   labels = c("<500", "500-1000", "1000-1500", "1500-2000", ">2000"),
                   include.lowest = TRUE)) |>
  group_by(bin) |>
  summarise(n = n(), RMSE = rmse(observed, predicted), MAE = mae(observed, predicted),
            bias = mean(predicted - observed), .groups = "drop")
print(bins, digits = 3)
write.csv(bins, file.path(CONFIG$out_dir, paste0("metrics_by_bin_", tag, ".csv")), row.names = FALSE)
