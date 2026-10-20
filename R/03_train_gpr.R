# =============================================================================
# 03 - Train and evaluate the final GPR models
# =============================================================================
# Fits the Matern 5/2 + Exponential GPR for one sensor on a stratified 80/20
# split, reports test-set metrics with bootstrap 95% CIs and saves the model,
# the scaling object and the test-set predictions.
#
# Output: models/gpr_<sensor>.rds
#         outputs/test_predictions_<sensor>.csv
#         outputs/test_metrics_<sensor>.csv
#         outputs/observed_vs_predicted_<sensor>.png
# Usage : Rscript R/03_train_gpr.R S2   (or PS)

source("R/00_utils.R")
library(ggplot2)

sensor <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(sensor)) sensor <- "S2"
preds <- PREDICTORS[[sensor]]

set.seed(CONFIG$seed)
df    <- load_agb_data(sensor)
split <- stratified_split(df)
s     <- scale_split(split, preds)

cat("Training", SENSOR_LABELS[sensor], "GPR on", nrow(s$X_train), "samples...\n")
gp <- fit_gpr(s$X_train, s$y_train)
save_model(sensor, gp, s$preproc, split$test$sample_id)

p <- predict_gpr(gp, s$X_test)
test_pred <- split$test |>
  select(sample_id, point_id, site, survey, lon, lat, observed = all_of(CONFIG$target)) |>
  mutate(predicted = p$mu, se = p$se)
write.csv(test_pred, file.path(CONFIG$out_dir, paste0("test_predictions_", tolower(sensor), ".csv")), row.names = FALSE)

m <- bootstrap_metrics(test_pred$observed, test_pred$predicted)
print(m, digits = 3)
write.csv(m, file.path(CONFIG$out_dir, paste0("test_metrics_", tolower(sensor), ".csv")), row.names = FALSE)

lim <- c(0, max(c(test_pred$observed, test_pred$predicted)) * 1.05)
fig <- ggplot(test_pred, aes(observed, predicted)) +
  geom_abline(linetype = "dashed", colour = "grey40") +
  geom_point(colour = SENSOR_COLOURS[sensor], alpha = 0.7, size = 2.5) +
  annotate("text", x = lim[2] * 0.03, y = lim[2] * 0.97, hjust = 0, vjust = 1, size = 5,
           label = sprintf("R² = %.3f\nRMSE = %.0f kg DM/ha\nMAE = %.0f kg DM/ha",
                           m$estimate[1], m$estimate[2], m$estimate[3])) +
  coord_fixed(xlim = lim, ylim = lim, expand = FALSE) +
  labs(x = "Observed AGB (kg DM/ha)", y = "Predicted AGB (kg DM/ha)",
       title = paste(SENSOR_LABELS[sensor], "- pixel level, test set")) +
  theme_classic(base_size = 13)
ggsave(file.path(CONFIG$out_dir, paste0("observed_vs_predicted_", tolower(sensor), ".png")),
       fig, width = 6, height = 6, dpi = 300)
