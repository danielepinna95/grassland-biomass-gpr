# =============================================================================
# 05 - SHAP feature importance
# =============================================================================
# Monte Carlo SHAP values (fastshap) on the test set, using 200 training rows
# as background and 400 simulations per feature. Global importance is the mean
# |SHAP| per band, normalised to percent.
#
# Note: adjacent bands are strongly correlated (r > 0.9 in red-edge/NIR), so
# read importance per spectral region rather than per single band.
#
# Requires: models/gpr_<sensor>.rds
# Output  : outputs/shap_importance_<sensor>.csv, outputs/shap_importance_<sensor>.png
# Usage   : Rscript R/05_shap.R S2   (or PS)

source("R/00_utils.R")
library(fastshap)
library(ggplot2)

sensor <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(sensor)) sensor <- "S2"
tag <- tolower(sensor)

obj   <- load_model(sensor)
preds <- obj$predictors
df    <- load_agb_data(sensor)

test  <- df |> filter(sample_id %in% obj$test_ids)
train <- df |> filter(!sample_id %in% obj$test_ids)

X_test <- as.data.frame(predict(obj$preproc, test[, preds]))
set.seed(CONFIG$seed)
X_bg   <- as.data.frame(predict(obj$preproc, train[sample(nrow(train), min(200, nrow(train))), preds]))

pred_wrapper <- function(object, newdata) predict(object, XX = as.matrix(newdata))

shap <- fastshap::explain(obj$model, X = X_bg, newdata = X_test,
                          pred_wrapper = pred_wrapper, nsim = 400)

imp <- data.frame(band = colnames(shap), mean_abs_shap = colMeans(abs(as.matrix(shap)))) |>
  mutate(share_pct = 100 * mean_abs_shap / sum(mean_abs_shap)) |>
  arrange(desc(share_pct))
print(imp, digits = 3)
write.csv(imp, file.path(CONFIG$out_dir, paste0("shap_importance_", tag, ".csv")), row.names = FALSE)

fig <- ggplot(imp, aes(reorder(band, share_pct), share_pct)) +
  geom_col(fill = SENSOR_COLOURS[sensor]) +
  coord_flip() +
  labs(x = NULL, y = "Share of mean |SHAP| (%)", title = paste(SENSOR_LABELS[sensor], "- SHAP importance")) +
  theme_classic(base_size = 13)
ggsave(file.path(CONFIG$out_dir, paste0("shap_importance_", tag, ".png")), fig, width = 6, height = 4.5, dpi = 300)
