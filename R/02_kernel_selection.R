# =============================================================================
# 02 - GPR kernel and nugget selection
# =============================================================================
# Every candidate kernel is scored by 5-fold cross-validation on the TRAINING
# set, so the test set never influences the choice. The best kernel is then
# re-scored over a grid of minimum nugget values.
#
# Output: outputs/kernel_selection_<sensor>.csv, outputs/nugget_selection_<sensor>.csv
# Usage : Rscript R/02_kernel_selection.R S2   (or PS)

source("R/00_utils.R")
library(GauPro)

sensor <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(sensor)) sensor <- "S2"
preds <- PREDICTORS[[sensor]]

set.seed(CONFIG$seed)
df    <- load_agb_data(sensor)
split <- stratified_split(df)
train <- split$train

kernels <- list(
  Matern52        = function(D) k_Matern52(D = D),
  Matern32        = function(D) k_Matern32(D = D),
  Exponential     = function(D) k_Exponential(D = D),
  PowerExp        = function(D) k_PowerExp(D = D),
  RationalQuad    = function(D) k_RatQuad(D = D),
  Triangle        = function(D) k_Triangle(D = D),
  Cubic           = function(D) k_Cubic(D = D),
  Gaussian        = function(D) k_Gaussian(D = D),
  `Matern52+Exp`  = function(D) k_Matern52(D = D) + k_Exponential(D = D),
  `Matern32+Exp`  = function(D) k_Matern32(D = D) + k_Exponential(D = D)
)

folds <- createFolds(train$agb_class, k = 5)

cv_score <- function(make_kernel, nug_min = CONFIG$nugget_min) {
  per_fold <- lapply(folds, function(val_idx) {
    tr  <- train[-val_idx, ]; va <- train[val_idx, ]
    pre <- preProcess(tr[, preds], method = c("center", "scale"))
    Xtr <- as.matrix(predict(pre, tr[, preds])); Xva <- as.matrix(predict(pre, va[, preds]))
    gp  <- gpkm(X = Xtr, Z = tr[[CONFIG$target]], kernel = make_kernel(ncol(Xtr)),
                parallel = FALSE, nug.min = nug_min)
    metrics(va[[CONFIG$target]], predict(gp, XX = Xva))
  })
  colMeans(do.call(rbind, per_fold))
}

kernel_res <- bind_rows(lapply(names(kernels), function(k) {
  cat("Kernel:", k, "\n")
  data.frame(kernel = k, t(cv_score(kernels[[k]])))
})) |> arrange(RMSE)
print(kernel_res, digits = 3)
write.csv(kernel_res, file.path(CONFIG$out_dir, paste0("kernel_selection_", tolower(sensor), ".csv")), row.names = FALSE)

best <- kernel_res$kernel[1]
nugget_res <- bind_rows(lapply(10^(-8:-4), function(n) {
  data.frame(kernel = best, nug_min = n, t(cv_score(kernels[[best]], nug_min = n)))
}))
print(nugget_res, digits = 3)
write.csv(nugget_res, file.path(CONFIG$out_dir, paste0("nugget_selection_", tolower(sensor), ".csv")), row.names = FALSE)
