# =============================================================================
# 01 - Benchmark of regression algorithms
# =============================================================================
# Compares Random Forest, Support Vector Regression, a multilayer perceptron,
# Lasso and Gaussian Process Regression on the same stratified 80/20 split.
# Hyperparameters are tuned by 5-fold cross-validation on the training set
# only; the test set is touched once, at the end.
#
# Output: outputs/benchmark_<sensor>.csv
# Usage : Rscript R/01_benchmark_models.R S2   (or PS)

source("R/00_utils.R")

sensor <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(sensor)) sensor <- "S2"
preds <- PREDICTORS[[sensor]]

set.seed(CONFIG$seed)
df    <- load_agb_data(sensor)
split <- stratified_split(df)
s     <- scale_split(split, preds)

ctrl <- trainControl(method = "cv", number = 5)

candidates <- list(
  RandomForest = list(method = "ranger",       tuneLength = 5),
  SVR          = list(method = "svmRadial",    tuneLength = 8),
  MLP          = list(method = "nnet",         tuneGrid = expand.grid(size = c(3, 5, 10), decay = c(1e-3, 1e-2, 1e-1)),
                      linout = TRUE, trace = FALSE, maxit = 1000),
  Lasso        = list(method = "glmnet",       tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-3, 3, length.out = 30))),
  GPR_RBF      = list(method = "gaussprRadial", tuneLength = 5)
)

results <- lapply(names(candidates), function(name) {
  args <- candidates[[name]]
  fit  <- do.call(caret::train, c(list(x = s$X_train, y = s$y_train, trControl = ctrl), args))
  yhat <- predict(fit, s$X_test)
  data.frame(model = name, t(metrics(s$y_test, yhat)))
})

# The final GPR (Matern 5/2 + Exponential, GauPro) for reference
if (requireNamespace("GauPro", quietly = TRUE)) {
  gp <- fit_gpr(s$X_train, s$y_train)
  results[[length(results) + 1]] <- data.frame(model = "GPR_Matern52+Exp",
                                               t(metrics(s$y_test, predict_gpr(gp, s$X_test)$mu)))
}

results <- bind_rows(results) |> arrange(desc(R2))
print(results, digits = 3)
write.csv(results, file.path(CONFIG$out_dir, paste0("benchmark_", tolower(sensor), ".csv")), row.names = FALSE)
