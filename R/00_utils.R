# =============================================================================
# Shared configuration and helper functions
# =============================================================================
# Sourced by every script in R/. Edit the CONFIG block to change sensor,
# paths or random seed; everything else reads from here.

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
})

CONFIG <- list(
  data_path  = "data/grassland_agb_s2_planetscope.csv",
  model_dir  = "models",
  out_dir    = "outputs",
  target     = "agb_kg_ha",
  seed       = 42,
  train_frac = 0.8,
  n_strata   = 5,        # quantile classes of AGB used for stratified splits
  nugget_min = 1e-6,     # selected from a 1e-8 ... 1e-4 grid (see paper)
  n_repeats  = 20,       # repeated stratified 80/20 splits
  n_boot     = 10000     # bootstrap iterations for metric CIs
)

PREDICTORS <- list(
  S2 = c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12"),
  PS = c("PLB2", "PLB3", "PLB4", "PLB5", "PLB6", "PLB7", "PLB8")
)

SENSOR_COLOURS <- c(S2 = "brown4", PS = "cadetblue4")
SENSOR_LABELS  <- c(S2 = "Sentinel-2", PS = "PlanetScope")

dir.create(CONFIG$model_dir, showWarnings = FALSE)
dir.create(CONFIG$out_dir,   showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------

#' Load the field dataset and keep rows usable for a given sensor.
#'
#' @param sensor "S2" or "PS".
#' @param trim   If TRUE, drop extreme AGB values with an IQR-like rule on the
#'               5th-95th percentiles (as described in the paper).
load_agb_data <- function(sensor, path = CONFIG$data_path, trim = TRUE) {
  preds <- PREDICTORS[[sensor]]
  y     <- CONFIG$target

  df <- read.csv(path, check.names = FALSE) |>
    mutate(across(all_of(y), as.numeric)) |>
    filter(!is.na(.data[[y]]), .data[[y]] > 0) |>
    filter(if_all(all_of(preds), ~ !is.na(.x)))

  if (trim) {
    q   <- quantile(df[[y]], c(0.05, 0.95))
    iqr <- q[2] - q[1]
    df  <- df |> filter(.data[[y]] > q[1] - 1.5 * iqr, .data[[y]] < q[2] + 1.5 * iqr)
  }

  df |> mutate(agb_class = cut(.data[[y]],
                               breaks = quantile(.data[[y]], probs = seq(0, 1, length.out = CONFIG$n_strata + 1)),
                               include.lowest = TRUE, labels = FALSE))
}

#' Stratified train/test split on AGB quantile classes.
stratified_split <- function(df, p = CONFIG$train_frac) {
  idx <- createDataPartition(df$agb_class, p = p, list = FALSE)
  list(train = df[idx, ], test = df[-idx, ])
}

#' Centre and scale predictors using training-set statistics only.
scale_split <- function(split, predictors) {
  pre <- preProcess(split$train[, predictors], method = c("center", "scale"))
  list(
    preproc = pre,
    X_train = as.matrix(predict(pre, split$train[, predictors])),
    X_test  = as.matrix(predict(pre, split$test[,  predictors])),
    y_train = split$train[[CONFIG$target]],
    y_test  = split$test[[CONFIG$target]]
  )
}

# -----------------------------------------------------------------------------
# GPR
# -----------------------------------------------------------------------------

#' Fit the final GPR: additive Matern 5/2 + Exponential kernel (GauPro).
fit_gpr <- function(X, y, nug_min = CONFIG$nugget_min) {
  D <- ncol(X)
  kernel <- GauPro::k_Matern52(D = D) + GauPro::k_Exponential(D = D)
  GauPro::gpkm(X = X, Z = y, kernel = kernel, parallel = FALSE, nug.min = nug_min)
}

#' Predict mean and posterior standard error.
predict_gpr <- function(model, X) {
  out <- predict(model, XX = as.matrix(X), se.fit = TRUE)
  data.frame(mu = as.numeric(out$mean), se = pmax(as.numeric(out$se), 1e-9))
}

# -----------------------------------------------------------------------------
# Metrics
# -----------------------------------------------------------------------------

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
mae  <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)
r2   <- function(y, yhat) 1 - sum((y - yhat)^2, na.rm = TRUE) / sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)

metrics <- function(y, yhat) c(R2 = r2(y, yhat), RMSE = rmse(y, yhat), MAE = mae(y, yhat))

#' Bootstrap confidence intervals for R2, RMSE and MAE.
bootstrap_metrics <- function(y, yhat, n_boot = CONFIG$n_boot, conf = 0.95, seed = CONFIG$seed) {
  set.seed(seed)
  boot <- t(replicate(n_boot, {
    i <- sample(seq_along(y), replace = TRUE)
    metrics(y[i], yhat[i])
  }))
  a <- (1 - conf) / 2
  data.frame(
    metric   = colnames(boot),
    estimate = metrics(y, yhat),
    ci_low   = apply(boot, 2, quantile, a),
    ci_high  = apply(boot, 2, quantile, 1 - a),
    row.names = NULL
  )
}

#' Paths for saved model objects.
model_path <- function(sensor) file.path(CONFIG$model_dir, paste0("gpr_", tolower(sensor), ".rds"))

save_model <- function(sensor, model, preproc, split_ids) {
  saveRDS(list(model = model, preproc = preproc, predictors = PREDICTORS[[sensor]],
               test_ids = split_ids, created = Sys.time()),
          model_path(sensor))
}

load_model <- function(sensor) readRDS(model_path(sensor))
