# Reproduce the full analysis for both sensors.
# Run from the repository root:  Rscript run_all.R
#
# Fitting GPRs is the slow part (20 repeats x 2 sensors in step 04, and
# 5-fold CV x 10 kernels in step 02); expect tens of minutes on a laptop.

steps <- c("01_benchmark_models.R", "02_kernel_selection.R", "03_train_gpr.R",
           "04_robustness_uncertainty.R", "05_shap.R", "06_residual_diagnostics.R")

for (sensor in c("S2", "PS")) {
  for (step in steps) {
    message("\n========== ", step, " [", sensor, "] ==========")
    status <- system2("Rscript", c(file.path("R", step), sensor))
    if (status != 0) stop(step, " failed for ", sensor)
  }
}
