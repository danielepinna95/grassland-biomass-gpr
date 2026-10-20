# Grassland above-ground biomass from Sentinel-2 and PlanetScope with Gaussian Process Regression

Code for:

> Pinna, D., Basso, E., Pornaro, C., Pullanagari, R., Macolino, S., Pezzuolo, A., & Marinello, F. (2025).
> **Optimising grassland above-ground biomass estimation for managed grasslands: a Gaussian process regression approach for Sentinel-2 and PlanetScope in Northern Italy.**
> *Precision Agriculture*, 26, 98. https://doi.org/10.1007/s11119-025-10298-2

The repository trains Gaussian Process Regression (GPR) models that predict grassland above-ground biomass (AGB, kg DM/ha) from multispectral reflectance. It reports calibrated uncertainty for every prediction and turns satellite scenes into biomass and uncertainty maps.

| | Sentinel-2 | PlanetScope |
|---|---|---|
| Pixel level, test set | R² 0.52, MAE ≈ 400 kg DM/ha | R² 0.51, MAE ≈ 400 kg DM/ha |
| Field level (site × survey) | R² 0.97, MAE 60–120 kg DM/ha | R² 0.97, MAE 60–120 kg DM/ha |

*Values from the paper.*

## What's inside

```
data/      expected input format (the field dataset is not public yet, see data/README.md)
python/    extract_reflectance.py   – sample raster stacks at survey points
R/
  00_utils.R                   – configuration, data loading, splits, GPR and metric helpers
  01_benchmark_models.R        – RF, SVR, MLP, Lasso and GPR on the same split
  02_kernel_selection.R        – 10 kernels + nugget grid, scored by CV on the training set
  03_train_gpr.R               – final Matern 5/2 + Exponential GPR, test metrics with bootstrap CIs
  04_robustness_uncertainty.R  – 20 repeated splits, conformal-calibrated intervals, field-level accuracy
  05_shap.R                    – SHAP band importance (fastshap)
  06_residual_diagnostics.R    – Shapiro-Wilk, Breusch-Pagan, QQ plot, Moran's I
  07_predict_maps.R            – pixel-wise AGB mean + SE GeoTIFFs for a field
run_all.R  – runs 01–06 for both sensors
```

## Method in brief

1. **Field data.** 950 cut-and-dry AGB samples collected over 18 months (May 2023 – Nov 2024) in three grassland systems near Schio and on the Asiago plateau (Vicenza, Italy): a lowland meadow, a lowland pasture and an alpine pasture (Malga Carriola).
2. **Spectra.** Surface reflectance from the closest cloud-free Sentinel-2 L2A and PlanetScope SuperDove scene, within 10 days of each survey. Coastal, water-vapour and cirrus bands are excluded.
3. **Pre-processing.** Extreme AGB values are trimmed with an IQR-like rule on the 5th–95th percentiles. The data is split 80/20, stratified on AGB quintiles, and predictors are centred and scaled on the training set.
4. **Model.** GPR (GauPro) with an additive Matern 5/2 + Exponential kernel and a minimum nugget of 1e-6, both selected by cross-validation.
5. **Evaluation.** Bootstrap CIs (10,000 resamples), 20 repeated stratified splits, cross-conformal calibration of the posterior SE, aggregation to field level, SHAP importance and residual spatial autocorrelation.
6. **Deployment.** Per-pixel mean and SE maps for any field polygon.

## Quick start

```r
# R >= 4.2
install.packages(c("GauPro", "caret", "dplyr", "ggplot2", "fastshap", "lmtest",
                   "sf", "spdep", "terra", "gridExtra",
                   "ranger", "kernlab", "nnet", "glmnet"))   # last four: benchmark only
```

```bash
# from the repository root
Rscript R/03_train_gpr.R S2                 # train + evaluate the Sentinel-2 model
Rscript R/04_robustness_uncertainty.R S2    # robustness, calibrated intervals, field level
Rscript run_all.R                           # everything, both sensors (slow)
```

Use `PS` instead of `S2` for PlanetScope. Outputs go to `outputs/` and fitted models to `models/`.

The scripts expect the field dataset at `data/grassland_agb_s2_planetscope.csv`. It is not included yet, so use your own data in the format described in [`data/README.md`](data/README.md).

### Mapping your own fields

```bash
# 1. sample reflectance at your survey points (optional, for new training data)
pip install geopandas rasterio pandas
python python/extract_reflectance.py points.csv scenes/ samples.csv --bands B2 B3 B4 B5 B6 B7 B8 B8A B11 B12

# 2. predict AGB + uncertainty maps with a trained model
Rscript R/07_predict_maps.R S2 scenes/ my_field.shp outputs/maps
```

Input stacks must have the model's bands in order (`B2 … B12` for Sentinel-2, `PLB2 … PLB8` for PlanetScope). Reflectance must be scaled ×10,000, as in the training data. The models were trained on grasslands in north-east Italy. Before using them elsewhere, check them against local field data.

## Data

The field dataset (950 samples) is not publicly available yet. For access, contact the corresponding author. [`data/README.md`](data/README.md) describes the input format, so the code can be run on other datasets.

## Citation

If you use this code, please cite the paper above (see also [`CITATION.cff`](CITATION.cff)).

## License

Code: MIT (see [`LICENSE`](LICENSE)).

## Acknowledgements

This study was carried out within the Agritech National Research Center and received funding from the European Union Next-GenerationEU (PNRR – Mission 4, Component 2, Investment 1.4 – D.D. 1032 17/06/2022, CN00000022). University of Padova (TESAF and DAFNAE departments).
