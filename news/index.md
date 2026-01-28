# Changelog

## cfperformance 0.5.0

Adds factual (non-counterfactual) prediction model transportability,
configurable propensity score trimming, improved standard error methods,
and critical bug fixes.

### New Features

#### Factual Prediction Model Transportability

- All transportability functions
  ([`tr_mse()`](https://boyercb.github.io/cfperformance/reference/tr_mse.md),
  [`tr_auc()`](https://boyercb.github.io/cfperformance/reference/tr_auc.md),
  [`tr_sensitivity()`](https://boyercb.github.io/cfperformance/reference/tr_sensitivity.md),
  [`tr_specificity()`](https://boyercb.github.io/cfperformance/reference/tr_specificity.md),
  [`tr_fpr()`](https://boyercb.github.io/cfperformance/reference/tr_fpr.md),
  [`tr_roc()`](https://boyercb.github.io/cfperformance/reference/tr_roc.md),
  [`tr_calibration()`](https://boyercb.github.io/cfperformance/reference/tr_calibration.md))
  now support **factual mode** when `treatment = NULL`.
- Factual mode estimates prediction model performance in the target
  population for observed (factual) outcomes, without requiring
  treatment/intervention data.
- This enables standard prediction model transportability analysis
  (covariate shift correction) without counterfactual assumptions.
- Print/summary methods display “Factual” or “Counterfactual” mode
  labels.

#### Propensity Score Trimming

- New `ps_trim` parameter for all `cf_*` and `tr_*` functions provides
  configurable propensity score trimming:
  - `NULL` (default): Absolute bounds `c(0.01, 0.99)`
  - `"none"`: No trimming
  - `"quantile"`: Quantile-based trimming
  - `"absolute"`: Explicit absolute bounds
  - Numeric vector `c(lower, upper)` for custom bounds
  - List with `method` and `bounds` for full control

#### Influence Function Standard Errors

- Exposed influence function SE via `se_method = "influence"` for:
  - [`cf_sensitivity()`](https://boyercb.github.io/cfperformance/reference/cf_sensitivity.md),
    [`cf_specificity()`](https://boyercb.github.io/cfperformance/reference/cf_specificity.md),
    [`cf_fpr()`](https://boyercb.github.io/cfperformance/reference/cf_fpr.md)
  - [`tr_sensitivity()`](https://boyercb.github.io/cfperformance/reference/tr_sensitivity.md),
    [`tr_specificity()`](https://boyercb.github.io/cfperformance/reference/tr_specificity.md),
    [`tr_fpr()`](https://boyercb.github.io/cfperformance/reference/tr_fpr.md)
- Requires `cross_fit = TRUE` for valid inference.

#### Outcome Type Detection

- Added `outcome_type` parameter to transportability MSE functions.
- Auto-detection of binary vs continuous outcomes.
- Binary outcomes use efficient transformation; continuous outcomes
  model loss directly.

#### Cross-Fitting for Transportability

- Extended cross-fitting support to transportability estimators for
  valid inference with flexible ML methods.

### Bug Fixes

- **Fixed bootstrap confidence interval coverage** for all estimators.
  Bootstrap now correctly preserves user-specified model formulas (e.g.,
  with quadratic terms like `I(X^2)`) when refitting models during
  resampling. Previously, bootstrap used `Y ~ .` which only fit main
  effects, causing severely undercovered confidence intervals when
  models included non-linear terms.

### Documentation

- Renamed “Traditional” mode to “Factual” mode throughout documentation
  and print output. This better reflects the distinction between factual
  (observed) outcomes and counterfactual outcomes under hypothetical
  interventions.
- Added `boot_ci_type` parameter documentation for bootstrap CI method
  selection.
- Added comprehensive simulation study script for benchmarking.
- Added benchmark tests verifying estimators against standard
  implementations (WeightedROC, pROC).

------------------------------------------------------------------------

## cfperformance 0.4.0

Adds sensitivity, specificity, and ROC curve functions for both
counterfactual and transportability settings.

### New Features

#### Counterfactual Sensitivity/Specificity

- [`cf_sensitivity()`](https://boyercb.github.io/cfperformance/reference/cf_sensitivity.md) -
  Counterfactual sensitivity (true positive rate)
- [`cf_specificity()`](https://boyercb.github.io/cfperformance/reference/cf_specificity.md) -
  Counterfactual specificity (true negative rate)
- [`cf_fpr()`](https://boyercb.github.io/cfperformance/reference/cf_fpr.md) -
  Counterfactual false positive rate (1 - specificity)
- [`cf_tpr()`](https://boyercb.github.io/cfperformance/reference/cf_sensitivity.md) -
  Alias for
  [`cf_sensitivity()`](https://boyercb.github.io/cfperformance/reference/cf_sensitivity.md)
- [`cf_tnr()`](https://boyercb.github.io/cfperformance/reference/cf_specificity.md) -
  Alias for
  [`cf_specificity()`](https://boyercb.github.io/cfperformance/reference/cf_specificity.md)
- Supports CL, IPW, DR, and naive estimators
- Vectorized threshold parameter for efficient ROC curve computation

#### Transportable Sensitivity/Specificity

- [`tr_sensitivity()`](https://boyercb.github.io/cfperformance/reference/tr_sensitivity.md) -
  Transportable sensitivity for target population
- [`tr_specificity()`](https://boyercb.github.io/cfperformance/reference/tr_specificity.md) -
  Transportable specificity for target population
- [`tr_fpr()`](https://boyercb.github.io/cfperformance/reference/tr_fpr.md) -
  Transportable false positive rate
- [`tr_tpr()`](https://boyercb.github.io/cfperformance/reference/tr_sensitivity.md) -
  Alias for
  [`tr_sensitivity()`](https://boyercb.github.io/cfperformance/reference/tr_sensitivity.md)
- [`tr_tnr()`](https://boyercb.github.io/cfperformance/reference/tr_specificity.md) -
  Alias for
  [`tr_specificity()`](https://boyercb.github.io/cfperformance/reference/tr_specificity.md)
- Supports OM, IPW, DR, and naive estimators
- Works with both “transport” and “joint” analysis types

#### ROC Curves

- [`tr_roc()`](https://boyercb.github.io/cfperformance/reference/tr_roc.md) -
  Compute transportable ROC curve in target population
- [`cf_roc()`](https://boyercb.github.io/cfperformance/reference/cf_roc.md) -
  Compute counterfactual ROC curve
- [`plot.tr_roc()`](https://boyercb.github.io/cfperformance/reference/plot.tr_roc.md)
  /
  [`plot.cf_roc()`](https://boyercb.github.io/cfperformance/reference/plot.tr_roc.md) -
  Plot ROC curves with AUC in legend
- [`as.data.frame.tr_roc()`](https://boyercb.github.io/cfperformance/reference/as.data.frame.tr_roc.md)
  /
  [`as.data.frame.cf_roc()`](https://boyercb.github.io/cfperformance/reference/as.data.frame.tr_roc.md) -
  Convert to data frame for ggplot2
- AUC computed via trapezoidal integration
- Option to include naive ROC curve for comparison

#### Documentation

- Updated introduction vignette with ROC curve examples
- Updated transportability vignette with ROC curve examples

------------------------------------------------------------------------

## cfperformance 0.3.0

Adds machine learning integration for flexible nuisance model estimation
with automatic cross-fitting for valid inference.

### New Features

#### ML Learner Interface

- [`ml_learner()`](https://boyercb.github.io/cfperformance/reference/ml_learner.md) -
  Specify ML methods for propensity score and outcome models
- Supports: `ranger`, `xgboost`, `grf`, `glmnet`, `superlearner`, and
  `custom`
- Automatic cross-fitting when `ml_learner` specs are detected
- Seamlessly integrates with existing `propensity_model`/`outcome_model`
  arguments

#### Cross-Fitting Support

- [`cf_mse()`](https://boyercb.github.io/cfperformance/reference/cf_mse.md) -
  Full support for ML learners with cross-fitting
- [`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md) -
  Full support for ML learners with cross-fitting (DR estimator)

#### Supported Learners

- **ranger** - Fast random forest implementation
- **xgboost** - Gradient boosting (XGBoost)
- **grf** - Generalized random forests with honest estimation
- **glmnet** - Elastic net regularization with CV-selected λ
- **superlearner** - Ensemble learning
- **custom** - User-supplied fit/predict functions

#### Usage Example

``` r
# MSE with ML learners
cf_mse(
  predictions = pred, outcomes = y, treatment = a, covariates = df,
  propensity_model = ml_learner("ranger", num.trees = 500),
  outcome_model = ml_learner("xgboost", nrounds = 100),
  cross_fit = TRUE
)

# AUC with ML learners
cf_auc(
  predictions = pred, outcomes = y, treatment = a, covariates = df,
  propensity_model = ml_learner("ranger", num.trees = 500),
  outcome_model = ml_learner("ranger", num.trees = 500),
  cross_fit = TRUE
)
```

#### Documentation

- New vignette: “Machine Learning Integration”
- Updated README with ML integration examples

### References

Chernozhukov, V., et al. (2018). Double/debiased machine learning for
treatment and structural parameters. *The Econometrics Journal*, 21(1),
C1-C68.

Li, B., Gatsonis, C., Dahabreh, I. J., & Steingrimsson, J. A. (2022).
Estimating the area under the ROC curve when transporting a prediction
model to a target population. *Biometrics*, 79(3), 2343-2356.

------------------------------------------------------------------------

## cfperformance 0.2.0

Major release adding transportability estimators from Voter et
al. (2025) for evaluating prediction model performance when transporting
from a source population (e.g., RCT) to a target population.

### New Features

#### Transportability Functions

- [`tr_mse()`](https://boyercb.github.io/cfperformance/reference/tr_mse.md) -
  Transportable MSE estimation with naive, om, ipw, dr estimators
- [`tr_auc()`](https://boyercb.github.io/cfperformance/reference/tr_auc.md) -
  Transportable AUC estimation with all estimators
- [`tr_calibration()`](https://boyercb.github.io/cfperformance/reference/tr_calibration.md) -
  Transportable calibration curves with ICI, E50, E90, Emax

#### Analysis Modes

- **Transport analysis** (`analysis = "transport"`): Use source/RCT
  outcomes to estimate performance in target population
- **Joint analysis** (`analysis = "joint"`): Pool source and target data
  for potentially more efficient estimation

#### Inference

- Bootstrap standard errors with stratified sampling option
  (`stratified_boot = TRUE`) to preserve source/target ratio
- Influence function-based standard errors for
  [`tr_mse()`](https://boyercb.github.io/cfperformance/reference/tr_mse.md)
  (all estimators)

#### New Data

- `transport_sim` - Simulated dataset with source (RCT) and target
  populations for transportability examples (n = 1,500)

#### Documentation

- New vignette: “Transportability Analysis with cfperformance”
- Complete roxygen documentation for all new functions
- Updated README with transportability examples

### S3 Methods for tr\_\* Functions

- [`print.tr_performance()`](https://boyercb.github.io/cfperformance/reference/print.tr_performance.md) -
  Print method for transportability results
- [`summary.tr_performance()`](https://boyercb.github.io/cfperformance/reference/summary.tr_performance.md) -
  Detailed summary
- [`coef.tr_performance()`](https://boyercb.github.io/cfperformance/reference/coef.tr_performance.md) -
  Extract point estimates
- [`confint.tr_performance()`](https://boyercb.github.io/cfperformance/reference/confint.tr_performance.md) -
  Confidence intervals
- [`plot.tr_calibration()`](https://boyercb.github.io/cfperformance/reference/plot.tr_calibration.md) -
  Calibration curve visualization

### References

Voter SR, et al. Transportability of machine learning-based
counterfactual prediction models with application to CASS. *Diagnostic
and Prognostic Research*. 2025; 9(4). <doi:10.1186/s41512-025-00201-y>

------------------------------------------------------------------------

## cfperformance 0.1.0

Initial release implementing methods from Boyer, Dahabreh &
Steingrimsson (2025), “Estimating and evaluating counterfactual
prediction models.”

### Features

#### Performance Metrics

- [`cf_mse()`](https://boyercb.github.io/cfperformance/reference/cf_mse.md) -
  Counterfactual MSE/Brier score estimation
- [`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md) -
  Counterfactual AUC estimation  
- [`cf_calibration()`](https://boyercb.github.io/cfperformance/reference/cf_calibration.md) -
  Counterfactual calibration curves with ICI, E50, E90, Emax

#### Estimators

- Naive estimator (subset-based)
- Conditional Loss (CL) / Outcome modeling estimator
- Inverse Probability Weighting (IPW) estimator
- Doubly Robust (DR) estimator

#### Inference

- Bootstrap standard errors with parallel support
- Influence function-based analytic standard errors
- Cross-fitting (sample splitting) for DR estimation
- Confidence intervals via normal approximation or percentile bootstrap

#### Model Selection

- [`cf_cv()`](https://boyercb.github.io/cfperformance/reference/cf_cv.md) -
  K-fold cross-validation with counterfactual metrics
- [`cf_compare()`](https://boyercb.github.io/cfperformance/reference/cf_compare.md) -
  Compare multiple prediction models

#### Data

- `cvd_sim` - Simulated cardiovascular disease dataset for examples

### References

Boyer CB, Dahabreh IJ, Steingrimsson JA. Counterfactual prediction model
performance. *Statistics in Medicine*. 2025; 44(23-24):e70287.
<doi:10.1002/sim.70287>
