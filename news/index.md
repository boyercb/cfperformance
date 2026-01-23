# Changelog

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

#### Supported Learners

- **ranger** - Fast random forest implementation
- **xgboost** - Gradient boosting (XGBoost)
- **grf** - Generalized random forests with honest estimation
- **glmnet** - Elastic net regularization with CV-selected λ
- **superlearner** - Ensemble learning
- **custom** - User-supplied fit/predict functions

#### Usage Example

``` r
cf_mse(
  predictions = pred, outcomes = y, treatment = a, covariates = df,
  propensity_model = ml_learner("ranger", num.trees = 500),
  outcome_model = ml_learner("xgboost", nrounds = 100),
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
