# cfperformance 0.3.0

Adds machine learning integration for flexible nuisance model estimation with
automatic cross-fitting for valid inference.

## New Features

### ML Learner Interface
* `ml_learner()` - Specify ML methods for propensity score and outcome models
* Supports: `ranger`, `xgboost`, `grf`, `glmnet`, `superlearner`, and `custom`
* Automatic cross-fitting when `ml_learner` specs are detected
* Seamlessly integrates with existing `propensity_model`/`outcome_model` arguments

### Cross-Fitting Support
* `cf_mse()` - Full support for ML learners with cross-fitting
* `cf_auc()` - Full support for ML learners with cross-fitting (DR estimator)

### Supported Learners
* **ranger** - Fast random forest implementation
* **xgboost** - Gradient boosting (XGBoost)
* **grf** - Generalized random forests with honest estimation
* **glmnet** - Elastic net regularization with CV-selected λ
* **superlearner** - Ensemble learning
* **custom** - User-supplied fit/predict functions

### Usage Example
```r
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

### Documentation
* New vignette: "Machine Learning Integration"
* Updated README with ML integration examples

## References

Chernozhukov, V., et al. (2018). Double/debiased machine learning for treatment
and structural parameters. *The Econometrics Journal*, 21(1), C1-C68.

Li, B., Gatsonis, C., Dahabreh, I. J., & Steingrimsson, J. A. (2022).
Estimating the area under the ROC curve when transporting a prediction
model to a target population. *Biometrics*, 79(3), 2343-2356.

---

# cfperformance 0.2.0

Major release adding transportability estimators from Voter et al. (2025) for 
evaluating prediction model performance when transporting from a source 
population (e.g., RCT) to a target population.

## New Features

### Transportability Functions
* `tr_mse()` - Transportable MSE estimation with naive, om, ipw, dr estimators
* `tr_auc()` - Transportable AUC estimation with all estimators
* `tr_calibration()` - Transportable calibration curves with ICI, E50, E90, Emax

### Analysis Modes
* **Transport analysis** (`analysis = "transport"`): Use source/RCT outcomes to 
  estimate performance in target population
* **Joint analysis** (`analysis = "joint"`): Pool source and target data for 
  potentially more efficient estimation

### Inference
* Bootstrap standard errors with stratified sampling option (`stratified_boot = TRUE`)
  to preserve source/target ratio
* Influence function-based standard errors for `tr_mse()` (all estimators)

### New Data
* `transport_sim` - Simulated dataset with source (RCT) and target populations
  for transportability examples (n = 1,500)

### Documentation
* New vignette: "Transportability Analysis with cfperformance"
* Complete roxygen documentation for all new functions
* Updated README with transportability examples

## S3 Methods for tr_* Functions
* `print.tr_performance()` - Print method for transportability results
* `summary.tr_performance()` - Detailed summary
* `coef.tr_performance()` - Extract point estimates
* `confint.tr_performance()` - Confidence intervals
* `plot.tr_calibration()` - Calibration curve visualization

## References

Voter SR, et al. Transportability of machine learning-based counterfactual 
prediction models with application to CASS. *Diagnostic and Prognostic Research*. 
2025; 9(4). doi:10.1186/s41512-025-00201-y

---

# cfperformance 0.1.0

Initial release implementing methods from Boyer, Dahabreh & Steingrimsson (2025),
"Estimating and evaluating counterfactual prediction models."

## Features

### Performance Metrics
* `cf_mse()` - Counterfactual MSE/Brier score estimation
* `cf_auc()` - Counterfactual AUC estimation  
* `cf_calibration()` - Counterfactual calibration curves with ICI, E50, E90, Emax

### Estimators
* Naive estimator (subset-based)
* Conditional Loss (CL) / Outcome modeling estimator
* Inverse Probability Weighting (IPW) estimator
* Doubly Robust (DR) estimator

### Inference
* Bootstrap standard errors with parallel support
* Influence function-based analytic standard errors
* Cross-fitting (sample splitting) for DR estimation
* Confidence intervals via normal approximation or percentile bootstrap

### Model Selection
* `cf_cv()` - K-fold cross-validation with counterfactual metrics
* `cf_compare()` - Compare multiple prediction models

### Data
* `cvd_sim` - Simulated cardiovascular disease dataset for examples

## References

Boyer CB, Dahabreh IJ, Steingrimsson JA. Counterfactual prediction model 
performance. *Statistics in Medicine*. 2025; 44(23-24):e70287. 
doi:10.1002/sim.70287
