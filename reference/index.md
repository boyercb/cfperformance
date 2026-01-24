# Package index

## Counterfactual Performance (Single Population)

Estimate prediction model performance under counterfactual interventions

- [`cf_mse()`](https://boyercb.github.io/cfperformance/reference/cf_mse.md)
  : Estimate Counterfactual Mean Squared Error
- [`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md)
  : Estimate Counterfactual Area Under the ROC Curve
- [`cf_calibration()`](https://boyercb.github.io/cfperformance/reference/cf_calibration.md)
  : Estimate Counterfactual Calibration Curve
- [`cf_sensitivity()`](https://boyercb.github.io/cfperformance/reference/cf_sensitivity.md)
  [`cf_tpr()`](https://boyercb.github.io/cfperformance/reference/cf_sensitivity.md)
  : Estimate Counterfactual Sensitivity
- [`cf_specificity()`](https://boyercb.github.io/cfperformance/reference/cf_specificity.md)
  [`cf_tnr()`](https://boyercb.github.io/cfperformance/reference/cf_specificity.md)
  : Estimate Counterfactual Specificity
- [`cf_fpr()`](https://boyercb.github.io/cfperformance/reference/cf_fpr.md)
  : Estimate Counterfactual False Positive Rate
- [`cf_roc()`](https://boyercb.github.io/cfperformance/reference/cf_roc.md)
  : Compute Counterfactual ROC Curve

## Transportability (Two Populations)

Transport prediction model performance from source to target population

- [`tr_mse()`](https://boyercb.github.io/cfperformance/reference/tr_mse.md)
  : Estimate (Counterfactual) Mean Squared Error in the Target
  Population
- [`tr_auc()`](https://boyercb.github.io/cfperformance/reference/tr_auc.md)
  : Estimate (Counterfactual) Area Under the ROC Curve in the Target
  Population
- [`tr_calibration()`](https://boyercb.github.io/cfperformance/reference/tr_calibration.md)
  : Estimate (Counterfactual) Calibration in the Target Population
- [`tr_sensitivity()`](https://boyercb.github.io/cfperformance/reference/tr_sensitivity.md)
  [`tr_tpr()`](https://boyercb.github.io/cfperformance/reference/tr_sensitivity.md)
  : Estimate (Counterfactual) Sensitivity in the Target Population
- [`tr_specificity()`](https://boyercb.github.io/cfperformance/reference/tr_specificity.md)
  [`tr_tnr()`](https://boyercb.github.io/cfperformance/reference/tr_specificity.md)
  : Estimate (Counterfactual) Specificity in the Target Population
- [`tr_fpr()`](https://boyercb.github.io/cfperformance/reference/tr_fpr.md)
  : Estimate (Counterfactual) False Positive Rate in the Target
  Population
- [`tr_roc()`](https://boyercb.github.io/cfperformance/reference/tr_roc.md)
  : Compute Transportable ROC Curve

## Model Selection

Cross-validation and model comparison

- [`cf_cv()`](https://boyercb.github.io/cfperformance/reference/cf_cv.md)
  : Cross-Validation with Counterfactual Performance Metrics
- [`cf_compare()`](https://boyercb.github.io/cfperformance/reference/cf_compare.md)
  : Compare Multiple Prediction Models

## Machine Learning

ML learner specifications for nuisance models

- [`ml_learner()`](https://boyercb.github.io/cfperformance/reference/ml_learner.md)
  : Specify a Machine Learning Learner for Nuisance Models
- [`print(`*`<ml_learner>`*`)`](https://boyercb.github.io/cfperformance/reference/print.ml_learner.md)
  : Print Method for ml_learner Objects

## Nuisance Models

Fit propensity and outcome models

- [`fit_nuisance()`](https://boyercb.github.io/cfperformance/reference/fit_nuisance.md)
  : Fit Nuisance Models for Counterfactual Performance Estimation

## S3 Methods

Methods for result objects

- [`as.data.frame(`*`<tr_roc>`*`)`](https://boyercb.github.io/cfperformance/reference/as.data.frame.tr_roc.md)
  [`as.data.frame(`*`<cf_roc>`*`)`](https://boyercb.github.io/cfperformance/reference/as.data.frame.tr_roc.md)
  : Convert ROC Curve to Data Frame
- [`coef(`*`<cf_performance>`*`)`](https://boyercb.github.io/cfperformance/reference/coef.cf_performance.md)
  : Coefficient Method for cf_performance Objects
- [`coef(`*`<tr_performance>`*`)`](https://boyercb.github.io/cfperformance/reference/coef.tr_performance.md)
  : Coefficient Method for tr_performance Objects
- [`confint(`*`<cf_performance>`*`)`](https://boyercb.github.io/cfperformance/reference/confint.cf_performance.md)
  : Confidence Interval Method for cf_performance Objects
- [`confint(`*`<tr_performance>`*`)`](https://boyercb.github.io/cfperformance/reference/confint.tr_performance.md)
  : Confidence Interval Method for tr_performance Objects
- [`plot(`*`<cf_auc>`*`)`](https://boyercb.github.io/cfperformance/reference/plot.cf_auc.md)
  : Plot Method for cf_auc Objects
- [`plot(`*`<cf_calibration>`*`)`](https://boyercb.github.io/cfperformance/reference/plot.cf_calibration.md)
  : Plot Method for cf_calibration Objects
- [`plot(`*`<tr_calibration>`*`)`](https://boyercb.github.io/cfperformance/reference/plot.tr_calibration.md)
  : Plot Method for tr_calibration Objects
- [`plot(`*`<tr_roc>`*`)`](https://boyercb.github.io/cfperformance/reference/plot.tr_roc.md)
  [`plot(`*`<cf_roc>`*`)`](https://boyercb.github.io/cfperformance/reference/plot.tr_roc.md)
  : Plot ROC Curve
- [`print(`*`<cf_compare>`*`)`](https://boyercb.github.io/cfperformance/reference/print.cf_compare.md)
  : Print Method for cf_compare Objects
- [`print(`*`<cf_cv>`*`)`](https://boyercb.github.io/cfperformance/reference/print.cf_cv.md)
  : Print Method for cf_cv Objects
- [`print(`*`<cf_nuisance>`*`)`](https://boyercb.github.io/cfperformance/reference/print.cf_nuisance.md)
  : Print Method for cf_nuisance Objects
- [`print(`*`<cf_performance>`*`)`](https://boyercb.github.io/cfperformance/reference/print.cf_performance.md)
  : Print Method for cf_performance Objects
- [`print(`*`<tr_performance>`*`)`](https://boyercb.github.io/cfperformance/reference/print.tr_performance.md)
  : Print Method for tr_performance Objects
- [`summary(`*`<cf_compare>`*`)`](https://boyercb.github.io/cfperformance/reference/summary.cf_compare.md)
  : Summary Method for cf_compare Objects
- [`summary(`*`<cf_cv>`*`)`](https://boyercb.github.io/cfperformance/reference/summary.cf_cv.md)
  : Summary Method for cf_cv Objects
- [`summary(`*`<cf_performance>`*`)`](https://boyercb.github.io/cfperformance/reference/summary.cf_performance.md)
  : Summary Method for cf_performance Objects
- [`summary(`*`<tr_performance>`*`)`](https://boyercb.github.io/cfperformance/reference/summary.tr_performance.md)
  : Summary Method for tr_performance Objects

## Example Data

- [`cvd_sim`](https://boyercb.github.io/cfperformance/reference/cvd_sim.md)
  : Simulated Cardiovascular Disease Data
- [`transport_sim`](https://boyercb.github.io/cfperformance/reference/transport_sim.md)
  : Simulated Transportability Data
