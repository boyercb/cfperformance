# Compute Counterfactual ROC Curve

Computes a receiver operating characteristic (ROC) curve under a
hypothetical intervention where treatment is set to a specific level.

## Usage

``` r
cf_roc(
  predictions,
  outcomes,
  treatment,
  covariates,
  treatment_level = 0,
  estimator = c("dr", "om", "ipw", "naive"),
  propensity_model = NULL,
  outcome_model = NULL,
  n_thresholds = 201,
  thresholds = NULL,
  include_naive = TRUE,
  ...
)
```

## Arguments

- predictions:

  Numeric vector of model predictions.

- outcomes:

  Numeric vector of observed outcomes.

- treatment:

  Numeric vector of treatment indicators (0/1).

- covariates:

  A matrix or data frame of baseline covariates (confounders).

- treatment_level:

  The counterfactual treatment level (default: 0).

- estimator:

  Character string specifying the estimator:

  - `"naive"`: Naive estimator (biased)

  - `"cl"`: Conditional loss estimator

  - `"ipw"`: Inverse probability weighting estimator

  - `"dr"`: Doubly robust estimator (default)

- propensity_model:

  Optional fitted propensity score model. If NULL, a logistic regression
  model is fit using the covariates.

- outcome_model:

  Optional fitted outcome model. If NULL, a regression model is fit
  using the covariates among treated/untreated. For binary outcomes,
  this should be a model for E\[Y\|X,A\] (binomial family). For
  continuous outcomes, this should be a model for E\[L\|X,A\] (gaussian
  family).

- n_thresholds:

  Integer specifying the number of thresholds to evaluate. Thresholds
  are evenly spaced between 0 and 1. Default is 201.

- thresholds:

  Optional numeric vector of specific thresholds to use. If provided,
  overrides `n_thresholds`.

- include_naive:

  Logical indicating whether to also compute the naive ROC curve for
  comparison. Default is TRUE.

- ...:

  Additional arguments passed to internal functions.

## Value

An object of class `c("cf_roc", "roc_curve")` containing:

- thresholds:

  Thresholds used

- sensitivity:

  Sensitivity (TPR) at each threshold

- fpr:

  False positive rate at each threshold

- specificity:

  Specificity at each threshold

- naive_sensitivity:

  Naive sensitivity (if include_naive=TRUE)

- naive_fpr:

  Naive FPR (if include_naive=TRUE)

- auc:

  Area under the ROC curve (computed via trapezoidal rule)

- naive_auc:

  Naive AUC (if include_naive=TRUE)

- estimator:

  Estimator used

- n_obs:

  Number of observations

## Details

The ROC curve plots sensitivity (true positive rate) against the false
positive rate (1 - specificity) at various classification thresholds.

This function computes counterfactual sensitivity and FPR at multiple
thresholds using the estimators from
[`cf_sensitivity()`](https://boyercb.github.io/cfperformance/reference/cf_sensitivity.md)
and
[`cf_fpr()`](https://boyercb.github.io/cfperformance/reference/cf_fpr.md).
The area under the curve (AUC) is computed using the trapezoidal rule on
the discrete threshold grid. For exact AUC estimation, use
[`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md)
which employs the Wilcoxon-Mann-Whitney statistic.

## References

Coston, A., Mishler, A., Kennedy, E. H., & Chouldechova, A. (2020).
"Counterfactual risk assessments, evaluation, and fairness."
*Proceedings of the 2020 Conference on Fairness, Accountability, and
Transparency*, 582-593.

## See also

[`cf_sensitivity()`](https://boyercb.github.io/cfperformance/reference/cf_sensitivity.md),
[`cf_specificity()`](https://boyercb.github.io/cfperformance/reference/cf_specificity.md),
[`cf_fpr()`](https://boyercb.github.io/cfperformance/reference/cf_fpr.md),
[`plot.cf_roc()`](https://boyercb.github.io/cfperformance/reference/plot.tr_roc.md),
[`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md)

## Examples

``` r
# Generate example data
set.seed(123)
n <- 500
x <- rnorm(n)
a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
pred <- plogis(-1 + 0.8 * x)

# Compute counterfactual ROC curve
roc <- cf_roc(
  predictions = pred,
  outcomes = y,
  treatment = a,
  covariates = data.frame(x = x),
  n_thresholds = 51
)
print(roc)
#> 
#> Counterfactual ROC Curve
#> ========================
#> 
#> Estimator: DR 
#> Treatment level: 0 
#> N: 500 
#> Thresholds evaluated: 51 
#> 
#> AUC: 0.7326 
#> Naive AUC: 0.7283 
#> 
#> Use plot() to visualize the ROC curve.
#> 

# Plot the ROC curve
plot(roc)
```
