# Estimate Counterfactual False Positive Rate

Estimates the false positive rate (FPR) of a binary classifier, which is
the complement of specificity (FPR = 1 - specificity).

## Usage

``` r
cf_fpr(
  predictions,
  outcomes,
  treatment,
  covariates,
  threshold = 0.5,
  treatment_level = 0,
  estimator = c("dr", "cl", "ipw", "naive"),
  propensity_model = NULL,
  outcome_model = NULL,
  se_method = c("none", "bootstrap"),
  n_boot = 200,
  conf_level = 0.95,
  cross_fit = FALSE,
  n_folds = 5,
  parallel = FALSE,
  ncores = NULL,
  ps_trim = NULL,
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

- threshold:

  Numeric vector of classification thresholds. Predictions above this
  value are classified as positive. Can be a single value or a vector
  for computing sensitivity at multiple thresholds simultaneously.
  Default is 0.5.

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

- se_method:

  Method for standard error estimation:

  - `"bootstrap"`: Bootstrap standard errors (default)

  - `"influence"`: Influence function-based standard errors

  - `"none"`: No standard error estimation

- n_boot:

  Number of bootstrap replications (default: 500).

- conf_level:

  Confidence level for intervals (default: 0.95).

- cross_fit:

  Logical indicating whether to use cross-fitting for nuisance model
  estimation (default: FALSE). Cross-fitting enables valid inference
  when using flexible machine learning estimators.

- n_folds:

  Number of folds for cross-fitting (default: 5).

- parallel:

  Logical indicating whether to use parallel processing for bootstrap
  (default: FALSE).

- ncores:

  Number of cores for parallel processing (default: NULL, which uses all
  available cores minus one).

- ps_trim:

  Propensity score trimming specification. Controls how extreme
  propensity scores are handled. Can be:

  - `NULL` (default): Uses absolute bounds `c(0.01, 0.99)`

  - `"none"`: No trimming applied

  - `"quantile"`: Quantile-based trimming with default `c(0.01, 0.99)`

  - `"absolute"`: Explicit absolute bounds with default `c(0.01, 0.99)`

  - A numeric vector of length 2: `c(lower, upper)` absolute bounds

  - A single numeric: Symmetric bounds `c(x, 1-x)`

  - A list with `method` ("absolute"/"quantile"/"none") and `bounds`

- ...:

  Additional arguments passed to internal functions.

## Value

An object of class `c("cf_fpr", "cf_performance")` containing:

- estimate:

  Point estimate(s) of counterfactual FPR

- se:

  Standard error(s) (if computed)

- ci_lower:

  Lower confidence interval bound(s)

- ci_upper:

  Upper confidence interval bound(s)

- threshold:

  Threshold value(s) used

- estimator:

  Estimator used

- naive_estimate:

  Naive FPR for comparison

## Details

False positive rate is defined as: \$\$FPR(c) = P(\hat{Y} \> c \|
Y^{(a)} = 0) = 1 - Specificity(c)\$\$

This function is provided as a convenience for ROC curve construction,
where FPR is typically plotted on the x-axis.

## See also

[`cf_sensitivity()`](https://boyercb.github.io/cfperformance/reference/cf_sensitivity.md),
[`cf_specificity()`](https://boyercb.github.io/cfperformance/reference/cf_specificity.md),
[`cf_tpr()`](https://boyercb.github.io/cfperformance/reference/cf_sensitivity.md)

## Examples

``` r
set.seed(123)
n <- 500
x <- rnorm(n)
a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
pred <- plogis(-1 + 0.8 * x)

result <- cf_fpr(
  predictions = pred,
  outcomes = y,
  treatment = a,
  covariates = data.frame(x = x),
  se_method = "none"
)
print(result)
#> 
#> Counterfactual False Positive Rate Estimate
#> ============================================
#> 
#> Estimator: DR 
#> Treatment level: 0 
#> N: 500 
#> 
#> Threshold: 0.5 
#> Estimate: 0.0206 
#> Naive estimate: 0.0466 
#> 
```
