# Estimate Counterfactual Specificity

Estimates the specificity (true negative rate) of a binary classifier at
one or more thresholds under a hypothetical intervention where treatment
is set to a specific level.

## Usage

``` r
cf_specificity(
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
  ...
)

cf_tnr(
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

  Optional fitted outcome model. If NULL, a logistic regression model is
  fit using the covariates among treated/untreated.

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

- ...:

  Additional arguments passed to internal functions.

## Value

An object of class `c("cf_specificity", "cf_performance")` containing:

- estimate:

  Point estimate(s) of counterfactual specificity

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

  Naive specificity for comparison

- n_obs:

  Number of observations

- treatment_level:

  Counterfactual treatment level

## Details

Specificity (also known as true negative rate) is defined as:
\$\$Specificity(c) = P(\hat{Y} \leq c \| Y^{(a)} = 0)\$\$

where \\Y^{(a)}\\ is the potential outcome under treatment level \\a\\.
The estimators mirror those for sensitivity (see
[`cf_sensitivity()`](https://boyercb.github.io/cfperformance/reference/cf_sensitivity.md)).

## References

Coston, A., Mishler, A., Kennedy, E. H., & Chouldechova, A. (2020).
"Counterfactual risk assessments, evaluation, and fairness."
*Proceedings of the 2020 Conference on Fairness, Accountability, and
Transparency*, 582-593.

## See also

[`cf_sensitivity()`](https://boyercb.github.io/cfperformance/reference/cf_sensitivity.md),
[`cf_fpr()`](https://boyercb.github.io/cfperformance/reference/cf_fpr.md),
[`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md)

## Examples

``` r
# Generate example data
set.seed(123)
n <- 1000
x <- rnorm(n)
a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
pred <- plogis(-1 + 0.8 * x)

# Estimate counterfactual specificity at default threshold (0.5)
result <- cf_specificity(
  predictions = pred,
  outcomes = y,
  treatment = a,
  covariates = data.frame(x = x),
  treatment_level = 0,
  estimator = "dr",
  se_method = "none"
)
print(result)
#> 
#> Counterfactual Specificity Estimate
#> ====================================
#> 
#> Estimator: DR 
#> Treatment level: 0 
#> N: 1000 
#> 
#> Threshold: 0.5 
#> Estimate: 0.9543 
#> Naive estimate: 0.9476 
#> 
```
