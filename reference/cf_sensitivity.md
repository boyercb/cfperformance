# Estimate Counterfactual Sensitivity

Estimates the sensitivity (true positive rate) of a binary classifier at
one or more thresholds under a hypothetical intervention where treatment
is set to a specific level.

## Usage

``` r
cf_sensitivity(
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

cf_tpr(
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

An object of class `c("cf_sensitivity", "cf_performance")` containing:

- estimate:

  Point estimate(s) of counterfactual sensitivity

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

  Naive sensitivity for comparison

- n_obs:

  Number of observations

- treatment_level:

  Counterfactual treatment level

## Details

Sensitivity (also known as true positive rate or recall) is defined as:
\$\$Sensitivity(c) = P(\hat{Y} \> c \| Y^{(a)} = 1)\$\$

where \\Y^{(a)}\\ is the potential outcome under treatment level \\a\\.

The function implements three estimators following Coston et al. (2020):

**Conditional Loss (CL) Estimator**: Weights by the predicted
probability of being a case under the counterfactual:
\$\$\hat{\psi}\_{sens,cl} = \frac{\sum_i I(\hat{h}(X_i) \> c)
\hat{m}(X_i)}{\sum_i \hat{m}(X_i)}\$\$ where \\\hat{m}(X) = P(Y=1\|X,
A=a)\\.

**IPW Estimator**: Weights by the inverse probability of treatment:
\$\$\hat{\psi}\_{sens,ipw} = \frac{\sum_i I(\hat{h}(X_i) \> c, Y_i=1,
A_i=a) / \hat{e}(X_i)}{\sum_i I(Y_i=1, A_i=a) / \hat{e}(X_i)}\$\$ where
\\\hat{e}(X) = P(A=a\|X)\\.

**Doubly Robust (DR) Estimator**: Combines CL and IPW for protection
against misspecification of either model.

## References

Coston, A., Mishler, A., Kennedy, E. H., & Chouldechova, A. (2020).
"Counterfactual risk assessments, evaluation, and fairness."
*Proceedings of the 2020 Conference on Fairness, Accountability, and
Transparency*, 582-593.

Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
"Estimating and evaluating counterfactual prediction models."
*Statistics in Medicine*, 44(23-24), e70287.
[doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

## See also

[`cf_specificity()`](https://boyercb.github.io/cfperformance/reference/cf_specificity.md),
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

# Estimate counterfactual sensitivity at default threshold (0.5)
result <- cf_sensitivity(
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
#> Counterfactual Sensitivity Estimate
#> ====================================
#> 
#> Estimator: DR 
#> Treatment level: 0 
#> N: 1000 
#> 
#> Threshold: 0.5 
#> Estimate: 0.21 
#> Naive estimate: 0.2291 
#> 

# Estimate at multiple thresholds
result_multi <- cf_sensitivity(
  predictions = pred,
  outcomes = y,
  treatment = a,
  covariates = data.frame(x = x),
  threshold = c(0.3, 0.5, 0.7),
  se_method = "none"
)
print(result_multi)
#> 
#> Counterfactual Sensitivity Estimate
#> ====================================
#> 
#> Estimator: DR 
#> Treatment level: 0 
#> N: 1000 
#> 
#> Results by threshold:
#>  Threshold Estimate  Naive
#>        0.3   0.6649 0.6545
#>        0.5   0.2100 0.2291
#>        0.7   0.0416 0.0400
#> 
```
