# Estimate (Counterfactual) False Positive Rate in the Target Population

Estimates the false positive rate (1 - specificity) of a binary
classifier. This is a convenience wrapper around
[`tr_specificity()`](https://boyercb.github.io/cfperformance/reference/tr_specificity.md).

## Usage

``` r
tr_fpr(
  predictions,
  outcomes,
  treatment,
  source,
  covariates,
  threshold = 0.5,
  treatment_level = 0,
  analysis = c("transport", "joint"),
  estimator = c("dr", "om", "ipw", "naive"),
  selection_model = NULL,
  propensity_model = NULL,
  outcome_model = NULL,
  se_method = c("none", "bootstrap"),
  n_boot = 200,
  conf_level = 0.95,
  stratified_boot = TRUE,
  ps_trim = NULL,
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

- source:

  Numeric vector of population indicators (1=source/RCT, 0=target).

- covariates:

  A matrix or data frame of baseline covariates.

- threshold:

  Numeric vector of classification thresholds. Predictions above this
  value are classified as positive. Can be a single value or a vector
  for computing sensitivity at multiple thresholds simultaneously.
  Default is 0.5.

- treatment_level:

  The treatment level of interest (default: 0).

- analysis:

  Character string specifying the type of analysis:

  - `"transport"`: Use source outcomes for target estimation (default)

  - `"joint"`: Pool source and target data

- estimator:

  Character string specifying the estimator:

  - `"naive"`: Naive estimator (biased)

  - `"om"`: Outcome model estimator

  - `"ipw"`: Inverse probability weighting estimator

  - `"dr"`: Doubly robust estimator (default)

- selection_model:

  Optional fitted selection model for P(S=0\|X). If NULL, a logistic
  regression model is fit using the covariates.

- propensity_model:

  Optional fitted propensity score model for P(A=1\|X,S=1). If NULL, a
  logistic regression model is fit using source data.

- outcome_model:

  Optional fitted outcome model for E\[L(Y,g)\|X,A,S\]. If NULL, a
  regression model is fit using the relevant data. For binary outcomes,
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

- stratified_boot:

  Logical indicating whether to use stratified bootstrap that preserves
  the source/target ratio (default: TRUE). Recommended for
  transportability analysis.

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

- parallel:

  Logical indicating whether to use parallel processing for bootstrap
  (default: FALSE).

- ncores:

  Number of cores for parallel processing (default: NULL, which uses all
  available cores minus one).

- ...:

  Additional arguments passed to internal functions.

## Value

An object of class `c("tr_fpr", "tr_performance")` with the same
structure as
[`tr_specificity()`](https://boyercb.github.io/cfperformance/reference/tr_specificity.md),
but with `estimate` containing 1 - specificity.

## See also

[`tr_specificity()`](https://boyercb.github.io/cfperformance/reference/tr_specificity.md),
[`tr_sensitivity()`](https://boyercb.github.io/cfperformance/reference/tr_sensitivity.md)

## Examples

``` r
set.seed(123)
n <- 500
x <- rnorm(n)
s <- rbinom(n, 1, 0.6)
a <- rbinom(n, 1, 0.5)
y <- rbinom(n, 1, plogis(-1 + x))
pred <- plogis(-1 + 0.8 * x)

result <- tr_fpr(
  predictions = pred, outcomes = y, treatment = a,
  source = s, covariates = data.frame(x = x),
  se_method = "none"
)
print(result)
#> 
#> Transportable False Positive Rate Estimate
#> ==========================================
#> 
#> Estimator: DR 
#> Analysis: transport 
#> Treatment level: 0 
#> N (source): 312 
#> N (target): 188 
#> 
#> Threshold: 0.5 
#> Estimate: 0.0782 
#> Naive estimate: 0.0833 
```
