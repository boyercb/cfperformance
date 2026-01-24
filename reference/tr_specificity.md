# Estimate (Counterfactual) Specificity in the Target Population

Estimates the specificity (true negative rate) of a binary classifier at
one or more thresholds in a target population using data transported
from a source population (typically an RCT).

## Usage

``` r
tr_specificity(
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
  parallel = FALSE,
  ncores = NULL,
  ...
)

tr_tnr(
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
  regression model is fit using the relevant data.

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

- parallel:

  Logical indicating whether to use parallel processing for bootstrap
  (default: FALSE).

- ncores:

  Number of cores for parallel processing (default: NULL, which uses all
  available cores minus one).

- ...:

  Additional arguments passed to internal functions.

## Value

An object of class `c("tr_specificity", "tr_performance")` containing:

- estimate:

  Point estimate(s) of transportable specificity

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

- analysis:

  Analysis type

- naive_estimate:

  Naive specificity for comparison

- n_target:

  Number of target observations

- n_source:

  Number of source observations

- treatment_level:

  Treatment level

## Details

Specificity (also known as true negative rate) is defined as:
\$\$Specificity(c) = P(\hat{Y} \le c \| Y = 0)\$\$

In the transportability setting, we estimate specificity in the target
population using outcome data from the source population. The estimators
mirror those for sensitivity (see
[`tr_sensitivity()`](https://boyercb.github.io/cfperformance/reference/tr_sensitivity.md)).

## References

Steingrimsson, J. A., Wen, L., Voter, S., & Dahabreh, I. J. (2024).
"Interpretable meta-analysis of model or marker performance." *arXiv
preprint arXiv:2409.13458*.

## See also

[`tr_sensitivity()`](https://boyercb.github.io/cfperformance/reference/tr_sensitivity.md),
[`tr_fpr()`](https://boyercb.github.io/cfperformance/reference/tr_fpr.md),
[`tr_auc()`](https://boyercb.github.io/cfperformance/reference/tr_auc.md)

## Examples

``` r
# Generate example data
set.seed(123)
n <- 1000
x <- rnorm(n)
s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
pred <- plogis(-1 + 0.8 * x)

# Estimate transportable specificity
result <- tr_specificity(
  predictions = pred,
  outcomes = y,
  treatment = a,
  source = s,
  covariates = data.frame(x = x),
  treatment_level = 0,
  estimator = "dr",
  se_method = "none"
)
print(result)
#> 
#> Transportable Specificity Estimate
#> ===================================
#> 
#> Estimator: DR 
#> Analysis: transport 
#> Treatment level: 0 
#> N (source): 624 
#> N (target): 376 
#> 
#> Threshold: 0.5 
#> Estimate: 0.9196 
#> Naive estimate: 0.9397 
```
