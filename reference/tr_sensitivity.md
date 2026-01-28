# Estimate Transportable Sensitivity in the Target Population

Estimates the sensitivity (true positive rate) of a binary classifier at
one or more thresholds in a target population using data transported
from a source population. Supports both **counterfactual** (under
hypothetical intervention) and **factual** (observational) prediction
model transportability.

## Usage

``` r
tr_sensitivity(
  predictions,
  outcomes,
  treatment = NULL,
  source,
  covariates,
  threshold = 0.5,
  treatment_level = NULL,
  analysis = c("transport", "joint"),
  estimator = c("dr", "om", "ipw", "naive"),
  selection_model = NULL,
  propensity_model = NULL,
  outcome_model = NULL,
  se_method = c("none", "bootstrap", "influence"),
  n_boot = 200,
  conf_level = 0.95,
  stratified_boot = TRUE,
  cross_fit = FALSE,
  n_folds = 5,
  ps_trim = NULL,
  parallel = FALSE,
  ncores = NULL,
  ...
)

tr_tpr(
  predictions,
  outcomes,
  treatment = NULL,
  source,
  covariates,
  threshold = 0.5,
  treatment_level = NULL,
  analysis = c("transport", "joint"),
  estimator = c("dr", "om", "ipw", "naive"),
  selection_model = NULL,
  propensity_model = NULL,
  outcome_model = NULL,
  se_method = c("none", "bootstrap", "influence"),
  n_boot = 200,
  conf_level = 0.95,
  stratified_boot = TRUE,
  cross_fit = FALSE,
  n_folds = 5,
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

  Numeric vector of treatment indicators (0/1), or `NULL` for factual
  prediction model transportability (no treatment/intervention). When
  `NULL`, only the selection model is used for weighting.

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

  The treatment level of interest (default: `NULL`). Required when
  `treatment` is provided; should be `NULL` when `treatment` is `NULL`
  (factual mode).

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

- cross_fit:

  Logical indicating whether to use cross-fitting for nuisance model
  estimation (default: FALSE).

- n_folds:

  Number of folds for cross-fitting (default: 5).

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

An object of class `c("tr_sensitivity", "tr_performance")` containing:

- estimate:

  Point estimate(s) of transportable sensitivity

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

  Naive sensitivity for comparison

- n_target:

  Number of target observations

- n_source:

  Number of source observations

- treatment_level:

  Treatment level (NULL for factual mode)

## Details

Sensitivity (also known as true positive rate or recall) is defined as:
\$\$Sensitivity(c) = P(\hat{Y} \> c \| Y = 1)\$\$

### Counterfactual Mode (treatment provided)

When `treatment` is specified, estimates sensitivity for counterfactual
outcomes under a hypothetical intervention. Requires selection,
propensity, and outcome models.

### Factual Mode (treatment = NULL)

When `treatment` is `NULL`, estimates sensitivity for observed outcomes
in the target population using only the selection model for inverse-odds
weighting. This is appropriate for factual prediction model
transportability.

### Estimators

**Outcome Model (OM) Estimator**: \$\$\hat{\psi}\_{sens,om} =
\frac{\sum_i I(S_i=0) I(\hat{h}(X_i) \> c) \hat{m}(X_i)}{\sum_i I(S_i=0)
\hat{m}(X_i)}\$\$ where \\\hat{m}(X) \approx P(Y=1\|X, R=1)\\.

**IPW Estimator**: \$\$\hat{\psi}\_{sens,ipw} = \frac{\sum_i
I(\hat{h}(X_i) \> c, Y_i=1, R_i=1) \hat{w}(X_i)}{\sum_i I(Y_i=1, R_i=1)
\hat{w}(X_i)}\$\$ where \\\hat{w}(X) \approx P(R=0\|X)/P(R=1\|X)\\.

**Doubly Robust (DR) Estimator**: Combines OM and IPW for protection
against misspecification of either model.

## References

Steingrimsson, J. A., et al. (2023). "Transporting a Prediction Model
for Use in a New Target Population." *American Journal of Epidemiology*,
192(2), 296-304.
[doi:10.1093/aje/kwac128](https://doi.org/10.1093/aje/kwac128)

Steingrimsson, J. A., Wen, L., Voter, S., & Dahabreh, I. J. (2024).
"Interpretable meta-analysis of model or marker performance." *arXiv
preprint arXiv:2409.13458*.

Coston, A., Mishler, A., Kennedy, E. H., & Chouldechova, A. (2020).
"Counterfactual risk assessments, evaluation, and fairness."
*Proceedings of the 2020 Conference on Fairness, Accountability, and
Transparency*, 582-593.

## See also

[`tr_specificity()`](https://boyercb.github.io/cfperformance/reference/tr_specificity.md),
[`tr_fpr()`](https://boyercb.github.io/cfperformance/reference/tr_fpr.md),
[`tr_auc()`](https://boyercb.github.io/cfperformance/reference/tr_auc.md)

## Examples

``` r
# Generate example data with source (RCT) and target populations
set.seed(123)
n <- 1000
x <- rnorm(n)
s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
pred <- plogis(-1 + 0.8 * x)

# Estimate transportable sensitivity at default threshold (0.5)
result <- tr_sensitivity(
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
#> Transportable Sensitivity Estimate
#> ===================================
#> 
#> Estimator: DR 
#> Analysis: transport 
#> Treatment level: 0 
#> N (source): 624 
#> N (target): 376 
#> 
#> Threshold: 0.5 
#> Estimate: 0.2177 
#> Naive estimate: 0.2254 

# Estimate at multiple thresholds
result_multi <- tr_sensitivity(
  predictions = pred,
  outcomes = y,
  treatment = a,
  source = s,
  covariates = data.frame(x = x),
  threshold = c(0.3, 0.5, 0.7),
  se_method = "none"
)
print(result_multi)
#> 
#> Transportable Sensitivity Estimate
#> ===================================
#> 
#> Estimator: DR 
#> Analysis: transport 
#> Treatment level: 1 
#> N (source): 624 
#> N (target): 376 
#> 
#> Results by threshold:
#>  Threshold Estimate  Naive
#>        0.3   0.8167 0.7097
#>        0.5   0.3252 0.2903
#>        0.7   0.0756 0.0484
```
