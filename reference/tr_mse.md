# Estimate Transportable Mean Squared Error in the Target Population

Estimates the mean squared error (MSE) of a prediction model in a target
population using data transported from a source population. Supports
both **counterfactual** (under hypothetical intervention) and
**factual** (observational) prediction model transportability.

## Usage

``` r
tr_mse(
  predictions,
  outcomes,
  treatment = NULL,
  source,
  covariates,
  treatment_level = NULL,
  analysis = c("transport", "joint"),
  estimator = c("dr", "om", "ipw", "naive"),
  selection_model = NULL,
  propensity_model = NULL,
  outcome_model = NULL,
  outcome_type = c("auto", "binary", "continuous"),
  se_method = c("bootstrap", "influence", "none"),
  n_boot = 500,
  conf_level = 0.95,
  boot_ci_type = c("percentile", "normal", "basic"),
  stratified_boot = TRUE,
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

  Numeric vector of treatment indicators (0/1), or `NULL` for factual
  prediction model transportability (no treatment/intervention). When
  `NULL`, only the selection model is used for weighting.

- source:

  Numeric vector of population indicators (1=source/RCT, 0=target).

- covariates:

  A matrix or data frame of baseline covariates.

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

- outcome_type:

  Character string specifying the outcome type:

  - `"auto"`: Auto-detect from data (default)

  - `"binary"`: Binary outcome (0/1) - uses efficient transformation

  - `"continuous"`: Continuous outcome - models loss directly

- se_method:

  Method for standard error estimation:

  - `"bootstrap"`: Bootstrap standard errors (default)

  - `"influence"`: Influence function-based standard errors

  - `"none"`: No standard error estimation

- n_boot:

  Number of bootstrap replications (default: 500).

- conf_level:

  Confidence level for intervals (default: 0.95).

- boot_ci_type:

  Type of bootstrap confidence interval to compute:

  - `"percentile"`: Percentile method (default)

  - `"normal"`: Normal approximation using bootstrap SE

  - `"basic"`: Basic bootstrap interval

- stratified_boot:

  Logical indicating whether to use stratified bootstrap that preserves
  the source/target ratio (default: TRUE). Recommended for
  transportability analysis.

- cross_fit:

  Logical indicating whether to use cross-fitting for nuisance model
  estimation (default: FALSE).

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

An object of class `c("tr_mse", "tr_performance")` containing:

- estimate:

  Point estimate of transportable MSE

- se:

  Standard error (if computed)

- ci_lower:

  Lower confidence interval bound

- ci_upper:

  Upper confidence interval bound

- estimator:

  Estimator used

- analysis:

  Analysis type

- naive_estimate:

  Naive MSE for comparison

- n_target:

  Number of target observations

- n_source:

  Number of source observations

- treatment_level:

  Treatment level (NULL for factual mode)

## Details

This function implements estimators for transporting prediction model
performance from a source population to a target population. It supports
two modes:

### Counterfactual Mode (treatment provided)

When `treatment` is specified, estimates the counterfactual MSE under a
hypothetical intervention, EL(Y^a, g(X)) \| S=0. This requires:

- Selection model: P(S=0\|X)

- Propensity model in source: P(A=a\|X, S=1)

- Outcome model trained on source data

### Factual Mode (treatment = NULL)

When `treatment` is `NULL`, estimates the MSE of observed outcomes,
EL(Y, g(X)) \| S=0. This is appropriate for factual prediction model
transportability (no causal/counterfactual interpretation). Only
requires:

- Selection model: P(S=0\|X)

- Outcome model trained on source data

### Analysis Types

**Transportability Analysis** (`analysis = "transport"`): Uses outcome
data from the source population to estimate performance in the target
population.

**Joint Analysis** (`analysis = "joint"`): Pools source and target data
to estimate performance in the target population. More efficient when
both populations have outcome data.

For observational analysis (single population), use
[`cf_mse()`](https://boyercb.github.io/cfperformance/reference/cf_mse.md)
instead.

## References

Steingrimsson, J. A., et al. (2023). "Transporting a Prediction Model
for Use in a New Target Population." *American Journal of Epidemiology*,
192(2), 296-304.
[doi:10.1093/aje/kwac128](https://doi.org/10.1093/aje/kwac128)

Li, S., et al. (2023). "Efficient estimation of the expected prediction
error under covariate shift." *Biometrics*, 79(1), 295-307.
[doi:10.1111/biom.13583](https://doi.org/10.1111/biom.13583)

Voter, S. R., et al. (2025). "Transportability of machine learning-based
counterfactual prediction models with application to CASS." *Diagnostic
and Prognostic Research*, 9(4).
[doi:10.1186/s41512-025-00201-y](https://doi.org/10.1186/s41512-025-00201-y)

Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
"Estimating and evaluating counterfactual prediction models."
*Statistics in Medicine*, 44(23-24), e70287.
[doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

## See also

[`cf_mse()`](https://boyercb.github.io/cfperformance/reference/cf_mse.md),
[`tr_auc()`](https://boyercb.github.io/cfperformance/reference/tr_auc.md)

## Examples

``` r
# Generate example data with source (RCT) and target populations
set.seed(123)
n <- 1000
# Covariates
x <- rnorm(n)
# Source indicator (S=1 for RCT, S=0 for target)
s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
# Treatment (randomized in source, confounded in target)
a <- ifelse(s == 1,
            rbinom(n, 1, 0.5),  # Randomized in RCT
            rbinom(n, 1, plogis(-0.5 + 0.5 * x)))  # Confounded in target
# Outcome
y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
# Predictions from some model
pred <- plogis(-1 + 0.8 * x)

# Estimate transportable MSE
result <- tr_mse(
  predictions = pred,
  outcomes = y,
  treatment = a,
  source = s,
  covariates = data.frame(x = x),
  treatment_level = 0,
  analysis = "transport",
  estimator = "dr",
  se_method = "none"  # Skip SE for speed
)
print(result)
#> 
#> Counterfactual Transportable MSE Estimation
#> --------------------------------------------- 
#> Analysis: transport 
#> Estimator: dr 
#> Treatment level: 0 
#> N target: 376  | N source: 624 
#> 
#> Estimate: 0.1779
#> 
#> Naive estimate: 0.1577 
#> 

# Factual prediction model transportability (no treatment)
result_trad <- tr_mse(
  predictions = pred,
  outcomes = y,
  source = s,
  covariates = data.frame(x = x),
  analysis = "transport",
  estimator = "dr",
  se_method = "none"
)
```
