# Estimate (Counterfactual) Calibration in the Target Population

Estimates the calibration of a prediction model in a target population
using data transported from a source population (typically an RCT).

## Usage

``` r
tr_calibration(
  predictions,
  outcomes,
  treatment,
  source,
  covariates,
  treatment_level = 0,
  analysis = c("transport", "joint"),
  estimator = c("ipw", "om"),
  selection_model = NULL,
  propensity_model = NULL,
  outcome_model = NULL,
  smoother = c("loess", "binned"),
  n_bins = 10,
  span = 0.75,
  se_method = c("bootstrap", "none"),
  n_boot = 500,
  conf_level = 0.95,
  stratified_boot = TRUE,
  ...
)
```

## Arguments

- predictions:

  Numeric vector of model predictions (typically probabilities).

- outcomes:

  Numeric vector of observed outcomes (must be binary 0/1).

- treatment:

  Numeric vector of treatment indicators (0/1).

- source:

  Numeric vector of population indicators (1=source/RCT, 0=target).

- covariates:

  A matrix or data frame of baseline covariates.

- treatment_level:

  The treatment level of interest (default: 0).

- analysis:

  Character string specifying the type of analysis:

  - `"transport"`: Use source outcomes for target estimation (default)

  - `"joint"`: Pool source and target data

- estimator:

  Character string specifying the estimator:

  - `"ipw"`: Inverse probability weighting estimator (default)

  - `"om"`: Outcome model estimator

- selection_model:

  Optional fitted selection model for P(S=0\|X). If NULL, a logistic
  regression model is fit using the covariates.

- propensity_model:

  Optional fitted propensity score model for P(A=1\|X,S=1). If NULL, a
  logistic regression model is fit using source data.

- outcome_model:

  Optional fitted outcome model for E\[Y\|X,A,S\]. If NULL, a regression
  model is fit using the relevant data.

- smoother:

  Smoothing method for the calibration curve:

  - `"loess"`: Local polynomial regression (default)

  - `"binned"`: Binned calibration

- n_bins:

  Number of bins for binned calibration (default: 10).

- span:

  Span parameter for LOESS smoothing (default: 0.75).

- se_method:

  Method for standard error estimation:

  - `"bootstrap"`: Bootstrap standard errors (default)

  - `"none"`: No standard error estimation

- n_boot:

  Number of bootstrap replications (default: 500).

- conf_level:

  Confidence level for intervals (default: 0.95).

- stratified_boot:

  Logical indicating whether to use stratified bootstrap that preserves
  the source/target ratio (default: TRUE).

- ...:

  Additional arguments passed to internal functions.

## Value

An object of class `c("tr_calibration", "tr_performance")` containing:

- predicted:

  Vector of predicted probabilities

- observed:

  Vector of smoothed observed probabilities

- smoother:

  Smoothing method used

- ici:

  Integrated calibration index

- e50:

  Median absolute calibration error

- e90:

  90th percentile absolute calibration error

- emax:

  Maximum absolute calibration error

- se:

  Named list of standard errors for calibration metrics (if computed)

- ci_lower:

  Named list of lower confidence bounds

- ci_upper:

  Named list of upper confidence bounds

- estimator:

  Estimator used

- analysis:

  Analysis type

- n_target:

  Number of target observations

- n_source:

  Number of source observations

## Details

This function implements estimators for transporting prediction model
calibration from a source population (typically an RCT) to a target
population.

**Transportability Analysis**: Uses outcome data from the source/RCT
population to estimate calibration in the target population. The IPW
estimator weights source observations to represent the target
population.

**Joint Analysis**: Pools source and target data to estimate calibration
in the target population.

The function computes several calibration metrics:

- **ICI** (Integrated Calibration Index): Mean absolute difference
  between predicted and observed probabilities

- **E50**: Median absolute calibration error

- **E90**: 90th percentile absolute calibration error

- **Emax**: Maximum absolute calibration error

For observational analysis (single population), use
[`cf_calibration()`](https://boyercb.github.io/cfperformance/reference/cf_calibration.md)
instead.

## References

Voter, S. R., et al. (2025). "Transportability of machine learning-based
counterfactual prediction models with application to CASS." *Diagnostic
and Prognostic Research*, 9(4).
[doi:10.1186/s41512-025-00201-y](https://doi.org/10.1186/s41512-025-00201-y)

Austin, P. C., & Steyerberg, E. W. (2019). "The Integrated Calibration
Index (ICI) and related metrics for quantifying the calibration of
logistic regression models." *Statistics in Medicine*, 38(21),
4051-4065.

## See also

[`cf_calibration()`](https://boyercb.github.io/cfperformance/reference/cf_calibration.md),
[`tr_mse()`](https://boyercb.github.io/cfperformance/reference/tr_mse.md),
[`tr_auc()`](https://boyercb.github.io/cfperformance/reference/tr_auc.md),
[`plot.tr_calibration()`](https://boyercb.github.io/cfperformance/reference/plot.tr_calibration.md)

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
# Outcome (binary)
y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
# Prediction model
pred <- plogis(-1 + 0.8 * x)

# Estimate transportable calibration with IPW
result <- tr_calibration(
  predictions = pred,
  outcomes = y,
  treatment = a,
  source = s,
  covariates = data.frame(x = x),
  treatment_level = 0,
  analysis = "transport",
  estimator = "ipw",
  se_method = "none"  # Use "bootstrap" for SEs
)
print(result)
#> 
#> Transportable CALIBRATION Estimation
#> --------------------------------------------- 
#> Analysis: transport 
#> Estimator: ipw 
#> Treatment level: 0 
#> N target: 376  | N source: 624 
#> 
#> Calibration Metrics:
#>   ICI (Integrated Calibration Index): 0.0493 
#>   E50 (Median absolute error): 0.0551 
#>   E90 (90th percentile error): 0.096 
#>   Emax (Maximum error): 0.1107 
#> 
```
