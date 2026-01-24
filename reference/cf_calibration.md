# Estimate Counterfactual Calibration Curve

Estimates the calibration curve of a prediction model under a
hypothetical intervention where treatment is set to a specific level.

## Usage

``` r
cf_calibration(
  predictions,
  outcomes,
  treatment,
  covariates,
  treatment_level = 0,
  estimator = c("dr", "ipw", "cl"),
  propensity_model = NULL,
  outcome_model = NULL,
  smoother = c("loess", "binned"),
  n_bins = 10,
  span = 0.75,
  se_method = c("none", "bootstrap"),
  n_boot = 200,
  conf_level = 0.95,
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

  - `"none"`: No standard errors (default, fastest)

  - `"bootstrap"`: Bootstrap standard errors

- n_boot:

  Number of bootstrap replications (default: 200).

- conf_level:

  Confidence level for intervals (default: 0.95).

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

An object of class `c("cf_calibration", "cf_performance")` containing:

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

  List of standard errors (if se_method = "bootstrap")

- ci_lower:

  List of lower CI bounds (if se_method = "bootstrap")

- ci_upper:

  List of upper CI bounds (if se_method = "bootstrap")

- boot_curves:

  Bootstrap calibration curves for CI bands (if se_method = "bootstrap")

## Details

The counterfactual calibration curve estimates the relationship between
predicted risk and observed risk under the counterfactual intervention.

The function implements three estimators:

**IPW Estimator**: Weights observations by the inverse probability of
receiving the counterfactual treatment. Requires a correctly specified
propensity score model.

**Conditional Loss (CL) Estimator**: Uses the fitted outcome model
\\E\[Y \| X, A=a\]\\ to estimate calibration over all observations.
Requires a correctly specified outcome model.

**Doubly Robust (DR) Estimator**: Combines CL and IPW approaches.
Consistent if either the propensity or outcome model is correctly
specified.

## References

Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
"Estimating and evaluating counterfactual prediction models."
*Statistics in Medicine*, 44(23-24), e70287.
[doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

Steingrimsson, J. A., Gatsonis, C., Li, B., & Dahabreh, I. J. (2023).
"Transporting a prediction model for use in a new target population."
*American Journal of Epidemiology*, 192(2), 296-304.

## See also

[`cf_mse()`](https://boyercb.github.io/cfperformance/reference/cf_mse.md),
[`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md),
[`plot.cf_calibration()`](https://boyercb.github.io/cfperformance/reference/plot.cf_calibration.md)

## Examples

``` r
# Generate example data
set.seed(123)
n <- 500
x <- rnorm(n)
a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
pred <- plogis(-1 + 0.8 * x)

# Estimate counterfactual calibration curve with different estimators
result_ipw <- cf_calibration(
  predictions = pred,
  outcomes = y,
  treatment = a,
  covariates = data.frame(x = x),
  treatment_level = 0,
  estimator = "ipw"
)

result_dr <- cf_calibration(
  predictions = pred,
  outcomes = y,
  treatment = a,
  covariates = data.frame(x = x),
  treatment_level = 0,
  estimator = "dr"
)
print(result_dr)
#> 
#> Counterfactual CALIBRATION Estimation
#> ---------------------------------------- 
#> Estimator: dr 
#> Treatment level: 0 
#> N observations: 500 
#> 
#> Calibration Metrics:
#>   ICI (Integrated Calibration Index): 0.0444 
#>   E50 (Median absolute error): 0.0184 
#>   E90 (90th percentile error): 0.0996 
#>   Emax (Maximum error): 0.3784 
#> 
# plot(result_dr)  # If ggplot2 is available

# With bootstrap confidence bands
# result_boot <- cf_calibration(
#   predictions = pred, outcomes = y, treatment = a,
#   covariates = data.frame(x = x), treatment_level = 0,
#   se_method = "bootstrap", n_boot = 200
# )
# plot(result_boot)  # Shows confidence bands
```
