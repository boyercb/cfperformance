# Estimate Counterfactual Mean Squared Error

Estimates the mean squared error (MSE) of a prediction model under a
hypothetical intervention where treatment is set to a specific level.

## Usage

``` r
cf_mse(
  predictions,
  outcomes,
  treatment,
  covariates,
  treatment_level = 0,
  estimator = c("dr", "cl", "ipw", "naive"),
  propensity_model = NULL,
  outcome_model = NULL,
  se_method = c("bootstrap", "influence", "none"),
  n_boot = 500,
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

An object of class `c("cf_mse", "cf_performance")` containing:

- estimate:

  Point estimate of counterfactual MSE

- se:

  Standard error (if computed)

- ci_lower:

  Lower confidence interval bound

- ci_upper:

  Upper confidence interval bound

- estimator:

  Estimator used

- naive_estimate:

  Naive MSE for comparison

- n_obs:

  Number of observations

- treatment_level:

  Counterfactual treatment level

## Details

The function implements four estimators for the counterfactual MSE:

**Naive Estimator**: Simply computes the empirical MSE using observed
outcomes. This is biased for the counterfactual estimand when treatment
affects outcomes.

**Conditional Loss (CL) Estimator**: Models the expected loss
conditional on covariates and treatment, then marginalizes. Requires a
correctly specified outcome model.

**IPW Estimator**: Weights observations by the inverse probability of
receiving the counterfactual treatment. Requires a correctly specified
propensity score model.

**Doubly Robust (DR) Estimator**: Combines CL and IPW approaches.
Consistent if either the propensity or outcome model is correctly
specified.

## References

Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
"Estimating and evaluating counterfactual prediction models."
*Statistics in Medicine*, 44(23-24), e70287.
[doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

## See also

[`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md),
[`cf_calibration()`](https://boyercb.github.io/cfperformance/reference/cf_calibration.md),
[`fit_nuisance()`](https://boyercb.github.io/cfperformance/reference/fit_nuisance.md)

## Examples

``` r
# Generate example data
set.seed(123)
n <- 500
x <- rnorm(n)
a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
pred <- plogis(-1 + 0.8 * x)  # Predictions from some model

# Estimate counterfactual MSE under no treatment
result <- cf_mse(
  predictions = pred,
  outcomes = y,
  treatment = a,
  covariates = data.frame(x = x),
  treatment_level = 0,
  estimator = "dr",
  se_method = "none"  # Skip SE for speed in example
)
print(result)
#> 
#> Counterfactual MSE Estimation
#> ---------------------------------------- 
#> Estimator: dr 
#> Treatment level: 0 
#> N observations: 500 
#> 
#> Estimate: 0.1754
#> 
#> Naive estimate: 0.1689 
#> 
```
