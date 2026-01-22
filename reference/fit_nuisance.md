# Fit Nuisance Models for Counterfactual Performance Estimation

Fits propensity score and outcome models for use in counterfactual
performance estimation.

## Usage

``` r
fit_nuisance(
  propensity_formula,
  outcome_formula,
  data,
  treatment_level = 0,
  propensity_method = c("glm", "gam"),
  outcome_method = c("glm", "gam"),
  ...
)
```

## Arguments

- propensity_formula:

  Formula for the propensity score model.

- outcome_formula:

  Formula for the outcome model.

- data:

  Data frame containing the variables.

- treatment_level:

  The counterfactual treatment level (default: 0).

- propensity_method:

  Method for propensity score estimation:

  - `"glm"`: Logistic regression (default)

  - `"gam"`: Generalized additive model (requires mgcv)

- outcome_method:

  Method for outcome model estimation:

  - `"glm"`: Logistic regression (default)

  - `"gam"`: Generalized additive model (requires mgcv)

- ...:

  Additional arguments passed to model fitting functions.

## Value

An object of class `cf_nuisance` containing:

- propensity:

  Fitted propensity score model

- outcome:

  Fitted outcome model

- treatment_level:

  Counterfactual treatment level

## Details

This function provides a convenient way to fit nuisance models that can
then be passed to
[`cf_mse()`](https://boyercb.github.io/cfperformance/reference/cf_mse.md),
[`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md),
or
[`cf_calibration()`](https://boyercb.github.io/cfperformance/reference/cf_calibration.md).

For the outcome model, only observations with the counterfactual
treatment level are used for fitting, as required for the conditional
loss estimator.

## See also

[`cf_mse()`](https://boyercb.github.io/cfperformance/reference/cf_mse.md),
[`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md),
[`cf_calibration()`](https://boyercb.github.io/cfperformance/reference/cf_calibration.md)

## Examples

``` r
# Generate example data
set.seed(123)
n <- 500
df <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n)
)
df$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * df$x1))
df$y <- rbinom(n, 1, plogis(-1 + df$x1 + 0.5 * df$x2 - 0.5 * df$a))

# Fit nuisance models
nuisance <- fit_nuisance(
  propensity_formula = a ~ x1 + x2,
  outcome_formula = y ~ x1 + x2,
  data = df,
  treatment_level = 0
)

print(nuisance)
#> 
#> Nuisance Models for Counterfactual Performance
#> --------------------------------------------- 
#> Treatment level: 0 
#> 
#> Propensity Score Model:
#>   Method: glm 
#>   Formula: a ~ x1 + x2 
#> 
#> Outcome Model:
#>   Method: glm 
#>   Formula: y ~ x1 + x2 
#>   Fitted on n = 323 observations with A = 0 
#> 
```
