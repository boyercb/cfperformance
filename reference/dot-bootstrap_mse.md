# Bootstrap Standard Errors for MSE

Computes bootstrap standard errors and confidence intervals for
counterfactual MSE estimators.

## Usage

``` r
.bootstrap_mse(
  predictions,
  outcomes,
  treatment,
  covariates,
  treatment_level,
  estimator,
  n_boot = 500,
  conf_level = 0.95,
  boot_ci_type = c("percentile", "normal", "basic"),
  parallel = FALSE,
  ncores = NULL,
  ps_trim = NULL,
  outcome_type = "binary",
  propensity_model = NULL,
  outcome_model = NULL,
  point_estimate = NULL,
  ...
)
```

## Arguments

- predictions:

  Numeric vector of model predictions.

- outcomes:

  Numeric vector of observed outcomes.

- treatment:

  Numeric vector of treatment indicators.

- covariates:

  Matrix or data frame of covariates.

- treatment_level:

  Counterfactual treatment level.

- estimator:

  Character string specifying estimator type.

- n_boot:

  Number of bootstrap replications.

- conf_level:

  Confidence level.

- parallel:

  Logical indicating whether to use parallel processing.

- ncores:

  Number of cores for parallel processing.

- ...:

  Additional arguments.

## Value

List containing:

- se:

  Bootstrap standard error

- ci_lower:

  Lower confidence interval bound

- ci_upper:

  Upper confidence interval bound

- boot_estimates:

  Vector of bootstrap estimates
