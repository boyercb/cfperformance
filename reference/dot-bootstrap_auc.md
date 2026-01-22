# Bootstrap Standard Errors for AUC

Computes bootstrap standard errors and confidence intervals for
counterfactual AUC estimators.

## Usage

``` r
.bootstrap_auc(
  predictions,
  outcomes,
  treatment,
  covariates,
  treatment_level,
  estimator,
  n_boot = 500,
  conf_level = 0.95,
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

List containing bootstrap results.
