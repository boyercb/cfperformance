# Compute DR MSE with Cross-Fitting

Computes the doubly robust MSE estimator using cross-fitted nuisance
functions.

## Usage

``` r
.compute_mse_crossfit(
  predictions,
  outcomes,
  treatment,
  covariates,
  treatment_level,
  K = 5,
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

- K:

  Number of folds for cross-fitting.

- ...:

  Additional arguments passed to internal functions.

## Value

Doubly robust MSE estimate with cross-fitting.
