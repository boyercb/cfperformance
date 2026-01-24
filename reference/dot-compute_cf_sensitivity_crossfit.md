# Compute DR sensitivity with cross-fitted nuisance

Computes the doubly robust sensitivity estimator using cross-fitted
nuisance functions and returns both point estimate and influence
function-based standard error.

## Usage

``` r
.compute_cf_sensitivity_crossfit(
  predictions,
  outcomes,
  treatment,
  threshold,
  treatment_level,
  ps,
  m_hat,
  return_se = FALSE
)
```

## Arguments

- predictions:

  Numeric vector of model predictions.

- outcomes:

  Numeric vector of observed outcomes.

- treatment:

  Numeric vector of treatment indicators (0/1).

- threshold:

  Numeric vector of classification thresholds. Predictions above this
  value are classified as positive. Can be a single value or a vector
  for computing sensitivity at multiple thresholds simultaneously.
  Default is 0.5.

- treatment_level:

  The counterfactual treatment level (default: 0).

- ps:

  Cross-fitted propensity scores.

- m_hat:

  Cross-fitted outcome probabilities.

- return_se:

  Logical; if TRUE, returns list with estimate and SE.

## Value

If return_se=FALSE, numeric estimate. If TRUE, list with estimate and
se.
