# Compute transportable MSE with cross-fitting

Computes the DR transportable MSE estimator using cross-fitted nuisance
functions.

## Usage

``` r
.compute_tr_mse_crossfit(
  predictions,
  outcomes,
  treatment,
  source,
  covariates,
  treatment_level,
  analysis,
  K = 5,
  selection_learner = NULL,
  propensity_learner = NULL,
  outcome_learner = NULL,
  outcome_type = "binary",
  parallel = FALSE,
  ps_trim_spec = NULL,
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

- K:

  Number of folds for cross-fitting.

- selection_learner:

  Optional ml_learner for selection model.

- propensity_learner:

  Optional ml_learner for propensity model.

- outcome_learner:

  Optional ml_learner for outcome model.

- outcome_type:

  Either "binary" or "continuous".

- parallel:

  Logical for parallel processing.

- ...:

  Additional arguments.

## Value

List with estimate, SE, and cross-fitted nuisance values.
