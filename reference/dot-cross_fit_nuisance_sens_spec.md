# Cross-fit nuisance models for counterfactual sens/spec

Implements K-fold cross-fitting for nuisance model estimation.

## Usage

``` r
.cross_fit_nuisance_sens_spec(
  treatment,
  outcomes,
  covariates,
  treatment_level,
  K = 5,
  propensity_learner = NULL,
  outcome_learner = NULL,
  ps_trim_spec = NULL
)
```

## Arguments

- treatment:

  Numeric vector of treatment indicators (0/1).

- outcomes:

  Numeric vector of observed outcomes.

- covariates:

  A matrix or data frame of baseline covariates (confounders).

- treatment_level:

  The counterfactual treatment level (default: 0).

- K:

  Number of folds for cross-fitting.

- propensity_learner:

  Optional ml_learner for propensity model.

- outcome_learner:

  Optional ml_learner for outcome model.

## Value

List containing cross-fitted nuisance function predictions.
