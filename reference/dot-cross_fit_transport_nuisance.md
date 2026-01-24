# Cross-fit nuisance models for transportability analysis

Implements K-fold cross-fitting for nuisance model estimation in
transportability settings.

## Usage

``` r
.cross_fit_transport_nuisance(
  treatment,
  outcomes,
  source,
  covariates,
  predictions,
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

- treatment:

  Numeric vector of treatment indicators.

- outcomes:

  Numeric vector of observed outcomes.

- source:

  Numeric vector indicating source (1) or target (0) population.

- covariates:

  Matrix or data frame of covariates.

- predictions:

  Numeric vector of model predictions.

- treatment_level:

  Counterfactual treatment level.

- analysis:

  Either "transport" or "joint".

- K:

  Number of folds for cross-fitting (default: 5).

- selection_learner:

  Optional ml_learner for selection model.

- propensity_learner:

  Optional ml_learner for propensity model.

- outcome_learner:

  Optional ml_learner for outcome model.

- outcome_type:

  Either "binary" or "continuous" (default: "binary").

- parallel:

  Logical for parallel processing.

- ps_trim_spec:

  Parsed propensity score trimming specification.

- ...:

  Additional arguments.

## Value

List containing cross-fitted nuisance function predictions.
