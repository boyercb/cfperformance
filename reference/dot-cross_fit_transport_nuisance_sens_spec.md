# Cross-fit nuisance models for transportability sensitivity/specificity

Implements K-fold cross-fitting for nuisance model estimation in
transportability settings for sensitivity and specificity estimation.

## Usage

``` r
.cross_fit_transport_nuisance_sens_spec(
  treatment,
  outcomes,
  source,
  covariates,
  treatment_level,
  analysis,
  K = 5,
  selection_learner = NULL,
  propensity_learner = NULL,
  outcome_learner = NULL,
  ps_trim_spec = NULL,
  parallel = FALSE,
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

- ps_trim_spec:

  Parsed propensity score trimming specification.

- parallel:

  Logical for parallel processing.

- ...:

  Additional arguments.

## Value

List containing cross-fitted nuisance function predictions.
