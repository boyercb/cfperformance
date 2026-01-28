# Bootstrap for cross-fitted transportable AUC

Bootstrap for cross-fitted transportable AUC

## Usage

``` r
.bootstrap_tr_auc_crossfit(
  predictions,
  outcomes,
  treatment,
  source,
  covariates,
  treatment_level,
  analysis,
  n_folds,
  n_boot,
  conf_level,
  boot_ci_type = c("percentile", "normal", "basic"),
  stratified = TRUE,
  selection_learner = NULL,
  propensity_learner = NULL,
  outcome_learner = NULL,
  parallel,
  ncores,
  ps_trim_spec = NULL,
  point_estimate = NULL,
  ...
)
```

## Arguments

- n_folds:

  Number of folds for cross-fitting.

- selection_learner:

  Optional ml_learner for selection model.

- propensity_learner:

  Optional ml_learner for propensity model.

- outcome_learner:

  Optional ml_learner for outcome model.

## Value

List with se, ci_lower, ci_upper, boot_estimates.
