# Bootstrap for cross-fitted sens/spec

Bootstrap for cross-fitted sens/spec

## Usage

``` r
.bootstrap_cf_sens_spec_crossfit(
  predictions,
  outcomes,
  treatment,
  covariates,
  threshold,
  treatment_level,
  metric,
  n_folds,
  n_boot,
  conf_level,
  propensity_learner = NULL,
  outcome_learner = NULL,
  parallel,
  ncores,
  ps_trim_spec = NULL
)
```
