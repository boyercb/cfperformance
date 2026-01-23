# Compute DR AUC with Cross-Fitting

Computes the doubly robust AUC estimator using cross-fitted nuisance
functions.

## Usage

``` r
.compute_auc_crossfit(
  predictions,
  outcomes,
  treatment,
  covariates,
  treatment_level,
  K = 5,
  propensity_learner = NULL,
  outcome_learner = NULL,
  parallel = FALSE,
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

- propensity_learner:

  Optional ml_learner for propensity model.

- outcome_learner:

  Optional ml_learner for outcome model.

- parallel:

  Logical for parallel processing (not yet implemented).

- ...:

  Additional arguments.

## Value

Doubly robust AUC estimate with cross-fitting.

## References

Li, B., Gatsonis, C., Dahabreh, I. J., & Steingrimsson, J. A. (2022).
"Estimating the area under the ROC curve when transporting a prediction
model to a target population." *Biometrics*, 79(3), 2343-2356.
[doi:10.1111/biom.13796](https://doi.org/10.1111/biom.13796)
