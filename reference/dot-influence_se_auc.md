# Influence Function-Based Standard Errors for AUC

Computes standard errors for counterfactual AUC using influence
functions.

## Usage

``` r
.influence_se_auc(
  predictions,
  outcomes,
  treatment,
  covariates,
  treatment_level,
  estimator,
  propensity_model,
  outcome_model
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

- propensity_model:

  Fitted propensity model.

- outcome_model:

  Fitted outcome model.

## Value

Standard error estimate.

## Details

The influence function for the AUC is based on the U-statistic
representation of the concordance probability. See Li et al. (2024) for
details on the influence function-based variance estimator for
counterfactual AUC.

## References

Li, B., Steingrimsson, J. A., & Dahabreh, I. J. (2024). "Efficient
inference for counterfactual area under the ROC curve."
