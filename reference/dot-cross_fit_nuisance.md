# Cross-Fitting for Nuisance Model Estimation

Implements K-fold cross-fitting for nuisance model estimation to enable
valid inference with flexible machine learning estimators.

## Usage

``` r
.cross_fit_nuisance(
  treatment,
  outcomes,
  covariates,
  treatment_level,
  predictions,
  K = 5,
  propensity_formula = NULL,
  outcome_formula = NULL,
  ...
)
```

## Arguments

- treatment:

  Numeric vector of treatment indicators.

- outcomes:

  Numeric vector of observed outcomes.

- covariates:

  Matrix or data frame of covariates.

- treatment_level:

  Counterfactual treatment level.

- predictions:

  Numeric vector of model predictions.

- K:

  Number of folds for cross-fitting (default: 5).

- propensity_formula:

  Formula for propensity model (optional).

- outcome_formula:

  Formula for outcome model (optional).

- ...:

  Additional arguments passed to model fitting functions.

## Value

List containing:

- ps:

  Cross-fitted propensity scores

- h:

  Cross-fitted conditional loss predictions

- folds:

  Fold assignments

## Details

Cross-fitting (sample splitting) allows the use of flexible machine
learning methods for nuisance function estimation while maintaining
valid inference. Each observation's nuisance function predictions are
made using models trained on data excluding that observation's fold.

## References

Chernozhukov, V., Chetverikov, D., Demirer, M., et al. (2018).
"Double/debiased machine learning for treatment and structural
parameters." *The Econometrics Journal*, 21(1), C1-C68.
