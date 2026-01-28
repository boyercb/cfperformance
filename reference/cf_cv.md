# Cross-Validation with Counterfactual Performance Metrics

Performs K-fold cross-validation to estimate out-of-sample
counterfactual model performance. This function trains and evaluates
prediction models while properly accounting for treatment effects.

## Usage

``` r
cf_cv(
  formula,
  data,
  treatment,
  treatment_level = 0,
  nuisance_covariates = NULL,
  metric = c("mse", "auc", "both"),
  estimator = c("dr", "cl", "om", "ipw", "naive"),
  K = 5,
  repeats = 1,
  stratify = TRUE,
  seed = NULL,
  ...
)
```

## Arguments

- formula:

  A formula specifying the prediction model (e.g., `Y ~ X1 + X2`).

- data:

  A data frame containing the variables in the formula plus `treatment`
  and any additional covariates for nuisance models.

- treatment:

  Character string naming the treatment variable in `data`.

- treatment_level:

  The counterfactual treatment level (default: 0).

- nuisance_covariates:

  Character vector of covariate names for nuisance models. If NULL, uses
  all predictors from the formula.

- metric:

  Character string specifying the performance metric: `"mse"` (default),
  `"auc"`, or `"both"`.

- estimator:

  Character string specifying the estimator: `"dr"` (default), `"cl"`
  (conditional loss, for MSE), `"om"` (outcome model, for AUC), `"ipw"`,
  or `"naive"`. Note: `"cl"` is automatically mapped to `"om"` for AUC
  metrics.

- K:

  Number of folds (default: 5).

- repeats:

  Number of times to repeat K-fold CV (default: 1).

- stratify:

  Logical indicating whether to stratify folds by outcome (default: TRUE
  for binary outcomes).

- seed:

  Random seed for reproducibility (default: NULL).

- ...:

  Additional arguments passed to internal functions.

## Value

An object of class `cf_cv` containing:

- results:

  Data frame with fold-level performance estimates

- summary:

  Summary statistics across folds

- metric:

  Performance metric used

- estimator:

  Estimator used

- K:

  Number of folds

- repeats:

  Number of repeats

- call:

  The matched call

## Details

Cross-validation for counterfactual prediction models requires special
care:

1.  **Nuisance model estimation**: Propensity and outcome models are
    re-fit in each training fold to avoid overfitting.

2.  **Sample splitting**: The prediction model is trained on the
    training fold and evaluated on the test fold using counterfactual
    estimators.

3.  **Stratification**: For binary outcomes, stratified sampling ensures
    each fold has similar outcome prevalence.

## References

Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
"Estimating and evaluating counterfactual prediction models."
*Statistics in Medicine*, 44(23-24), e70287.
[doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

## See also

[`cf_mse()`](https://boyercb.github.io/cfperformance/reference/cf_mse.md),
[`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md),
[`cf_compare()`](https://boyercb.github.io/cfperformance/reference/cf_compare.md)

## Examples

``` r
# Generate example data
set.seed(123)
n <- 300
data <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n)
)
data$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * data$x1))
data$y <- rbinom(n, 1, plogis(-1 + data$x1 + 0.5 * data$x2 - 0.3 * data$a))

# 5-fold cross-validation
cv_result <- cf_cv(
  formula = y ~ x1 + x2,
  data = data,
  treatment = "a",
  treatment_level = 0,
  metric = "mse",
  K = 5
)
print(cv_result)
#> 
#> Counterfactual Cross-Validation Results
#> --------------------------------------------- 
#> Folds: 5
#> Estimator: dr 
#> Treatment level: 0 
#> 
#> MSE:
#>   Counterfactual: 0.1568 (SE: 0.0117 )
#>   Naive:          0.1668 
#> 
```
