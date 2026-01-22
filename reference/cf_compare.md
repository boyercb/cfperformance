# Compare Multiple Prediction Models

Compares the counterfactual performance of multiple prediction models
using cross-validation or a held-out test set.

## Usage

``` r
cf_compare(
  models,
  data,
  treatment,
  treatment_level = 0,
  nuisance_covariates = NULL,
  metric = c("mse", "auc", "both"),
  estimator = c("dr", "cl", "ipw", "naive"),
  method = c("cv", "holdout"),
  K = 5,
  test_prop = 0.2,
  seed = NULL,
  ...
)
```

## Arguments

- models:

  A named list of model formulas or fitted model objects.

- data:

  A data frame containing all variables.

- treatment:

  Character string naming the treatment variable.

- treatment_level:

  The counterfactual treatment level (default: 0).

- nuisance_covariates:

  Character vector of covariate names for nuisance models. If NULL,
  inferred from model formulas.

- metric:

  Character string specifying the performance metric.

- estimator:

  Character string specifying the estimator.

- method:

  Comparison method: `"cv"` for cross-validation (default), or
  `"holdout"` for train/test split.

- K:

  Number of folds for CV (default: 5).

- test_prop:

  Proportion of data for test set if method = "holdout".

- seed:

  Random seed for reproducibility.

- ...:

  Additional arguments passed to cf_cv or internal functions.

## Value

An object of class `cf_compare` containing:

- results:

  Data frame with model performance comparisons

- best_model:

  Name of the best performing model

- metric:

  Performance metric used

- method:

  Comparison method used

## Examples

``` r
# Generate example data
set.seed(123)
n <- 300
data <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  x3 = rnorm(n)
)
data$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * data$x1))
data$y <- rbinom(n, 1, plogis(-1 + data$x1 + 0.5 * data$x2 - 0.3 * data$a))

# Compare models
models <- list(
  "Simple" = y ~ x1,
  "Full" = y ~ x1 + x2 + x3
)

comparison <- cf_compare(
  models = models,
  data = data,
  treatment = "a",
  metric = "mse",
  K = 3
)
print(comparison)
#> 
#> Counterfactual Model Comparison
#> --------------------------------------------- 
#> Method: cv (K = 3 )
#> Estimator: dr 
#> 
#>   model mse_mean mse_se mse_naive_mean
#>  Simple   0.1737 0.0179      0.1787998
#>    Full   0.1616 0.0028      0.1764572
#> 
#> Best model: Full 
#> 
```
