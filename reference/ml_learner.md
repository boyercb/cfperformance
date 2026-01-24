# Specify a Machine Learning Learner for Nuisance Models

Creates a learner specification that can be passed to `propensity_model`
or `outcome_model` arguments in
[`cf_mse()`](https://boyercb.github.io/cfperformance/reference/cf_mse.md),
[`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md),
[`cf_calibration()`](https://boyercb.github.io/cfperformance/reference/cf_calibration.md),
and their transportability variants. When an `ml_learner` specification
is provided, cross-fitting is automatically used for valid inference.

## Usage

``` r
ml_learner(
  method = c("ranger", "xgboost", "grf", "glmnet", "superlearner", "custom"),
  ...,
  fit_fun = NULL,
  predict_fun = NULL
)
```

## Arguments

- method:

  Character string specifying the learner type:

  - `"ranger"`: Random forest via
    [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md)

  - `"xgboost"`: Gradient boosting via
    [`xgboost::xgboost()`](https://rdrr.io/pkg/xgboost/man/xgboost.html)

  - `"grf"`: Generalized random forest via
    [`grf::regression_forest()`](https://rdrr.io/pkg/grf/man/regression_forest.html)
    or
    [`grf::probability_forest()`](https://rdrr.io/pkg/grf/man/probability_forest.html)

  - `"glmnet"`: Regularized regression via
    [`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html)

  - `"superlearner"`: Ensemble via
    [`SuperLearner::SuperLearner()`](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)

  - `"custom"`: User-supplied fit/predict functions

- ...:

  Additional arguments passed to the fitting function.

- fit_fun:

  For `method = "custom"`, a function with signature
  `function(formula, data, family, ...)` that returns a fitted model
  object.

- predict_fun:

  For `method = "custom"`, a function with signature
  `function(object, newdata, ...)` that returns predicted probabilities.

## Value

An object of class `ml_learner` containing the learner specification.

## Details

### Supported Learners

**ranger**: Fast random forest implementation. Key arguments:

- `num.trees`: Number of trees (default: 500)

- `mtry`: Number of variables to sample at each split

- `min.node.size`: Minimum node size

**xgboost**: Gradient boosting. Key arguments:

- `nrounds`: Number of boosting rounds (default: 100)

- `max_depth`: Maximum tree depth (default: 6)

- `eta`: Learning rate (default: 0.3)

**grf**: Generalized random forests with built-in honesty. Key
arguments:

- `num.trees`: Number of trees (default: 2000)

- `honesty`: Whether to use honest estimation (default: TRUE)

**glmnet**: Elastic net regularization with cross-validation. Key
arguments:

- `alpha`: Elastic net mixing parameter (0 = ridge, 1 = lasso, default:
  1)

- `nfolds`: Number of CV folds for lambda selection (default: 10)

**superlearner**: Ensemble of multiple learners. Key arguments:

- `SL.library`: Vector of learner names (default:
  `c("SL.glm", "SL.ranger")`)

## See also

[`cf_mse()`](https://boyercb.github.io/cfperformance/reference/cf_mse.md),
[`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md),
[`cf_calibration()`](https://boyercb.github.io/cfperformance/reference/cf_calibration.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Random forest for propensity score
cf_mse(
  Y = outcome, A = treatment, predictions = preds, data = df,
  propensity_formula = treatment ~ x1 + x2,
  propensity_model = ml_learner("ranger", num.trees = 500),
  cross_fit = TRUE
)

# XGBoost with custom parameters
cf_mse(
  Y = outcome, A = treatment, predictions = preds, data = df,
  propensity_formula = treatment ~ x1 + x2,
  propensity_model = ml_learner("xgboost", nrounds = 200, max_depth = 4),
  cross_fit = TRUE
)

# Custom learner
my_fit <- function(formula, data, family, ...) {
  glm(formula, data = data, family = binomial())
}
my_predict <- function(object, newdata, ...) {
  predict(object, newdata = newdata, type = "response")
}
cf_mse(
  ...,
  propensity_model = ml_learner("custom", fit_fun = my_fit,
                                 predict_fun = my_predict)
)
} # }
```
