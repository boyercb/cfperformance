# Estimate Counterfactual Area Under the ROC Curve

Estimates the area under the receiver operating characteristic curve
(AUC) of a prediction model under a hypothetical intervention where
treatment is set to a specific level.

## Usage

``` r
cf_auc(
  predictions,
  outcomes,
  treatment,
  covariates,
  treatment_level = 0,
  estimator = c("dr", "cl", "ipw", "naive"),
  propensity_model = NULL,
  outcome_model = NULL,
  se_method = c("bootstrap", "influence", "none"),
  n_boot = 500,
  conf_level = 0.95,
  cross_fit = FALSE,
  n_folds = 5,
  parallel = FALSE,
  ncores = NULL,
  ps_trim = NULL,
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

- estimator:

  Character string specifying the estimator:

  - `"naive"`: Naive estimator (biased)

  - `"cl"`: Conditional loss estimator

  - `"ipw"`: Inverse probability weighting estimator

  - `"dr"`: Doubly robust estimator (default)

- propensity_model:

  Optional fitted propensity score model. If NULL, a logistic regression
  model is fit using the covariates.

- outcome_model:

  Optional fitted outcome model. If NULL, a regression model is fit
  using the covariates among treated/untreated. For binary outcomes,
  this should be a model for E\[Y\|X,A\] (binomial family). For
  continuous outcomes, this should be a model for E\[L\|X,A\] (gaussian
  family).

- se_method:

  Method for standard error estimation:

  - `"bootstrap"`: Bootstrap standard errors (default)

  - `"influence"`: Influence function-based standard errors

  - `"none"`: No standard error estimation

- n_boot:

  Number of bootstrap replications (default: 500).

- conf_level:

  Confidence level for intervals (default: 0.95).

- cross_fit:

  Logical indicating whether to use cross-fitting for nuisance model
  estimation (default: FALSE). Cross-fitting enables valid inference
  when using flexible machine learning estimators.

- n_folds:

  Number of folds for cross-fitting (default: 5).

- parallel:

  Logical indicating whether to use parallel processing for bootstrap
  (default: FALSE).

- ncores:

  Number of cores for parallel processing (default: NULL, which uses all
  available cores minus one).

- ps_trim:

  Propensity score trimming specification. Controls how extreme
  propensity scores are handled. Can be:

  - `NULL` (default): Uses absolute bounds `c(0.01, 0.99)`

  - `"none"`: No trimming applied

  - `"quantile"`: Quantile-based trimming with default `c(0.01, 0.99)`

  - `"absolute"`: Explicit absolute bounds with default `c(0.01, 0.99)`

  - A numeric vector of length 2: `c(lower, upper)` absolute bounds

  - A single numeric: Symmetric bounds `c(x, 1-x)`

  - A list with `method` ("absolute"/"quantile"/"none") and `bounds`

- ...:

  Additional arguments passed to internal functions.

## Value

An object of class `c("cf_auc", "cf_performance")` containing:

- estimate:

  Point estimate of counterfactual AUC

- se:

  Standard error (if computed)

- ci_lower:

  Lower confidence interval bound

- ci_upper:

  Upper confidence interval bound

- estimator:

  Estimator used

- naive_estimate:

  Naive AUC for comparison

- n_obs:

  Number of observations

- treatment_level:

  Counterfactual treatment level

## Details

The counterfactual AUC is defined as the probability that a randomly
selected individual with the outcome under the counterfactual
intervention has a higher predicted risk than a randomly selected
individual without the outcome.

The function implements three estimators:

**Outcome Model (OM/CL) Estimator**: Weights concordant pairs by the
predicted probability of case/non-case status under the counterfactual.

**IPW Estimator**: Weights concordant pairs by the inverse probability
of treatment.

**Doubly Robust (DR) Estimator**: Combines OM and IPW for double
robustness. When `cross_fit = TRUE`, uses cross-fitting for valid
inference with flexible ML methods (see
[`ml_learner()`](https://boyercb.github.io/cfperformance/reference/ml_learner.md)).

## References

Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
"Estimating and evaluating counterfactual prediction models."
*Statistics in Medicine*, 44(23-24), e70287.
[doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

Li, B., Gatsonis, C., Dahabreh, I. J., & Steingrimsson, J. A. (2022).
"Estimating the area under the ROC curve when transporting a prediction
model to a target population." *Biometrics*, 79(3), 2343-2356.
[doi:10.1111/biom.13796](https://doi.org/10.1111/biom.13796)

## See also

[`cf_mse()`](https://boyercb.github.io/cfperformance/reference/cf_mse.md),
[`cf_calibration()`](https://boyercb.github.io/cfperformance/reference/cf_calibration.md),
[`ml_learner()`](https://boyercb.github.io/cfperformance/reference/ml_learner.md)

## Examples

``` r
# Generate example data
set.seed(123)
n <- 500
x <- rnorm(n)
a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
pred <- plogis(-1 + 0.8 * x)

# Estimate counterfactual AUC under no treatment
result <- cf_auc(
  predictions = pred,
  outcomes = y,
  treatment = a,
  covariates = data.frame(x = x),
  treatment_level = 0,
  estimator = "dr",
  se_method = "none"
)
print(result)
#> 
#> Counterfactual AUC Estimation
#> ---------------------------------------- 
#> Estimator: dr 
#> Treatment level: 0 
#> N observations: 500 
#> 
#> Estimate: 0.7331
#> 
#> Naive estimate: 0.7291 
#> 
```
