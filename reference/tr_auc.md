# Estimate Transportable Area Under the ROC Curve

Estimates the area under the receiver operating characteristic curve
(AUC) of a prediction model in a target population using data
transported from a source population (typically an RCT).

## Usage

``` r
tr_auc(
  predictions,
  outcomes,
  treatment,
  source,
  covariates,
  treatment_level = 0,
  analysis = c("transport", "joint"),
  estimator = c("dr", "om", "ipw", "naive"),
  selection_model = NULL,
  propensity_model = NULL,
  outcome_model = NULL,
  se_method = c("bootstrap", "influence", "none"),
  n_boot = 500,
  conf_level = 0.95,
  stratified_boot = TRUE,
  cross_fit = FALSE,
  n_folds = 5,
  parallel = FALSE,
  ncores = NULL,
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

- source:

  Numeric vector of population indicators (1=source/RCT, 0=target).

- covariates:

  A matrix or data frame of baseline covariates.

- treatment_level:

  The treatment level of interest (default: 0).

- analysis:

  Character string specifying the type of analysis:

  - `"transport"`: Use source outcomes for target estimation (default)

  - `"joint"`: Pool source and target data

- estimator:

  Character string specifying the estimator:

  - `"naive"`: Naive estimator (biased)

  - `"om"`: Outcome model estimator

  - `"ipw"`: Inverse probability weighting estimator

  - `"dr"`: Doubly robust estimator (default)

- selection_model:

  Optional fitted selection model for P(S=0\|X). If NULL, a logistic
  regression model is fit using the covariates.

- propensity_model:

  Optional fitted propensity score model for P(A=1\|X,S=1). If NULL, a
  logistic regression model is fit using source data.

- outcome_model:

  Optional fitted outcome model for E\[L(Y,g)\|X,A,S\]. If NULL, a
  regression model is fit using the relevant data.

- se_method:

  Method for standard error estimation:

  - `"bootstrap"`: Bootstrap standard errors (default)

  - `"influence"`: Influence function-based standard errors

  - `"none"`: No standard error estimation

- n_boot:

  Number of bootstrap replications (default: 500).

- conf_level:

  Confidence level for intervals (default: 0.95).

- stratified_boot:

  Logical indicating whether to use stratified bootstrap that preserves
  the source/target ratio (default: TRUE). Recommended for
  transportability analysis.

- cross_fit:

  Logical indicating whether to use cross-fitting for nuisance model
  estimation (default: FALSE).

- n_folds:

  Number of folds for cross-fitting (default: 5).

- parallel:

  Logical indicating whether to use parallel processing for bootstrap
  (default: FALSE).

- ncores:

  Number of cores for parallel processing (default: NULL, which uses all
  available cores minus one).

- ...:

  Additional arguments passed to internal functions.

## Value

An object of class `c("tr_auc", "tr_performance")` containing:

- estimate:

  Point estimate of transportable AUC

- se:

  Standard error (if computed)

- ci_lower:

  Lower confidence interval bound

- ci_upper:

  Upper confidence interval bound

- estimator:

  Estimator used

- analysis:

  Analysis type

- naive_estimate:

  Naive AUC for comparison

- n_target:

  Number of target observations

- n_source:

  Number of source observations

- treatment_level:

  Treatment level

## Details

This function implements estimators for transporting prediction model
AUC from a source population (typically an RCT) to a target population.
The AUC is defined as the probability that a randomly selected case has
a higher predicted risk than a randomly selected non-case.

**Transportability Analysis**: Uses outcome data from the source/RCT
population to estimate AUC in the target population. Requires:

- Selection model: P(S=0\|X)

- Propensity model in source: P(A=1\|X, S=1)

- Outcome model trained on source data: E\[Y\|X, A, S=1\]

**Joint Analysis**: Pools source and target data to estimate AUC in the
target population. More efficient when both populations have outcome
data.

For observational analysis (single population), use
[`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md)
instead.

The estimators use U-statistic formulations that weight concordant pairs
by their estimated probabilities of occurring in the target population
under the counterfactual treatment.

## References

Voter, S. R., et al. (2025). "Transportability of machine learning-based
counterfactual prediction models with application to CASS." *Diagnostic
and Prognostic Research*, 9(4).
[doi:10.1186/s41512-025-00201-y](https://doi.org/10.1186/s41512-025-00201-y)

Li, B., Gatsonis, C., Dahabreh, I. J., & Steingrimsson, J. A. (2022).
"Estimating the area under the ROC curve when transporting a prediction
model to a target population." *Biometrics*.

Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
"Estimating and evaluating counterfactual prediction models."
*Statistics in Medicine*, 44(23-24), e70287.
[doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

## See also

[`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md),
[`tr_mse()`](https://boyercb.github.io/cfperformance/reference/tr_mse.md)

## Examples

``` r
# Generate example data with source (RCT) and target populations
set.seed(123)
n <- 1000
# Covariates
x <- rnorm(n)
# Source indicator (S=1 for RCT, S=0 for target)
s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
# Treatment (randomized in source, confounded in target)
a <- ifelse(s == 1,
            rbinom(n, 1, 0.5),  # Randomized in RCT
            rbinom(n, 1, plogis(-0.5 + 0.5 * x)))  # Confounded in target
# Outcome (binary)
y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
# Predictions from some model
pred <- plogis(-1 + 0.8 * x)

# Estimate transportable AUC
result <- tr_auc(
  predictions = pred,
  outcomes = y,
  treatment = a,
  source = s,
  covariates = data.frame(x = x),
  treatment_level = 0,
  analysis = "transport",
  estimator = "dr",
  se_method = "none"  # Skip SE for speed
)
print(result)
#> 
#> Transportable AUC Estimation
#> --------------------------------------------- 
#> Analysis: transport 
#> Estimator: dr 
#> Treatment level: 0 
#> N target: 376  | N source: 624 
#> 
#> Estimate: 0.7051
#> 
#> Naive estimate: 0.7386 
#> 
```
