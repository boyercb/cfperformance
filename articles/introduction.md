# Introduction to cfperformance

## Overview

The `cfperformance` package provides methods for estimating prediction
model performance under hypothetical (counterfactual) interventions.
This is essential when:

1.  **Prediction models will be deployed in settings where treatment
    policies differ from training** - A model trained on patients who
    received a mixture of treatments may perform differently when
    deployed where everyone receives a specific treatment.

2.  **Predictions support treatment decisions** - When predictions
    inform who should receive treatment, naive performance estimates
    conflate model accuracy with treatment effects.

The methods implemented here are based on Boyer, Dahabreh &
Steingrimsson (2025), “Estimating and evaluating counterfactual
prediction models,” *Statistics in Medicine*, 44(23-24), e70287.
[doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

## Installation

``` r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("boyercb/cfperformance")
```

## Quick Start

``` r
# Load the included example dataset
data(cvd_sim)
head(cvd_sim)
#>          age          bp        chol treatment event risk_score
#> 1 -0.2078913  0.92673619  0.51441033         1     0 0.16456973
#> 2 -1.2517361 -0.07559820  1.85761302         1     0 0.06879905
#> 3  1.7957878  0.16615295 -1.09269567         0     0 0.23461394
#> 4 -1.2464064 -1.40349275  1.05642052         0     0 0.02765908
#> 5 -0.5880067 -0.08772247  1.18161254         1     0 0.08778086
#> 6 -0.9132198  0.16445812 -0.03630469         1     0 0.06799805
```

The `cvd_sim` dataset contains simulated cardiovascular data with: -
`age`, `bp`, `chol`: Patient covariates  
- `treatment`: Binary treatment indicator (confounded by covariates) -
`event`: Binary outcome (cardiovascular event) - `risk_score`:
Pre-computed predictions from a logistic regression model

### Estimating Counterfactual MSE

Now we can estimate how well the model would perform if everyone were
untreated (`treatment_level = 0`):

``` r
# Estimate MSE under counterfactual "no treatment" policy
mse_result <- cf_mse(
  predictions = cvd_sim$risk_score,
  outcomes = cvd_sim$event,
  treatment = cvd_sim$treatment,
  covariates = cvd_sim[, c("age", "bp", "chol")],
  treatment_level = 0,
  estimator = "dr"  # doubly robust estimator
)

mse_result
#> 
#> Counterfactual MSE Estimation
#> ---------------------------------------- 
#> Estimator: dr 
#> Treatment level: 0 
#> N observations: 1000 
#> 
#> Estimate: 0.1061 (SE: 0.0092 )
#> 95% CI: [0.0891, 0.125]
#> 
#> Naive estimate: 0.1061
```

The doubly robust estimator adjusts for confounding using both a
propensity score model and an outcome model, providing consistent
estimates even if one model is misspecified.

### Comparing Estimators

Let’s compare all available estimators:

``` r
estimators <- c("naive", "cl", "ipw", "dr")
results <- sapply(estimators, function(est) {
  cf_mse(
    predictions = cvd_sim$risk_score,
    outcomes = cvd_sim$event,
    treatment = cvd_sim$treatment,
    covariates = cvd_sim[, c("age", "bp", "chol")],
    treatment_level = 0,
    estimator = est
  )$estimate
})
names(results) <- estimators
round(results, 4)
#>  naive     cl    ipw     dr 
#> 0.1061 0.1055 0.1924 0.1061
```

- **naive**: Simply computes MSE on the subset with the target treatment
  level. Biased when treatment is confounded.
- **cl** (Conditional Loss): Models the outcome and integrates over the
  covariate distribution.
- **ipw** (Inverse Probability Weighting): Reweights observations to
  mimic the counterfactual population.
- **dr** (Doubly Robust): Combines outcome modeling and IPW; consistent
  if either model is correct.

### Estimating Counterfactual AUC

For discrimination (AUC), we can use similar methods:

``` r
auc_result <- cf_auc(
  predictions = cvd_sim$risk_score,
  outcomes = cvd_sim$event,
  treatment = cvd_sim$treatment,
  covariates = cvd_sim[, c("age", "bp", "chol")],
  treatment_level = 0,
  estimator = "dr"
)

auc_result
#> 
#> Counterfactual AUC Estimation
#> ---------------------------------------- 
#> Estimator: dr 
#> Treatment level: 0 
#> N observations: 1000 
#> 
#> Estimate: 0.6702 (SE: 0.0383 )
#> 95% CI: [0.5911, 0.7472]
#> 
#> Naive estimate: 0.7201
```

### Bootstrap Standard Errors

Both functions support bootstrap standard errors:

``` r
mse_with_se <- cf_mse(
  predictions = predictions,
  outcomes = outcome,
  treatment = treatment,
  covariates = data.frame(x1, x2),
  treatment_level = 0,
  estimator = "dr",
  se_method = "bootstrap",
  n_boot = 500
)
summary(mse_with_se)
```

## Calibration Curves

The package also supports counterfactual calibration assessment:

``` r
cal_result <- cf_calibration(
  predictions = predictions,
  outcomes = outcome,
  treatment = treatment,
  covariates = data.frame(x1, x2),
  treatment_level = 0
)

# Plot calibration curve
plot(cal_result)
```

## Cross-Validation for Model Selection

When comparing multiple prediction models, use counterfactual
cross-validation:

``` r
# Compare two models using counterfactual CV
models <- list(
  "Simple" = event ~ age,
  "Full" = event ~ age + bp + chol
)

comparison <- cf_compare(
  models = models,
  data = cvd_sim,
  treatment = "treatment",
  treatment_level = 0,
  metric = "mse",
  K = 5
)

comparison
#> 
#> Counterfactual Model Comparison
#> --------------------------------------------- 
#> Method: cv (K = 5 )
#> Estimator: dr 
#> 
#>   model mse_mean mse_se mse_naive_mean
#>  Simple   0.1085 0.0088      0.1121394
#>    Full   0.1070 0.0063      0.1072778
#> 
#> Best model: Full
```

## Key Concepts

### Why Counterfactual Performance?

Standard model performance evaluation computes metrics like MSE or AUC
on a test set. However, this answers: “How well does the model predict
outcomes *as they occurred*?”

When a model will be used to inform treatment decisions, we often need
to answer: “How well would the model predict outcomes *if everyone
received (or didn’t receive) treatment*?”

These can differ substantially when:

1.  Treatment is related to the outcome (treatment effects exist)
2.  Treatment is related to the covariates used for prediction
    (confounding)

### Assumptions

The methods in this package require:

1.  **Consistency**: Observed outcomes equal potential outcomes under
    the observed treatment.
2.  **Positivity**: All covariate patterns have positive probability of
    receiving each treatment level.
3.  **No unmeasured confounding**: Treatment is independent of potential
    outcomes given measured covariates.

These are standard causal inference assumptions. The package provides
warnings when positivity may be violated (extreme propensity scores).

### Choosing an Estimator

- Use **doubly robust (dr)** as the default - it’s consistent if either
  the propensity or outcome model is correct.
- Use **ipw** when you trust your propensity model but not your outcome
  model.
- Use **cl** when you trust your outcome model but not your propensity
  model.
- Use **naive** only as a baseline comparison.

## Further Reading

- Boyer CB, Dahabreh IJ, Steingrimsson JA. Counterfactual prediction
  model performance. *Statistics in Medicine*. 2025; 44(23-24):e70287.
  [doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

- Dahabreh IJ, Robertson SE, Steingrimsson JA. Extending inferences from
  a randomized trial to a new target population. *Statistics in
  Medicine*. 2020.

- Bang H, Robins JM. Doubly robust estimation in missing data and causal
  inference models. *Biometrics*. 2005.
