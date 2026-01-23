# Simulated Cardiovascular Disease Data

A simulated dataset for demonstrating counterfactual prediction model
performance estimation. The data mimics a scenario where patients
receive a cardiovascular treatment based on their risk factors, and we
want to evaluate how well a prediction model would perform under
different treatment policies.

## Usage

``` r
cvd_sim
```

## Format

A data frame with 1000 rows and 6 variables:

- age:

  Standardized age (mean 0, SD 1)

- bp:

  Standardized blood pressure (mean 0, SD 1)

- chol:

  Standardized cholesterol level (mean 0, SD 1)

- treatment:

  Binary treatment indicator (1 = treated, 0 = untreated). Treatment
  assignment is confounded by age and blood pressure.

- event:

  Binary outcome indicating cardiovascular event (1 = event, 0 = no
  event). Risk depends on age, bp, chol, and is reduced by treatment.

- risk_score:

  Predicted probability of event from a logistic regression model fit on
  observed data using age, bp, and chol as predictors.

## Source

Simulated data based on the framework in Boyer, Dahabreh & Steingrimsson
(2025). "Estimating and evaluating counterfactual prediction models."
Statistics in Medicine, 44(23-24), e70287.
[doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

## Details

The data generating process is:

- Covariates are independent standard normal

- Treatment probability: `plogis(-0.3 + 0.4*age + 0.3*bp + 0.1*chol)`

- Outcome probability:
  `plogis(-2 + 0.6*age + 0.5*bp + 0.3*chol - 0.4*treatment)`

This creates confounding because sicker patients (higher age, bp) are
more likely to receive treatment, but treatment also affects outcomes.

## Examples

``` r
data(cvd_sim)
head(cvd_sim)
#>          age          bp       chol treatment event risk_score
#> 1 -0.2078913 -0.43879526 -0.5697974         0     0 0.07548152
#> 2 -1.2517361  1.30171507  0.7798967         0     0 0.13491154
#> 3  1.7957878 -0.39076092 -0.1731313         1     0 0.18333022
#> 4 -1.2464064  0.08506276  0.0269594         1     0 0.07145644
#> 5 -0.5880067  0.10358176  0.8346190         1     0 0.11730651
#> 6 -0.9132198  0.88158838  0.6061392         0     0 0.12684642

# Estimate counterfactual MSE under no treatment
result <- cf_mse(
  predictions = cvd_sim$risk_score,
  outcomes = cvd_sim$event,
  treatment = cvd_sim$treatment,
  covariates = cvd_sim[, c("age", "bp", "chol")],
  treatment_level = 0,
  estimator = "dr"
)
result
#> 
#> Counterfactual MSE Estimation
#> ---------------------------------------- 
#> Estimator: dr 
#> Treatment level: 0 
#> N observations: 2500 
#> 
#> Estimate: 0.1186 (SE: 0.0064 )
#> 95% CI: [0.1059, 0.1322]
#> 
#> Naive estimate: 0.1086 
#> 
```
