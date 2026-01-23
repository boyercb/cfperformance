# Simulated Transportability Data

A simulated dataset for demonstrating transportability analysis of
prediction model performance. The data includes source (RCT) and target
populations, where treatment is randomized in the source but confounded
in the target.

## Usage

``` r
transport_sim
```

## Format

A data frame with 1500 rows and 7 variables:

- age:

  Standardized age (mean 0, SD 1)

- biomarker:

  Continuous biomarker value

- smoking:

  Binary smoking status (1 = smoker, 0 = non-smoker)

- source:

  Population indicator (1 = source/RCT, 0 = target). RCT patients tend
  to be younger and less likely to smoke.

- treatment:

  Binary treatment indicator (1 = treated, 0 = untreated). Randomized
  (50/50) in source, confounded by age and biomarker in target.

- event:

  Binary outcome indicating event (1 = event, 0 = no event). Risk
  depends on age, biomarker, smoking, and is reduced by treatment.

- risk_score:

  Predicted probability of event from a model trained on the source
  population, using age, biomarker, and smoking as predictors.

## Source

Simulated data based on the framework in Voter et al. (2025).
"Transportability of machine learning-based counterfactual prediction
models." Diagnostic and Prognostic Research, 9(4).
[doi:10.1186/s41512-025-00201-y](https://doi.org/10.1186/s41512-025-00201-y)

## Details

The data generating process creates realistic heterogeneity between
source (RCT) and target populations:

- **Selection into RCT**: RCT patients are younger on average and less
  likely to be smokers, reflecting typical trial enrollment patterns.

- **Treatment assignment**:

  - Source: Randomized with `P(A=1) = 0.5`

  - Target: Confounded with
    `P(A=1|X) = plogis(-0.3 + 0.015*age + 0.2*biomarker)`

- **Outcome model**:
  `P(Y=1|X,A) = plogis(-2.5 + 0.03*age + 0.4*biomarker + 0.3*smoking - 0.5*A)`

This creates a scenario where naive performance estimates from the RCT
will not generalize to the target population due to covariate shift.

## Examples

``` r
data(transport_sim)
head(transport_sim)
#>           age  biomarker smoking source treatment event risk_score
#> 1  0.61691842 -0.5700952       0      0         1     0  0.2603127
#> 2  0.04053456  0.3601832       0      1         1     1  0.2962904
#> 3  0.76704199  0.4024288       1      1         0     1  0.3855963
#> 4  1.25892520 -1.5801491       0      0         1     1  0.2253197
#> 5  0.37086452  0.1336976       0      0         0     1  0.2972428
#> 6 -0.15499910  0.5227082       0      1         0     0  0.2978064

# Population sizes
table(transport_sim$source)  # 0=target, 1=source
#> 
#>   0   1 
#> 698 802 

# Estimate transportable MSE under no treatment
result <- tr_mse(
  predictions = transport_sim$risk_score,
  outcomes = transport_sim$event,
  treatment = transport_sim$treatment,
  source = transport_sim$source,
  covariates = transport_sim[, c("age", "biomarker", "smoking")],
  treatment_level = 0,
  analysis = "transport",
  estimator = "dr"
)
result
#> 
#> Transportable MSE Estimation
#> --------------------------------------------- 
#> Analysis: transport 
#> Estimator: dr 
#> Treatment level: 0 
#> N target: 698  | N source: 802 
#> 
#> Estimate: 0.2432 (SE: 0.0102 )
#> 95% CI: [0.2232, 0.263]
#> 
#> Naive estimate: 0.227 
#> 
```
