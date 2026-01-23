# Transportability Analysis with cfperformance

## Overview

This vignette demonstrates how to use the **transportability** functions
in `cfperformance` to evaluate prediction model performance when
transporting from a source population (such as a randomized controlled
trial) to a target population.

The methods are based on Voter et al. (2025), “Transportability of
machine learning-based counterfactual prediction models with application
to CASS,” *Diagnostic and Prognostic Research*, 9(4).
[doi:10.1186/s41512-025-00201-y](https://doi.org/10.1186/s41512-025-00201-y)

## When to Use Transportability Analysis

Transportability analysis is appropriate when:

1.  **You have a prediction model** trained in one population (e.g., an
    RCT)
2.  **You want to evaluate performance** in a different target
    population
3.  **Treatment is randomized** in the source population
4.  **Covariates are measured** in both populations

Common scenarios include:

- Evaluating an RCT-derived risk score in a real-world patient
  population
- Assessing whether a clinical prediction rule “transports” to a new
  setting
- Understanding how model performance varies across populations

## The Setting

Consider two populations:

- **Source (S=1)**: Often an RCT where treatment `A` is randomized
- **Target (S=0)**: The population where we want to deploy the model

We observe: - Covariates `X` in both populations - Treatment `A` in both
populations  
- Outcomes `Y` in the source (and possibly target) - Model predictions
`g(X)` for all individuals

The goal is to estimate the counterfactual prediction performance
$E\left\lbrack L\left( Y^{a},g(X) \right)|S = 0 \right\rbrack$ — how
well would the model perform in the target population if everyone
received treatment level `a`?

## Using the Included Example Data

The package includes a simulated transportability dataset:

``` r
data(transport_sim)
head(transport_sim)
#>           age   biomarker smoking source treatment event risk_score
#> 1  0.63916439 -0.06780891       1      1         0     0  0.3388230
#> 2  0.06014812  0.81999552       0      1         0     0  0.3309058
#> 3  0.78997359 -1.00222911       1      0         1     0  0.2773807
#> 4  1.28410328  0.78155960       1      0         1     0  0.4481046
#> 5  0.39198673  1.24382106       0      1         1     0  0.3841622
#> 6 -0.13627856  0.61176294       1      1         1     0  0.3482827

# Population sizes
cat("Source (RCT) n =", sum(transport_sim$source == 1), "\n")
#> Source (RCT) n = 1276
cat("Target n =", sum(transport_sim$source == 0), "\n")
#> Target n = 1224
```

The `transport_sim` dataset contains: - `age`, `biomarker`, `smoking`:
Patient covariates - `source`: Population indicator (1 = source/RCT, 0 =
target) - `treatment`: Binary treatment (randomized in source,
confounded in target) - `event`: Binary outcome - `risk_score`:
Predictions from a model trained in the source population

## Transport Analysis

The **transport analysis** uses outcomes from the source/RCT to estimate
performance in the target population. This is useful when:

- Outcomes are only observed in the RCT
- You want to leverage randomization in the source

### MSE in the Target Population

``` r
mse_result <- tr_mse(
  predictions = transport_sim$risk_score,
  outcomes = transport_sim$event,
  treatment = transport_sim$treatment,
  source = transport_sim$source,
  covariates = transport_sim[, c("age", "biomarker", "smoking")],
  treatment_level = 0,  # Evaluate under no treatment
  analysis = "transport",
  estimator = "dr",
  se_method = "none"
)

print(mse_result)
#> 
#> Transportable MSE Estimation
#> --------------------------------------------- 
#> Analysis: transport 
#> Estimator: dr 
#> Treatment level: 0 
#> N target: 1224  | N source: 1276 
#> 
#> Estimate: 0.2211
#> 
#> Naive estimate: 0.2178
```

### Comparing Estimators

Let’s compare all four estimators:

``` r
estimators <- c("naive", "om", "ipw", "dr")

mse_estimates <- sapply(estimators, function(est) {
  result <- tr_mse(
    predictions = transport_sim$risk_score,
    outcomes = transport_sim$event,
    treatment = transport_sim$treatment,
    source = transport_sim$source,
    covariates = transport_sim[, c("age", "biomarker", "smoking")],
    treatment_level = 0,
    analysis = "transport",
    estimator = est,
    se_method = "none"
  )
  result$estimate
})

data.frame(
  Estimator = estimators,
  MSE = round(mse_estimates, 4)
)
#>       Estimator    MSE
#> naive     naive 0.2178
#> om           om 0.2212
#> ipw         ipw 0.2213
#> dr           dr 0.2211
```

- **naive**: Uses source data with A=0 only (ignores population
  differences)
- **om** (outcome model): Models E\[L\|X,A,S=1\] and averages over
  target X
- **ipw**: Reweights source observations to match target population
- **dr** (doubly robust): Combines outcome modeling and IPW

### AUC in the Target Population

For discrimination, we can estimate the transportable AUC:

``` r
auc_result <- tr_auc(
  predictions = transport_sim$risk_score,
  outcomes = transport_sim$event,
  treatment = transport_sim$treatment,
  source = transport_sim$source,
  covariates = transport_sim[, c("age", "biomarker", "smoking")],
  treatment_level = 0,
  analysis = "transport",
  estimator = "dr",
  se_method = "none"
)

print(auc_result)
#> 
#> Transportable AUC Estimation
#> --------------------------------------------- 
#> Analysis: transport 
#> Estimator: dr 
#> Treatment level: 0 
#> N target: 1224  | N source: 1276 
#> 
#> Estimate: 0.6204
#> 
#> Naive estimate: 0.6249
```

### Calibration in the Target Population

To assess calibration in the target population:

``` r
calib_result <- tr_calibration(
  predictions = transport_sim$risk_score,
  outcomes = transport_sim$event,
  treatment = transport_sim$treatment,
  source = transport_sim$source,
  covariates = transport_sim[, c("age", "biomarker", "smoking")],
  treatment_level = 0,
  analysis = "transport",
  estimator = "ipw",
  smoother = "loess",
  se_method = "none"
)

print(calib_result)
#> 
#> Transportable CALIBRATION Estimation
#> --------------------------------------------- 
#> Analysis: transport 
#> Estimator: ipw 
#> Treatment level: 0 
#> N target: 1224  | N source: 1276 
#> 
#> Calibration Metrics:
#>   ICI (Integrated Calibration Index): 0.0546 
#>   E50 (Median absolute error): 0.0558 
#>   E90 (90th percentile error): 0.0852 
#>   Emax (Maximum error): 0.1806
```

``` r
plot(calib_result)
#> Install 'patchwork' package for histogram subplot.
```

![Transportable calibration curve for the target
population](transportability_files/figure-html/plot-calibration-1.png)

Transportable calibration curve for the target population

## Joint Analysis

The **joint analysis** pools source and target data for potentially more
efficient estimation. This requires outcome data in both populations.

``` r
mse_joint <- tr_mse(
  predictions = transport_sim$risk_score,
  outcomes = transport_sim$event,
  treatment = transport_sim$treatment,
  source = transport_sim$source,
  covariates = transport_sim[, c("age", "biomarker", "smoking")],
  treatment_level = 0,
  analysis = "joint",  # Pool data
  estimator = "dr",
  se_method = "none"
)

cat("Transport MSE:", round(mse_result$estimate, 4), "\n")
#> Transport MSE: 0.2211
cat("Joint MSE:", round(mse_joint$estimate, 4), "\n")
#> Joint MSE: 0.2222
```

## Bootstrap Standard Errors

For inference, use bootstrap standard errors:

``` r
mse_with_se <- tr_mse(
  predictions = transport_sim$risk_score,
  outcomes = transport_sim$event,
  treatment = transport_sim$treatment,
  source = transport_sim$source,
  covariates = transport_sim[, c("age", "biomarker", "smoking")],
  treatment_level = 0,
  analysis = "transport",
  estimator = "dr",
  se_method = "bootstrap",
  n_boot = 500,
  stratified_boot = TRUE  # Preserve source/target ratio
)

summary(mse_with_se)
#> 
#> Summary: Transportable MSE Estimation
#> ======================================================= 
#> 
#> Call:
#> tr_mse(predictions = transport_sim$risk_score, outcomes = transport_sim$event, 
#>     treatment = transport_sim$treatment, source = transport_sim$source, 
#>     covariates = transport_sim[, c("age", "biomarker", "smoking")], 
#>     treatment_level = 0, analysis = "transport", estimator = "dr", 
#>     se_method = "bootstrap", n_boot = 500, stratified_boot = TRUE)
#> 
#> Settings:
#>   Analysis type: transport 
#>   Estimator: dr 
#>   Treatment level: 0 
#>   Target sample size: 1224 
#>   Source sample size: 1276 
#> 
#> Results:
#>      Estimator Estimate       SE CI_lower CI_upper
#>  Transportable   0.2211 0.008868   0.2043    0.239
#>          Naive   0.2178       NA       NA       NA
#> 
#> Difference (Transportable - Naive): 0.0033
confint(mse_with_se)
#>             2.5%     97.5%
#> tr_mse 0.2042813 0.2390305
```

The `stratified_boot = TRUE` option (default) ensures that bootstrap
samples preserve the ratio of source to target observations, which is
recommended for transportability analysis.

## Key Assumptions

The transportability estimators rely on several assumptions:

### 1. Consistency in the Source and Target Populations.

For all individuals $i$, we have $Y_{i}^{a} = Y_{i}$ if $A_{i} = a$.

The observed outcome equals the potential outcome under the received
treatment. Implies no interference and well-defined treatments.

### 2. Conditional Exchangeability in the Source Population (Trial)

$$\left( Y^{a}\bot A \right)|X,S = 1$$

Treatment is randomized in the source population, so there is no
confounding between treatment and outcome given covariates.

### 3. Positivity of Treatment in the Source Population (Trial)

$$P\left( A = a|X,S = 1 \right) > 0{\mspace{6mu}\text{for the treatment level of interest}}$$

In the source population, the treatment level must have positive
probability (guaranteed by randomization).

### 4. Conditional Exchangeability (Transportability)

$$\left( Y^{a}\bot S \right)|X$$

The potential outcome is independent of population membership given
covariates. This means the source population is “representative” of how
the outcome model would perform in the target, after conditioning on X.

### 5. Positivity of Selection

$$P\left( S = 0|X \right) > 0{\mspace{6mu}\text{and}\mspace{6mu}}P\left( S = 1|X \right) > 0{\mspace{6mu}\text{for all}\mspace{6mu}}X$$

Both populations must have positive probability for all covariate
values.

## When to Use Each Function

| Scenario                                        | Function                                                                                                                                                                                                                                        | Analysis      |
|-------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------|
| Single population, observational data           | [`cf_mse()`](https://boyercb.github.io/cfperformance/reference/cf_mse.md), [`cf_auc()`](https://boyercb.github.io/cfperformance/reference/cf_auc.md), [`cf_calibration()`](https://boyercb.github.io/cfperformance/reference/cf_calibration.md) | \-            |
| Source (RCT) to target, outcomes in source only | [`tr_mse()`](https://boyercb.github.io/cfperformance/reference/tr_mse.md), [`tr_auc()`](https://boyercb.github.io/cfperformance/reference/tr_auc.md), [`tr_calibration()`](https://boyercb.github.io/cfperformance/reference/tr_calibration.md) | `"transport"` |
| Source + target, outcomes in both               | [`tr_mse()`](https://boyercb.github.io/cfperformance/reference/tr_mse.md), [`tr_auc()`](https://boyercb.github.io/cfperformance/reference/tr_auc.md), [`tr_calibration()`](https://boyercb.github.io/cfperformance/reference/tr_calibration.md) | `"joint"`     |

## Comparison with cf\_\* Functions

The `cf_*` functions (counterfactual) are for single-population
analysis:

- Use observational data from one population
- Adjust for confounding between treatment and outcome
- Equivalent to “Observational (OBS)” estimators in Voter et al.

The `tr_*` functions (transportability) are for two-population analysis:

- Transport from source to target population
- Leverage randomization in source (if RCT)
- Adjust for differences in covariate distributions

## References

Voter, S. R., et al. (2025). “Transportability of machine learning-based
counterfactual prediction models with application to CASS.” *Diagnostic
and Prognostic Research*, 9(4).
[doi:10.1186/s41512-025-00201-y](https://doi.org/10.1186/s41512-025-00201-y)

Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
“Estimating and evaluating counterfactual prediction models.”
*Statistics in Medicine*, 44(23-24), e70287.
[doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

Dahabreh, I. J., Robertson, S. E., Tchetgen, E. J., Stuart, E. A., &
Hernán, M. A. (2019). “Generalizing causal inferences from randomized
trials: counterfactual and graphical identification.” *Biometrics*,
75(2), 685-694.
