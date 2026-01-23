# cfperformance

## Overview

`cfperformance` provides methods for estimating model performance
measures (MSE, AUC, calibration) under hypothetical/counterfactual
interventions. These methods are essential when:

- A prediction model will be deployed in settings where treatment
  policies differ from the training setting
- Predictions are meant to support decisions about treatment initiation
- You want to assess model performance after transporting from a source
  (e.g., RCT) to a target population

Based on Boyer, Dahabreh & Steingrimsson (2025). “Estimating and
evaluating counterfactual prediction models.” *Statistics in Medicine*,
44(23-24), e70287.
[doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

## Installation

``` r
# Install development version from GitHub
# install.packages("devtools")
devtools::install_github("boyercb/cfperformance")
```

## Quick Start

``` r
library(cfperformance)

# Load example data
data(cvd_sim)

# Estimate counterfactual MSE under no treatment
result <- cf_mse(
  predictions = cvd_sim$risk_score,
  outcomes = cvd_sim$event,
  treatment = cvd_sim$treatment,
  covariates = cvd_sim[, c("age", "bp", "chol")],
  treatment_level = 0,     # Evaluate under no treatment
  estimator = "dr"         # Doubly robust estimator
)

result
```

## Key Features

- **MSE/Brier Score:** Loss-based performance under counterfactual
  intervention
- **AUC:** Discrimination ability under counterfactual intervention  
- **Calibration:** Reliability of risk predictions under counterfactual
  intervention
- **Multiple Estimators:** Naive, Conditional Loss/Outcome Model, IPW,
  and Doubly Robust
- **Inference:** Bootstrap and influence function-based standard errors
- **Cross-validation:** Counterfactual-aware model selection with
  [`cf_cv()`](https://boyercb.github.io/cfperformance/reference/cf_cv.md)
  and
  [`cf_compare()`](https://boyercb.github.io/cfperformance/reference/cf_compare.md)
- **Transportability:** Evaluate model performance when transporting
  from a source population (e.g., RCT) to a target population

## Counterfactual Performance Estimation

Estimate how a prediction model would perform under a hypothetical
treatment policy (e.g., if everyone received or avoided treatment):

``` r
library(cfperformance)
data(cvd_sim)

# Compare estimators for counterfactual MSE
estimators <- c("naive", "cl", "ipw", "dr")
sapply(estimators, function(est) {
  cf_mse(
    predictions = cvd_sim$risk_score,
    outcomes = cvd_sim$event,
    treatment = cvd_sim$treatment,
    covariates = cvd_sim[, c("age", "bp", "chol")],
    treatment_level = 0,
    estimator = est
  )$estimate
})

# Estimate counterfactual AUC
cf_auc(
  predictions = cvd_sim$risk_score,
  outcomes = cvd_sim$event,
  treatment = cvd_sim$treatment,
  covariates = cvd_sim[, c("age", "bp", "chol")],
  treatment_level = 0,
  estimator = "dr"
)

# Cross-validation for model selection
cf_compare(
  models = list(
    "Simple" = event ~ age,
    "Full" = event ~ age + bp + chol
  ),
  data = cvd_sim,
  treatment = "treatment",
  treatment_level = 0,
  metric = "mse",
  K = 5
)
```

## Transportability Analysis

The package also implements transportability estimators from
Steingrimsson et al. (2022) and Voter et al. (2025) for evaluating
prediction model performance when transporting from a source population
(typically an RCT) to a target population:

``` r
# Load transportability example data
data(transport_sim)

# Estimate MSE in target population
tr_mse(
  predictions = transport_sim$risk_score,
  outcomes = transport_sim$event,
  treatment = transport_sim$treatment,
  source = transport_sim$source,  # 1=source/RCT, 0=target

  covariates = transport_sim[, c("age", "biomarker", "smoking")],
  treatment_level = 0,
  analysis = "transport",
  estimator = "dr"
)

# Estimate AUC in the target population
tr_auc(
  predictions = transport_sim$risk_score,
  outcomes = transport_sim$event,
  treatment = transport_sim$treatment,
  source = transport_sim$source,
  covariates = transport_sim[, c("age", "biomarker", "smoking")],
  treatment_level = 0,
  analysis = "transport",
  estimator = "dr"
)

# Estimate calibration in the target population
tr_calibration(
  predictions = transport_sim$risk_score,
  outcomes = transport_sim$event,
  treatment = transport_sim$treatment,
  source = transport_sim$source,
  covariates = transport_sim[, c("age", "biomarker", "smoking")],
  treatment_level = 0,
  analysis = "transport",
  estimator = "ipw"
)
```

See
[`vignette("transportability", package = "cfperformance")`](https://boyercb.github.io/cfperformance/articles/transportability.md)
for details.

## Machine Learning Integration

The package supports flexible ML methods for nuisance model estimation
with automatic cross-fitting for valid inference:

``` r
# Use random forest for propensity scores and outcome models
cf_mse(
  predictions = cvd_sim$risk_score,
  outcomes = cvd_sim$event,
  treatment = cvd_sim$treatment,
  covariates = cvd_sim[, c("age", "bp", "chol")],
  treatment_level = 0,
  estimator = "dr",
  propensity_model = ml_learner("ranger", num.trees = 500),
  outcome_model = ml_learner("xgboost", nrounds = 100),
  cross_fit = TRUE
)
```

**Supported learners:** - `ranger` - Fast random forests - `xgboost` -
Gradient boosting - `grf` - Generalized random forests (honest
estimation) - `glmnet` - Elastic net with cross-validated λ -
`superlearner` - Ensemble learning - `custom` - User-supplied
fit/predict functions

See
[`vignette("ml-integration", package = "cfperformance")`](https://boyercb.github.io/cfperformance/articles/ml-integration.md)
for details.

## Documentation

See
[`vignette("introduction", package = "cfperformance")`](https://boyercb.github.io/cfperformance/articles/introduction.md)
for a comprehensive introduction.

## Citation

If you use this package in your research, please cite:

Boyer CB, Dahabreh IJ, Steingrimsson JA. Estimating and evaluating
counterfactual prediction models. *Statistics in Medicine*. 2025;
44(23-24):e70287. doi:
[10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

For transportability methods, also cite:

Steingrimsson JA, Gatsonis C, Li B, Dahabreh IJ. Transporting a
Prediction Model for Use in a New Target Population. *American Journal
of Epidemiology*. 2022; 192(2):296-304. doi:
[10.1093/aje/kwac128](https://doi.org/10.1093/aje/kwac128)

Voter SR, et al. Transportability of machine learning-based
counterfactual prediction models with application to CASS. *Diagnostic
and Prognostic Research*. 2025; 9(4). doi:
[10.1186/s41512-025-00201-y](https://doi.org/10.1186/s41512-025-00201-y)

``` bibtex
@article{boyer2025estimating,
  title={Estimating and Evaluating Counterfactual Prediction Models},
  author={Boyer, Christopher B. and Dahabreh, Issa J. and Steingrimsson, Jon A.},
  journal={Statistics in Medicine},
  volume={44},
  number={23-24},
  pages={e70287},
  year={2025},
  doi={10.1002/sim.70287}
}

@article{10.1093/aje/kwac128,
    title = {Transporting a Prediction Model for Use in a New Target Population},
    author = {Steingrimsson, Jon A. and Gatsonis, Constantine and Li, Bing and Dahabreh, Issa J.},
    journal = {American Journal of Epidemiology},
    volume = {192},
    number = {2},
    pages = {296-304},
    year = {2022},
    doi = {10.1093/aje/kwac128}
}

@article{voter2025transportability,
  title={Transportability of machine learning-based counterfactual prediction models with application to CASS},
  author = {Voter, Sarah C. and Dahabreh, Issa J. and Boyer, Christopher B. and Rahbar, Habib and Kontos, Despina and Steingrimsson, Jon A.},
  journal={Diagnostic and Prognostic Research},
  volume={9},
  number={4},
  year={2025},
  doi={10.1186/s41512-025-00201-y}
}
```

## License

MIT License
