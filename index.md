# cfperformance

## Overview

`cfperformance` provides methods for estimating model performance
measures (MSE, AUC, calibration) under hypothetical/counterfactual
interventions. These methods are essential when:

- A prediction model will be deployed in settings where treatment
  policies differ from the training setting
- Predictions are meant to support decisions about treatment initiation
- You need valid performance estimates even when the prediction model is
  misspecified

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
- **Multiple Estimators:** Naive, Conditional Loss, IPW, and Doubly
  Robust
- **Inference:** Bootstrap and influence function-based standard errors
- **Cross-validation:** Counterfactual-aware model selection with
  [`cf_cv()`](https://boyercb.github.io/cfperformance/reference/cf_cv.md)
  and
  [`cf_compare()`](https://boyercb.github.io/cfperformance/reference/cf_compare.md)

## Example

``` r
library(cfperformance)
data(cvd_sim)

# Compare estimators
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

## Documentation

See
[`vignette("introduction", package = "cfperformance")`](https://boyercb.github.io/cfperformance/articles/introduction.md)
for a comprehensive introduction.

## Citation

If you use this package in your research, please cite:

Boyer CB, Dahabreh IJ, Steingrimsson JA. Estimating and evaluating
counterfactual prediction models. *Statistics in Medicine*. 2025;
44(23-24):e70287.
<doi:%5B10.1002/sim.70287>\](<https://doi.org/10.1002/sim.70287>)

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
```

## License

MIT License
