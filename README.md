# cfperformance: Counterfactual Prediction Model Performance

## Overview

`cfperformance` provides methods for estimating model performance measures (MSE, AUC, calibration) under hypothetical/counterfactual interventions. These methods are essential when:

- A prediction model will be deployed in settings where treatment policies differ from the training setting
- Predictions are meant to support decisions about treatment initiation
- You need valid performance estimates even when the prediction model is misspecified

Based on Boyer, Dahabreh & Steingrimsson (2025). "Estimating and evaluating counterfactual prediction models." *Statistics in Medicine*, 44(23-24), e70287.

## Installation

```r
# Install from CRAN (once available)
install.packages("cfperformance")

# Install development version from GitHub
# devtools::install_github("boyercb/cfperformance")
```
## Quick Start

```r
library(cfperformance)

# Estimate counterfactual MSE
result <- cf_mse(
  predictions = model_predictions,
  outcomes = observed_outcomes,
  treatment = treatment_indicator,
  covariates = covariate_matrix,
  treatment_level = 0,     # Evaluate under no treatment
  estimator = "dr"         # Doubly robust estimator
)

print(result)
# Counterfactual MSE: 0.152 (SE: 0.021)
# 95% CI: [0.111, 0.193]
# Estimator: Doubly Robust

# Compare with naive estimate
summary(result)
```

## Key Features
- **MSE/Brier Score:** Loss-based performance under counterfactual intervention
- **AUC:** Discrimination ability under counterfactual intervention  
- **Calibration:** Reliability of risk predictions under counterfactual intervention
- **Multiple Estimators:** Conditional Loss, IPW, and Doubly Robust
- **Inference:** Bootstrap and influence function-based standard errors

## Documentation

- [Introduction to Counterfactual Performance](vignettes/introduction.html)
- [MSE Estimation](vignettes/mse-estimation.html)
- [AUC Estimation](vignettes/auc-estimation.html)
- [Calibration Curves](vignettes/calibration.html)

## Citation

```bibtex
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
