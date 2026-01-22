# cfperformance: Counterfactual Prediction Model Performance Estimation

Provides methods for estimating model performance measures (mean squared
error, area under the ROC curve, calibration) under
hypothetical/counterfactual interventions. Implements conditional loss,
inverse probability weighting, and doubly robust estimators from Boyer,
Dahabreh & Steingrimsson (2025)
[doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287) . These
methods are essential when prediction models will be deployed in
settings where treatment policies differ from training, or when
predictions support treatment decisions.

Provides methods for estimating model performance measures (MSE, AUC,
calibration) under hypothetical/counterfactual interventions.

## Main Functions

- [`cf_mse`](https://boyercb.github.io/cfperformance/reference/cf_mse.md):
  Estimate counterfactual mean squared error

- [`cf_auc`](https://boyercb.github.io/cfperformance/reference/cf_auc.md):
  Estimate counterfactual AUC

- [`cf_calibration`](https://boyercb.github.io/cfperformance/reference/cf_calibration.md):
  Estimate counterfactual calibration curve

- [`fit_nuisance`](https://boyercb.github.io/cfperformance/reference/fit_nuisance.md):
  Fit nuisance models for estimation

## Estimators

Three main estimators are implemented:

- **Conditional Loss (CL)**: Relies on correct specification of the
  outcome model \\h_a(X) = E\[L(Y, \mu(X^\*)) \| X, A=a\]\\

- **Inverse Probability Weighting (IPW)**: Relies on correct
  specification of the propensity score \\e_a(X) = Pr\[A=a \| X\]\\

- **Doubly Robust (DR)**: Consistent if either the outcome model or
  propensity score is correctly specified

## References

Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
"Estimating and evaluating counterfactual prediction models." Statistics
in Medicine. [doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)

## See also

Useful links:

- <https://github.com/boyercb/cfperformance>

- Report bugs at <https://github.com/boyercb/cfperformance/issues>

## Author

**Maintainer**: Christopher Boyer <cboyer@hsph.harvard.edu>
([ORCID](https://orcid.org/0000-0003-0935-8722))

Authors:

- Issa Dahabreh ([ORCID](https://orcid.org/0000-0002-1195-1850))

- Jon Steingrimsson ([ORCID](https://orcid.org/0000-0002-8284-8877))
