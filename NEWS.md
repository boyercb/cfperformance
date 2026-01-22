# cfperformance 0.1.0

Initial release implementing counterfactual prediction model performance 
estimation from Boyer, Dahabreh & Steingrimsson (2025).

## Features

### Performance Metrics
* `cf_mse()` - Counterfactual MSE/Brier score estimation
* `cf_auc()` - Counterfactual AUC estimation  
* `cf_calibration()` - Counterfactual calibration curves with ICI, E50, E90, Emax

### Estimators
* Naive estimator (subset-based)
* Conditional Loss (CL) / Outcome modeling estimator
* Inverse Probability Weighting (IPW) estimator
* Doubly Robust (DR) estimator

### Inference
* Bootstrap standard errors with parallel support
* Influence function-based analytic standard errors
* Cross-fitting (sample splitting) for DR estimation
* Confidence intervals via normal approximation or percentile bootstrap

### Model Selection
* `cf_cv()` - K-fold cross-validation with counterfactual metrics
* `cf_compare()` - Compare multiple prediction models

### Data
* `cvd_sim` - Simulated cardiovascular disease dataset for examples

## References

Boyer CB, Dahabreh IJ, Steingrimsson JA. Counterfactual prediction model 
performance. *Statistics in Medicine*. 2025; 44(23-24):e70287. 
doi:10.1002/sim.70287
