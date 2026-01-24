# Influence Function-Based Standard Errors for MSE

Computes standard errors using the efficient influence function for the
doubly robust MSE estimator.

## Usage

``` r
.influence_se_mse(
  predictions,
  outcomes,
  treatment,
  covariates,
  treatment_level,
  estimator,
  propensity_model,
  outcome_model,
  ps_trim_spec = NULL
)
```

## Arguments

- predictions:

  Numeric vector of model predictions.

- outcomes:

  Numeric vector of observed outcomes.

- treatment:

  Numeric vector of treatment indicators.

- covariates:

  Matrix or data frame of covariates.

- treatment_level:

  Counterfactual treatment level.

- estimator:

  Character string specifying estimator type.

- propensity_model:

  Fitted propensity model.

- outcome_model:

  Fitted outcome model.

## Value

Standard error estimate.

## Details

The efficient influence function for the MSE under treatment level \\a\\
is: \$\$\chi(O) = \frac{I(A=a)}{e_a(X)}\[L(Y, \mu(X^\*)) - h_a(X)\] +
h_a(X) - \psi(a)\$\$

where \\e_a(X) = P(A=a\|X)\\ is the propensity score, \\h_a(X) =
E\[L\|X,A=a\]\\ is the conditional loss, and \\\psi(a)\\ is the target
estimand.

## References

Boyer, C. B., Dahabreh, I. J., & Steingrimsson, J. A. (2025).
"Estimating and evaluating counterfactual prediction models."
*Statistics in Medicine*, 44(23-24), e70287.
[doi:10.1002/sim.70287](https://doi.org/10.1002/sim.70287)
