# DeLong's Method for AUC Standard Error

Computes the standard error of the empirical AUC using DeLong's method.

## Usage

``` r
.delong_se(predictions, outcomes)
```

## Arguments

- predictions:

  Numeric vector of predictions.

- outcomes:

  Binary outcome vector.

## Value

Standard error of AUC.

## References

DeLong, E. R., DeLong, D. M., & Clarke-Pearson, D. L. (1988). "Comparing
the areas under two or more correlated receiver operating characteristic
curves: a nonparametric approach." *Biometrics*, 44(3), 837-845.
