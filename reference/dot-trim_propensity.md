# Trim propensity scores to avoid extreme values

Internal function for trimming propensity scores using either absolute
bounds or quantile-based trimming.

## Usage

``` r
.trim_propensity(
  ps,
  method = "absolute",
  bounds = c(0.01, 0.99),
  weights = NULL
)
```

## Arguments

- ps:

  Numeric vector of propensity scores (probabilities).

- method:

  Character string specifying trimming method:

  - `"absolute"`: Trim to fixed probability bounds (default)

  - `"quantile"`: Trim based on quantiles of the propensity distribution

  - `"none"`: No trimming

- bounds:

  Numeric vector of length 2 specifying trimming bounds. For
  `method = "absolute"`: probability bounds (default: c(0.01, 0.99)).
  For `method = "quantile"`: quantile bounds (default: c(0.01, 0.99)).

- weights:

  Optional numeric vector of weights for quantile calculation. Only used
  when `method = "quantile"`.

## Value

Trimmed propensity scores.

## Details

For `method = "absolute"`, propensity scores are clipped to
`[bounds[1], bounds[2]]`. This is the most common approach and ensures
that no propensity score falls below or above the specified thresholds.

For `method = "quantile"`, the bounds are computed as the `bounds[1]`
and `bounds[2]` quantiles of the propensity score distribution, and then
scores are clipped to these data-dependent bounds. This can be useful
when the propensity score distribution varies across datasets.
