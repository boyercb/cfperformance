# Plot Method for cf_calibration Objects

Creates a calibration plot showing predicted vs observed probabilities
under the counterfactual intervention.

## Usage

``` r
# S3 method for class 'cf_calibration'
plot(x, add_histogram = TRUE, add_rug = TRUE, ...)
```

## Arguments

- x:

  A `cf_calibration` object.

- add_histogram:

  Logical; add histogram of predictions (default: TRUE).

- add_rug:

  Logical; add rug plot (default: TRUE).

- ...:

  Additional arguments passed to plotting functions.

## Value

A ggplot object (if ggplot2 available) or base R plot.
