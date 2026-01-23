# Plot Method for tr_calibration Objects

Creates a calibration plot showing predicted vs observed probabilities
in the target population.

## Usage

``` r
# S3 method for class 'tr_calibration'
plot(x, add_reference = TRUE, show_metrics = TRUE, ...)
```

## Arguments

- x:

  A `tr_calibration` object.

- add_reference:

  Logical; add 45-degree reference line (default: TRUE).

- show_metrics:

  Logical; show calibration metrics on plot (default: TRUE).

- ...:

  Additional arguments passed to plotting functions.

## Value

A ggplot object (if ggplot2 available) or base R plot.
