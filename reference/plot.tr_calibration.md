# Plot Method for tr_calibration Objects

Creates a calibration plot showing predicted vs observed probabilities
in the target population, with optional confidence bands and histogram.

## Usage

``` r
# S3 method for class 'tr_calibration'
plot(
  x,
  add_reference = TRUE,
  show_metrics = TRUE,
  add_histogram = TRUE,
  add_rug = FALSE,
  add_ci = !is.null(x$boot_curves),
  ...
)
```

## Arguments

- x:

  A `tr_calibration` object.

- add_reference:

  Logical; add 45-degree reference line (default: TRUE).

- show_metrics:

  Logical; show calibration metrics on plot (default: TRUE).

- add_histogram:

  Logical; add histogram of predictions below the calibration curve
  (default: TRUE).

- add_rug:

  Logical; add rug plot to show individual predictions (default: FALSE).

- add_ci:

  Logical; add bootstrap confidence bands if available (default: TRUE if
  bootstrap was run).

- ...:

  Additional arguments passed to plotting functions.

## Value

A ggplot object (if ggplot2 available) or base R plot.
