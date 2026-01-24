# Parse ps_trim argument into method and bounds

Internal function to parse the `ps_trim` argument into method and
bounds.

## Usage

``` r
.parse_ps_trim(ps_trim)
```

## Arguments

- ps_trim:

  Either:

  - A list with `method` and `bounds` elements

  - A numeric vector of length 2 (interpreted as absolute bounds)

  - A single number (interpreted as symmetric absolute bounds: c(x,
    1-x))

  - "none" or NULL for no trimming

  - "quantile" for quantile-based trimming with default bounds

## Value

A list with `method` and `bounds` elements.
