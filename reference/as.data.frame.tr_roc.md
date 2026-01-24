# Convert ROC Curve to Data Frame

Converts an ROC curve object to a data frame suitable for use with
ggplot2.

## Usage

``` r
# S3 method for class 'tr_roc'
as.data.frame(x, ...)

# S3 method for class 'cf_roc'
as.data.frame(x, ...)
```

## Arguments

- x:

  An object of class `tr_roc` or `cf_roc`.

- ...:

  Additional arguments (ignored).

## Value

A data frame with columns:

- threshold:

  Classification threshold

- fpr:

  False positive rate

- sensitivity:

  Sensitivity (TPR)

- specificity:

  Specificity

- type:

  Either "adjusted" or "naive"

## Examples

``` r
set.seed(123)
n <- 500
x <- rnorm(n)
a <- rbinom(n, 1, 0.5)
y <- rbinom(n, 1, plogis(-1 + x))
pred <- plogis(-1 + 0.8 * x)

roc <- cf_roc(
  predictions = pred,
  outcomes = y,
  treatment = a,
  covariates = data.frame(x = x),
  n_thresholds = 21
)
df <- as.data.frame(roc)
head(df)
#>   threshold       fpr sensitivity specificity     type
#> 1      0.00 1.0000000   1.0000000 0.000000000 adjusted
#> 2      0.05 0.9918223   0.9998694 0.008177702 adjusted
#> 3      0.10 0.9239015   0.9981700 0.076098496 adjusted
#> 4      0.15 0.7785828   0.9759714 0.221417225 adjusted
#> 5      0.20 0.6416111   0.8795424 0.358388853 adjusted
#> 6      0.25 0.4724056   0.7960865 0.527594393 adjusted
```
