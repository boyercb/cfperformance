# Unified Prediction from Nuisance Models

Handles prediction from both standard glm/gam models and ml_fitted
objects.

## Usage

``` r
.predict_nuisance(model, newdata, type = "response")
```

## Arguments

- model:

  A fitted model (glm, gam, or ml_fitted).

- newdata:

  Data frame for prediction.

- type:

  Prediction type (default: "response").

## Value

Numeric vector of predictions.
