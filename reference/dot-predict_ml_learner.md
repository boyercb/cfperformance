# Predict from Fitted ML Learner

Predict from Fitted ML Learner

## Usage

``` r
.predict_ml_learner(object, newdata, type = "response", ...)
```

## Arguments

- object:

  An `ml_fitted` object.

- newdata:

  Data frame for prediction.

- type:

  Character: "response" for probabilities (default), "class" for class
  predictions.

- ...:

  Ignored.

## Value

Numeric vector of predicted probabilities.
