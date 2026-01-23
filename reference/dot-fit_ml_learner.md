# Fit an ML Learner to Data

Fit an ML Learner to Data

## Usage

``` r
.fit_ml_learner(learner, formula, data, family = "binomial")
```

## Arguments

- learner:

  An `ml_learner` object.

- formula:

  Model formula.

- data:

  Training data.

- family:

  Character: "binomial" for classification, "gaussian" for regression.

## Value

A fitted model object with class `ml_fitted`.
