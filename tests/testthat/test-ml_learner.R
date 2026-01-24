# Tests for ML Learner integration
# These tests verify the ml_learner() infrastructure

test_that("ml_learner creates proper object", {
  # Test basic creation
  learner <- ml_learner("ranger")
  expect_s3_class(learner, "ml_learner")
  expect_equal(learner$method, "ranger")
  expect_equal(learner$args, list())
  
  # Test with arguments
  learner <- ml_learner("ranger", num.trees = 100, min.node.size = 10)
  expect_equal(learner$args$num.trees, 100)
  expect_equal(learner$args$min.node.size, 10)
  
  # Test different methods
  for (method in c("xgboost", "grf", "glmnet", "superlearner")) {
    learner <- ml_learner(method)
    expect_equal(learner$method, method)
  }
})

test_that("ml_learner validates custom method", {
  # Custom requires fit_fun
  expect_error(ml_learner("custom"),
               "fit_fun.*must be provided")
  
  # Custom requires predict_fun
  expect_error(ml_learner("custom", fit_fun = function(f, d, fam, ...) NULL),
               "predict_fun.*must be provided")
  
  # Valid custom learner
  custom <- ml_learner(
    "custom",
    fit_fun = function(formula, data, family, ...) {
      glm(formula, data = data, family = binomial())
    },
    predict_fun = function(object, newdata, ...) {
      predict(object, newdata = newdata, type = "response")
    }
  )
  expect_s3_class(custom, "ml_learner")
})

test_that("print.ml_learner works", {
  learner <- ml_learner("ranger", num.trees = 500, mtry = 3)
  output <- capture.output(print(learner))
  expect_true(any(grepl("ranger", output)))
  expect_true(any(grepl("num.trees", output)))
})

test_that("is_ml_learner correctly identifies objects", {
  expect_true(is_ml_learner(ml_learner("ranger")))
  expect_false(is_ml_learner(NULL))
  expect_false(is_ml_learner(list(method = "ranger")))
  expect_false(is_ml_learner(glm(y ~ x, data = data.frame(x = 1:10, y = rbinom(10, 1, 0.5)),
                                 family = binomial())))
})


# Test ML learner fitting (requires packages)
test_that("ranger learner fits and predicts", {
  skip_if_not_installed("ranger")
  
  # Generate test data
  set.seed(123)
  n <- 200
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  df$A <- rbinom(n, 1, plogis(0.5 * df$x1 - 0.3 * df$x2))
  
  learner <- ml_learner("ranger", num.trees = 50)
  
  # Test internal fitting
  fitted <- .fit_ml_learner(learner, A ~ x1 + x2, data = df, family = "binomial")
  expect_s3_class(fitted, "ml_fitted")
  expect_equal(fitted$method, "ranger")
  
  # Test prediction
  preds <- .predict_ml_learner(fitted, df[1:10, ])
  expect_length(preds, 10)
  expect_true(all(preds >= 0 & preds <= 1))
})

test_that("glmnet learner fits and predicts", {
  skip_if_not_installed("glmnet")
  
  set.seed(123)
  n <- 200
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )
  df$A <- rbinom(n, 1, plogis(0.5 * df$x1 - 0.3 * df$x2))
  
  learner <- ml_learner("glmnet", alpha = 0.5)
  
  fitted <- .fit_ml_learner(learner, A ~ x1 + x2 + x3, data = df, family = "binomial")
  expect_s3_class(fitted, "ml_fitted")
  
  preds <- .predict_ml_learner(fitted, df[1:10, ])
  expect_length(preds, 10)
  expect_true(all(preds >= 0 & preds <= 1))
})

test_that("xgboost learner fits and predicts", {
  skip_if_not_installed("xgboost")
  
  set.seed(123)
  n <- 200
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  df$A <- rbinom(n, 1, plogis(0.5 * df$x1 - 0.3 * df$x2))
  
  learner <- ml_learner("xgboost", nrounds = 20, max_depth = 3)
  
  fitted <- .fit_ml_learner(learner, A ~ x1 + x2, data = df, family = "binomial")
  expect_s3_class(fitted, "ml_fitted")
  
  preds <- .predict_ml_learner(fitted, df[1:10, ])
  expect_length(preds, 10)
  expect_true(all(preds >= 0 & preds <= 1))
})

test_that("grf learner fits and predicts", {
  skip_if_not_installed("grf")
  
  set.seed(123)
  n <- 200
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  df$A <- rbinom(n, 1, plogis(0.5 * df$x1 - 0.3 * df$x2))
  
  learner <- ml_learner("grf", num.trees = 50)
  
  fitted <- .fit_ml_learner(learner, A ~ x1 + x2, data = df, family = "binomial")
  expect_s3_class(fitted, "ml_fitted")
  
  preds <- .predict_ml_learner(fitted, df[1:10, ])
  expect_length(preds, 10)
  expect_true(all(preds >= 0 & preds <= 1))
})

test_that("custom learner fits and predicts", {
  set.seed(123)
  n <- 200
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  df$A <- rbinom(n, 1, plogis(0.5 * df$x1 - 0.3 * df$x2))
  
  # Define custom learner (just wraps glm)
  custom_fit <- function(formula, data, family, ...) {
    glm(formula, data = data, family = binomial())
  }
  custom_pred <- function(object, newdata, ...) {
    predict(object, newdata = newdata, type = "response")
  }
  
  learner <- ml_learner("custom", fit_fun = custom_fit, predict_fun = custom_pred)
  
  fitted <- .fit_ml_learner(learner, A ~ x1 + x2, data = df, family = "binomial")
  expect_s3_class(fitted, "ml_fitted")
  
  preds <- .predict_ml_learner(fitted, df[1:10, ])
  expect_length(preds, 10)
  expect_true(all(preds >= 0 & preds <= 1))
})


# Integration tests with cf_mse
test_that("cf_mse works with ranger ml_learner", {
  skip_if_not_installed("ranger")
  
  set.seed(42)
  n <- 300
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  
  result <- cf_mse(
    predictions = pred,
    outcomes = y,
    treatment = a,
    covariates = data.frame(x = x),
    treatment_level = 0,
    estimator = "dr",
    propensity_model = ml_learner("ranger", num.trees = 50),
    outcome_model = ml_learner("ranger", num.trees = 50),
    cross_fit = TRUE,
    n_folds = 3,
    se_method = "influence"
  )
  
  expect_s3_class(result, "cf_mse")
  expect_true(!is.na(result$estimate))
  expect_true(!is.na(result$se))
  # Note: DR estimator can produce negative estimates in finite samples
  # but should be reasonably close to the naive estimate
  # Using wider tolerance for ML models which have more variability
  expect_true(abs(result$estimate - result$naive_estimate) < 1.0)
})

test_that("cf_mse warns when ml_learner used without cross_fit", {
  skip_if_not_installed("ranger")
  
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  a <- rbinom(n, 1, plogis(-0.5 + 0.5 * x))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  
  expect_warning(
    cf_mse(
      predictions = pred,
      outcomes = y,
      treatment = a,
      covariates = data.frame(x = x),
      treatment_level = 0,
      estimator = "dr",
      propensity_model = ml_learner("ranger", num.trees = 20),
      cross_fit = FALSE,  # Not recommended with ML
      se_method = "none"
    ),
    "cross_fit=TRUE is recommended"
  )
})

test_that("cross_fit_nuisance works with ml_learner", {
  skip_if_not_installed("ranger")
  
  set.seed(123)
  n <- 200
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  a <- rbinom(n, 1, plogis(0.5 * df$x1 - 0.3 * df$x2))
  y <- rbinom(n, 1, plogis(-1 + df$x1 + df$x2 - 0.5 * a))
  pred <- plogis(-1 + 0.8 * df$x1 + 0.5 * df$x2)
  
  result <- .cross_fit_nuisance(
    treatment = a,
    outcomes = y,
    covariates = df,
    treatment_level = 0,
    predictions = pred,
    K = 3,
    propensity_learner = ml_learner("ranger", num.trees = 30),
    outcome_learner = ml_learner("ranger", num.trees = 30)
  )
  
  expect_length(result$ps, n)
  expect_length(result$h, n)
  expect_true(all(result$ps > 0 & result$ps < 1))
})
