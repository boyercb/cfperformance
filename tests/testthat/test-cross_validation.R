# Tests for cross-validation functions

test_that("cf_cv runs with default settings", {
  set.seed(123)
  n <- 200

  data <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  data$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * data$x1))
  data$y <- rbinom(n, 1, plogis(-1 + data$x1 + 0.5 * data$x2 - 0.3 * data$a))

  result <- cf_cv(
    formula = y ~ x1 + x2,
    data = data,
    treatment = "a",
    treatment_level = 0,
    metric = "mse",
    K = 3
  )

  expect_s3_class(result, "cf_cv")
  expect_equal(nrow(result$results), 3)  # 3 folds
  expect_true("mse" %in% names(result$results))
  expect_true("mse_mean" %in% names(result$summary))
  expect_true(result$summary$mse_mean > 0)
})

test_that("cf_cv works with repeated CV", {
  set.seed(123)
  n <- 200

  data <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  data$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * data$x1))
  data$y <- rbinom(n, 1, plogis(-1 + data$x1 + 0.5 * data$x2 - 0.3 * data$a))

  result <- cf_cv(
    formula = y ~ x1 + x2,
    data = data,
    treatment = "a",
    K = 3,
    repeats = 2
  )

  expect_equal(nrow(result$results), 6)  # 3 folds x 2 repeats
  expect_equal(result$repeats, 2)
})

test_that("cf_cv computes AUC metric", {
  set.seed(123)
  n <- 200

  data <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  data$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * data$x1))
  data$y <- rbinom(n, 1, plogis(-1 + data$x1 + 0.5 * data$x2 - 0.3 * data$a))

  result <- cf_cv(
    formula = y ~ x1 + x2,
    data = data,
    treatment = "a",
    metric = "auc",
    K = 3
  )

  expect_true("auc" %in% names(result$results))
  expect_true("auc_mean" %in% names(result$summary))
  expect_true(result$summary$auc_mean > 0.5)  # Better than random
  expect_true(result$summary$auc_mean <= 1)
})

test_that("cf_cv computes both MSE and AUC", {
  set.seed(123)
  n <- 200

  data <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  data$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * data$x1))
  data$y <- rbinom(n, 1, plogis(-1 + data$x1 + 0.5 * data$x2 - 0.3 * data$a))

  result <- cf_cv(
    formula = y ~ x1 + x2,
    data = data,
    treatment = "a",
    metric = "both",
    K = 3
  )

  expect_true("mse" %in% names(result$results))
  expect_true("auc" %in% names(result$results))
  expect_true("mse_mean" %in% names(result$summary))
  expect_true("auc_mean" %in% names(result$summary))
})

test_that("cf_cv works with different estimators", {
  set.seed(123)
  n <- 200

  data <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  data$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * data$x1))
  data$y <- rbinom(n, 1, plogis(-1 + data$x1 + 0.5 * data$x2 - 0.3 * data$a))

  estimators <- c("naive", "cl", "ipw", "dr")

  for (est in estimators) {
    result <- cf_cv(
      formula = y ~ x1 + x2,
      data = data,
      treatment = "a",
      estimator = est,
      K = 3
    )
    expect_equal(result$estimator, est)
    expect_true(!is.na(result$summary$mse_mean))
  }
})

test_that("cf_compare works with CV method", {
  set.seed(123)
  n <- 200

  data <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )
  data$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * data$x1))
  data$y <- rbinom(n, 1, plogis(-1 + data$x1 + 0.5 * data$x2 - 0.3 * data$a))

  models <- list(
    "Simple" = y ~ x1,
    "Full" = y ~ x1 + x2 + x3
  )

  result <- cf_compare(
    models = models,
    data = data,
    treatment = "a",
    metric = "mse",
    method = "cv",
    K = 3
  )

  expect_s3_class(result, "cf_compare")
  expect_equal(nrow(result$results), 2)
  expect_true(result$best_model %in% c("Simple", "Full"))
  expect_true("mse_mean" %in% names(result$results))
})

test_that("cf_compare works with holdout method", {
  set.seed(123)
  n <- 200

  data <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n)
  )
  data$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * data$x1))
  data$y <- rbinom(n, 1, plogis(-1 + data$x1 + 0.5 * data$x2 - 0.3 * data$a))

  models <- list(
    "Simple" = y ~ x1,
    "Full" = y ~ x1 + x2 + x3
  )

  result <- cf_compare(
    models = models,
    data = data,
    treatment = "a",
    metric = "mse",
    method = "holdout",
    test_prop = 0.3
  )

  expect_s3_class(result, "cf_compare")
  expect_equal(result$method, "holdout")
  expect_true("mse_se" %in% names(result$results))
})

test_that("cf_compare selects best model correctly", {
  set.seed(42)
  n <- 300

  data <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    noise = rnorm(n)  # Pure noise
  )
  data$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * data$x1))
  data$y <- rbinom(n, 1, plogis(-1 + 2 * data$x1 + data$x2 - 0.3 * data$a))

  models <- list(
    "Good" = y ~ x1 + x2,
    "Noise" = y ~ noise
  )

  result <- cf_compare(
    models = models,
    data = data,
    treatment = "a",
    metric = "mse",
    K = 5,
    seed = 42
  )

  # Good model should have lower MSE
  good_mse <- result$results$mse_mean[result$results$model == "Good"]
  noise_mse <- result$results$mse_mean[result$results$model == "Noise"]

  expect_true(good_mse < noise_mse)
  expect_equal(result$best_model, "Good")
})

test_that("cf_compare selects by AUC correctly", {
  set.seed(42)
  n <- 300

  data <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    noise = rnorm(n)
  )
  data$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * data$x1))
  data$y <- rbinom(n, 1, plogis(-1 + 2 * data$x1 + data$x2 - 0.3 * data$a))

  models <- list(
    "Good" = y ~ x1 + x2,
    "Noise" = y ~ noise
  )

  result <- cf_compare(
    models = models,
    data = data,
    treatment = "a",
    metric = "auc",
    K = 5,
    seed = 42
  )

  # Good model should have higher AUC
  good_auc <- result$results$auc_mean[result$results$model == "Good"]
  noise_auc <- result$results$auc_mean[result$results$model == "Noise"]

  expect_true(good_auc > noise_auc)
  expect_equal(result$best_model, "Good")
})

test_that("print methods work for cf_cv", {
  set.seed(123)
  n <- 100

  data <- data.frame(x = rnorm(n))
  data$a <- rbinom(n, 1, 0.5)
  data$y <- rbinom(n, 1, plogis(data$x))

  result <- cf_cv(
    formula = y ~ x,
    data = data,
    treatment = "a",
    K = 3
  )

  # Should not error
  expect_output(print(result), "Counterfactual Cross-Validation")
  expect_output(summary(result), "Summary")
})

test_that("print methods work for cf_compare", {
  set.seed(123)
  n <- 100

  data <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  data$a <- rbinom(n, 1, 0.5)
  data$y <- rbinom(n, 1, plogis(data$x1))

  result <- cf_compare(
    models = list("M1" = y ~ x1, "M2" = y ~ x1 + x2),
    data = data,
    treatment = "a",
    K = 3
  )

  expect_output(print(result), "Counterfactual Model Comparison")
  expect_output(summary(result), "Summary")
})

test_that("cf_cv validates inputs", {
  data <- data.frame(x = 1:10, y = 1:10, a = rep(0:1, 5))

  # Missing outcome variable

  expect_error(
    cf_cv(z ~ x, data = data, treatment = "a"),
    "not found in data"
  )

  # Missing treatment variable
  expect_error(
    cf_cv(y ~ x, data = data, treatment = "missing"),
    "not found in data"
  )

  # Non-binary outcome for AUC
  expect_error(
    cf_cv(y ~ x, data = data, treatment = "a", metric = "auc"),
    "AUC requires binary outcome"
  )
})

test_that("stratified folds preserve class balance", {
  set.seed(123)
  n <- 100
  y <- c(rep(1, 20), rep(0, 80))  # 20% positive

  folds <- cfperformance:::.stratified_folds(y, K = 5)

  # Check each fold has similar proportion
  for (k in 1:5) {
    fold_y <- y[folds == k]
    prop <- mean(fold_y)
    expect_true(prop > 0.1 && prop < 0.3,
                info = paste("Fold", k, "proportion:", prop))
  }
})
