# Tests for tr_auc() function

test_that("tr_auc works with basic inputs", {
  skip_on_cran()

  set.seed(123)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.5),
              rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "none"
  )

  expect_s3_class(result, "tr_auc")
  expect_s3_class(result, "tr_performance")
  expect_true(result$estimate >= 0 && result$estimate <= 1)
  expect_equal(result$analysis, "transport")
  expect_equal(result$estimator, "dr")
})

test_that("tr_auc works with joint analysis", {
  skip_on_cran()

  set.seed(456)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.5),
              rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "joint",
    estimator = "dr",
    se_method = "none"
  )

  expect_s3_class(result, "tr_auc")
  expect_equal(result$analysis, "joint")
  expect_true(result$estimate >= 0 && result$estimate <= 1)
})

test_that("tr_auc naive estimator works", {
  skip_on_cran()

  set.seed(789)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.5),
              rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "naive",
    se_method = "none"
  )

  expect_s3_class(result, "tr_auc")
  expect_equal(result$estimator, "naive")
  expect_true(result$estimate >= 0 && result$estimate <= 1)
})

test_that("tr_auc om estimator works for transport", {
  skip_on_cran()

  set.seed(111)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.5),
              rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "om",
    se_method = "none"
  )

  expect_s3_class(result, "tr_auc")
  expect_equal(result$estimator, "om")
  expect_true(result$estimate >= 0 && result$estimate <= 1)
})

test_that("tr_auc ipw estimator works for transport", {
  skip_on_cran()

  set.seed(222)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.5),
              rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "ipw",
    se_method = "none"
  )

  expect_s3_class(result, "tr_auc")
  expect_equal(result$estimator, "ipw")
  expect_true(result$estimate >= 0 && result$estimate <= 1)
})

test_that("tr_auc om estimator works for joint", {
  skip_on_cran()

  set.seed(333)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.5),
              rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "joint",
    estimator = "om",
    se_method = "none"
  )

  expect_s3_class(result, "tr_auc")
  expect_equal(result$analysis, "joint")
  expect_equal(result$estimator, "om")
  expect_true(result$estimate >= 0 && result$estimate <= 1)
})

test_that("tr_auc ipw estimator works for joint", {
  skip_on_cran()

  set.seed(444)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.5),
              rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "joint",
    estimator = "ipw",
    se_method = "none"
  )

  expect_s3_class(result, "tr_auc")
  expect_equal(result$analysis, "joint")
  expect_equal(result$estimator, "ipw")
  expect_true(result$estimate >= 0 && result$estimate <= 1)
})

test_that("tr_auc bootstrap SE works", {
  skip_on_cran()

  set.seed(555)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.5),
              rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "bootstrap",
    n_boot = 50  # Small for testing
  )

  expect_s3_class(result, "tr_auc")
  expect_true(!is.null(result$se))
  expect_true(!is.na(result$se))
  expect_true(result$se > 0)
  expect_true(!is.null(result$ci_lower))
  expect_true(!is.null(result$ci_upper))
  expect_true(result$ci_lower < result$ci_upper)
})

test_that("tr_auc works with treatment_level = 1", {
  skip_on_cran()

  set.seed(666)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.5),
              rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 1,
    analysis = "transport",
    estimator = "dr",
    se_method = "none"
  )

  expect_s3_class(result, "tr_auc")
  expect_equal(result$treatment_level, 1)
  expect_true(result$estimate >= 0 && result$estimate <= 1)
})

test_that("tr_auc requires binary outcomes", {
  skip_on_cran()

  set.seed(777)
  n <- 100
  x <- rnorm(n)
  s <- rbinom(n, 1, 0.5)
  a <- rbinom(n, 1, 0.5)
  y <- rnorm(n)  # Continuous outcomes
  pred <- plogis(x)
  covars <- data.frame(x = x)

  expect_error(
    tr_auc(
      predictions = pred,
      outcomes = y,
      treatment = a,
      source = s,
      covariates = covars,
      treatment_level = 0,
      estimator = "dr",
      se_method = "none"
    ),
    "binary outcomes"
  )
})

test_that("tr_auc validates inputs", {
  skip_on_cran()

  set.seed(888)
  n <- 100
  x <- rnorm(n)
  s <- rbinom(n, 1, 0.5)
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, 0.5)
  pred <- plogis(x)
  covars <- data.frame(x = x)

  # Mismatched lengths
  expect_error(
    tr_auc(
      predictions = pred[1:50],
      outcomes = y,
      treatment = a,
      source = s,
      covariates = covars,
      estimator = "dr",
      se_method = "none"
    )
  )

  # Invalid source values
  expect_error(
    tr_auc(
      predictions = pred,
      outcomes = y,
      treatment = a,
      source = rep(2, n),  # Invalid
      covariates = covars,
      estimator = "dr",
      se_method = "none"
    )
  )
})

test_that("tr_auc print method works", {
  skip_on_cran()

  set.seed(999)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "none"
  )

  expect_output(print(result), "Transportable AUC")
  expect_output(print(result), "Estimate:")
})

test_that("tr_auc summary method works", {
  skip_on_cran()

  set.seed(1010)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "bootstrap",
    n_boot = 30
  )

  expect_output(summary(result), "Transportable AUC")
  expect_output(summary(result), "Source sample size")
  expect_output(summary(result), "Target sample size")
})

test_that("tr_auc coef method works", {
  skip_on_cran()

  set.seed(1111)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "none"
  )

  coefs <- coef(result)
  expect_true(is.numeric(coefs))
  expect_true(coefs >= 0 && coefs <= 1)
})

test_that("tr_auc returns correct structure", {
  skip_on_cran()

  set.seed(1212)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "none"
  )

  # Check structure
  expect_true("estimate" %in% names(result))
  expect_true("estimator" %in% names(result))
  expect_true("analysis" %in% names(result))
  expect_true("metric" %in% names(result))
  expect_true("treatment_level" %in% names(result))
  expect_true("n_target" %in% names(result))
  expect_true("n_source" %in% names(result))
  expect_true("naive_estimate" %in% names(result))

  expect_equal(result$metric, "auc")
  expect_true(result$n_target > 0)
  expect_true(result$n_source > 0)
})

test_that("tr_auc naive estimate is computed correctly", {
  skip_on_cran()

  set.seed(1313)
  n <- 500
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "none"
  )

  # Naive estimate should be stored and valid
  expect_true(!is.null(result$naive_estimate))
  expect_true(result$naive_estimate >= 0 && result$naive_estimate <= 1)
})

test_that("tr_auc with multiple covariates", {
  skip_on_cran()

  set.seed(1414)
  n <- 500
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.2 * x1 + 0.1 * x2))
  a <- ifelse(s == 1,
              rbinom(n, 1, 0.5),
              rbinom(n, 1, plogis(-0.5 + 0.3 * x1 + 0.2 * x2)))
  y <- rbinom(n, 1, plogis(-1 + 0.5 * x1 + 0.3 * x2 - 0.5 * a))
  pred <- plogis(-1 + 0.4 * x1 + 0.2 * x2)
  covars <- data.frame(x1 = x1, x2 = x2)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "none"
  )

  expect_s3_class(result, "tr_auc")
  expect_true(result$estimate >= 0 && result$estimate <= 1)
})

test_that("tr_auc confint method works", {
  skip_on_cran()

  set.seed(1515)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "bootstrap",
    n_boot = 30
  )

  ci <- confint(result)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 1)
  expect_equal(ncol(ci), 2)
  expect_true(ci[1, 1] < ci[1, 2])
})


# =============================================================================
# Tests for stratified bootstrap
# =============================================================================

test_that("tr_auc stratified bootstrap works", {
  skip_on_cran()

  set.seed(1616)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "bootstrap",
    n_boot = 30,
    stratified_boot = TRUE
  )

  expect_s3_class(result, "tr_auc")
  expect_true(!is.na(result$se))
  expect_true(result$se > 0)
})

test_that("tr_auc non-stratified bootstrap works", {
  skip_on_cran()

  set.seed(1717)
  n <- 300
  x <- rnorm(n)
  s <- rbinom(n, 1, plogis(0.5 - 0.3 * x))
  a <- ifelse(s == 1, rbinom(n, 1, 0.5), rbinom(n, 1, plogis(-0.5 + 0.5 * x)))
  y <- rbinom(n, 1, plogis(-1 + x - 0.5 * a))
  pred <- plogis(-1 + 0.8 * x)
  covars <- data.frame(x = x)

  result <- tr_auc(
    predictions = pred,
    outcomes = y,
    treatment = a,
    source = s,
    covariates = covars,
    treatment_level = 0,
    analysis = "transport",
    estimator = "dr",
    se_method = "bootstrap",
    n_boot = 30,
    stratified_boot = FALSE
  )

  expect_s3_class(result, "tr_auc")
  expect_true(!is.na(result$se))
  expect_true(result$se > 0)
})
