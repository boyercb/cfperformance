# Tests for nuisance model fitting

test_that("fit_nuisance returns correct structure", {
  set.seed(42)
  n <- 200
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  df$a <- rbinom(n, 1, plogis(-0.5 + 0.5 * df$x1))
  df$y <- rbinom(n, 1, plogis(-1 + df$x1 + 0.5 * df$x2))

  result <- fit_nuisance(
    propensity_formula = a ~ x1 + x2,
    outcome_formula = y ~ x1 + x2,
    data = df,
    treatment_level = 0
  )

  expect_s3_class(result, "cf_nuisance")
  expect_true("propensity" %in% names(result))
  expect_true("outcome" %in% names(result))
  expect_s3_class(result$propensity, "glm")
  expect_s3_class(result$outcome, "glm")
})


test_that("fit_nuisance outcome model uses correct subset", {
  set.seed(42)
  n <- 200
  df <- data.frame(
    x = rnorm(n),
    a = rbinom(n, 1, 0.5),
    y = rbinom(n, 1, 0.3)
  )

  result <- fit_nuisance(
    propensity_formula = a ~ x,
    outcome_formula = y ~ x,
    data = df,
    treatment_level = 0
  )

  # Outcome model should be fit only on A=0 observations
  n_untreated <- sum(df$a == 0)
  expect_equal(nobs(result$outcome), n_untreated)
})


test_that("fit_nuisance print method works", {
  set.seed(42)
  n <- 100
  df <- data.frame(
    x = rnorm(n),
    a = rbinom(n, 1, 0.5),
    y = rbinom(n, 1, 0.3)
  )

  result <- fit_nuisance(
    propensity_formula = a ~ x,
    outcome_formula = y ~ x,
    data = df
  )

  expect_output(print(result), "Nuisance Models")
  expect_output(print(result), "Propensity")
  expect_output(print(result), "Outcome")
})
