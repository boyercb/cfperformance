# Script to create example dataset for cfperformance package
# This simulates data similar to the motivating example in Boyer et al. (2025)

set.seed(20250122)

n <- 2500

# Covariates: age (scaled), blood pressure, cholesterol
age <- rnorm(n, mean = 0, sd = 1)
bp <- rnorm(n, mean = 0, sd = 1)
chol <- rnorm(n, mean = 0, sd = 1)

# Treatment assignment (confounded by covariates)
# Higher age and blood pressure -> more likely to receive treatment
propensity <- plogis(-0.3 + 0.4 * age + 0.3 * bp + 0.1 * chol)
treatment <- rbinom(n, 1, propensity)

# Outcome (cardiovascular event)
# Risk depends on age, bp, chol and treatment reduces risk
linear_pred <- -2 + 0.6 * age + 0.5 * bp + 0.3 * chol - 0.4 * treatment
prob_event <- plogis(linear_pred)
event <- rbinom(n, 1, prob_event)

# Create the dataset
cvd_sim <- data.frame(
  age = age,
  bp = bp,
  chol = chol,
  treatment = treatment,
  event = event
)

# Add prediction model scores (from a model fit on observed data)
# This mimics a scenario where we have a prediction model
pred_model <- glm(event ~ age + bp + chol, data = cvd_sim, family = binomial)
cvd_sim$risk_score <- predict(pred_model, type = "response")

# Save the dataset
usethis::use_data(cvd_sim, overwrite = TRUE)
