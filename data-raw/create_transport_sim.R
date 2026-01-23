# Create transportability example dataset
set.seed(2025)
n <- 1500

# Covariates
age <- rnorm(n, mean = 60, sd = 10)
biomarker <- rnorm(n)
smoking <- rbinom(n, 1, 0.3)

# Source indicator: RCT patients are younger, less smoking
# S=1 for source/RCT, S=0 for target
source_prob <- plogis(-0.3 + 0.025 * (age - 60) + 0.5 * smoking)
source <- rbinom(n, 1, 1 - source_prob)

# Treatment assignment
# Randomized in source (RCT), confounded in target
treatment <- ifelse(
  source == 1,
  rbinom(n, 1, 0.5),
  rbinom(n, 1, plogis(-0.3 + 0.015 * age + 0.2 * biomarker))
)

# Outcome
true_risk <- plogis(-2.5 + 0.03 * age + 0.4 * biomarker + 0.3 * smoking - 0.5 * treatment)
event <- rbinom(n, 1, true_risk)

# Predictions (model trained on RCT, assumes A=0)
risk_score <- plogis(-2.5 + 0.025 * age + 0.35 * biomarker + 0.2 * smoking)

# Create data frame with standardized covariates
transport_sim <- data.frame(
  age = as.numeric(scale(age)),
  biomarker = biomarker,
  smoking = smoking,
  source = source,
  treatment = treatment,
  event = event,
  risk_score = risk_score
)

# Save
save(transport_sim, file = "data/transport_sim.rda", compress = "xz")
cat("Created transport_sim.rda\n")
cat("Source n =", sum(source == 1), "\n")
cat("Target n =", sum(source == 0), "\n")
