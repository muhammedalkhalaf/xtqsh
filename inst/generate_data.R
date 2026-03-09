# Script to generate sample data for xtqsh package
# Run this from the package root directory

set.seed(12345)

# Panel dimensions
n_panels <- 20   # Number of cross-sections
T_periods <- 100 # Number of time periods

# True parameters (homogeneous slopes)
beta1 <- 1.5
beta2 <- -0.8

# Generate panel data
qsh_sample <- data.frame(
  id = rep(1:n_panels, each = T_periods),
  time = rep(1:T_periods, times = n_panels)
)

# Panel-specific fixed effects
alpha <- rnorm(n_panels, mean = 0, sd = 2)

# Generate independent variables
qsh_sample$x1 <- rnorm(nrow(qsh_sample), mean = 0, sd = 1)
qsh_sample$x2 <- rnorm(nrow(qsh_sample), mean = 0, sd = 1)

# Generate errors (heteroskedastic across panels for realism)
sigma_i <- runif(n_panels, min = 0.5, max = 1.5)
errors <- rnorm(nrow(qsh_sample), mean = 0, sd = sigma_i[qsh_sample$id])

# Generate dependent variable
qsh_sample$y <- alpha[qsh_sample$id] + beta1 * qsh_sample$x1 +
                beta2 * qsh_sample$x2 + errors

# Save to data directory
save(qsh_sample, file = "data/qsh_sample.rda", compress = "xz")

cat("Generated qsh_sample with", nrow(qsh_sample), "observations\n")
cat("True coefficients: beta1 =", beta1, ", beta2 =", beta2, "\n")
