# Generate sample panel data for xtqsh package
# Data is generated under slope homogeneity (common beta across panels)

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

# Save data
usethis::use_data(qsh_sample, overwrite = TRUE)

cat("Generated qsh_sample with", nrow(qsh_sample), "observations\n")
cat("True coefficients: beta1 =", beta1, ", beta2 =", beta2, "\n")
