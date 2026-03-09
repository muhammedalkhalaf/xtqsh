# xtqsh: Panel Quantile Regression Slope Homogeneity Test

An R package for testing slope homogeneity across panels in quantile regression
models. Implements the methodology of Galvao, Juhl, Montes-Rojas, and Olmo (2017).

## Installation

```r
# Install from source
install.packages("xtqsh", repos = NULL, type = "source")
```

## Overview

The `xtqsh` package tests whether slope coefficients are equal across panels
at various quantiles of the conditional distribution. This is useful for
determining whether pooled quantile regression is appropriate or whether
panel-specific slopes are needed.

### Null Hypothesis

$$H_0: \beta_1(\tau) = \beta_2(\tau) = \cdots = \beta_n(\tau)$$

### Alternative Hypothesis

$$H_1: \beta_i(\tau) \neq \beta_j(\tau) \text{ for some } i \neq j$$

## Test Statistics

Two test statistics are provided:

- **S-hat statistic**: Chi-squared asymptotics for fixed n, large T
- **D-hat statistic**: Standard normal asymptotics for large n and T

## Usage

```r
library(xtqsh)

# Load example data
data(qsh_sample)

# Test slope homogeneity at multiple quantiles
fit <- xtqsh(
  formula = y ~ x1 + x2,
  data = qsh_sample,
  id = "id",
  time = "time",
  tau = c(0.25, 0.50, 0.75)
)

# View results
summary(fit)

# Include marginal tests for each variable
fit_marg <- xtqsh(
  formula = y ~ x1 + x2,
  data = qsh_sample,
  id = "id",
  time = "time",
  tau = c(0.25, 0.50, 0.75),
  marginal = TRUE
)
summary(fit_marg)
```

## Minimum Distance Estimator

Under slope homogeneity, the package computes the optimal minimum distance
estimator that combines individual panel estimates:

$$\hat{\beta}_{MD}(\tau) = \left(\sum_{i=1}^{n} \hat{V}_i^{-1}\right)^{-1} \sum_{i=1}^{n} \hat{V}_i^{-1} \hat{\beta}_i(\tau)$$

## Bandwidth Selection

Two methods are available for estimating the sparsity function:

- `hallsheather`: Hall and Sheather (1988) bandwidth (default)
- `bofinger`: Bofinger (1975) bandwidth

## References

Galvao AF, Juhl T, Montes-Rojas G, Olmo J (2017). "Testing Slope Homogeneity
in Quantile Regression Panel Data with an Application to the Cross-Section
of Stock Returns." *Journal of Financial Econometrics*, 16(2), 211-243.
[doi:10.1093/jjfinec/nbx027](https://doi.org/10.1093/jjfinec/nbx027)

Koenker R, Bassett G (1978). "Regression Quantiles." *Econometrica*,
46(1), 33-50. [doi:10.2307/1913643](https://doi.org/10.2307/1913643)

## License

GPL-3
