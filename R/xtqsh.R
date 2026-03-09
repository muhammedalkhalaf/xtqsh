#' Panel Quantile Regression Slope Homogeneity Test
#'
#' Tests for slope homogeneity across panels in quantile regression models.
#' Implements the S-hat statistic (chi-squared asymptotics) and D-hat statistic
#' (standard normal asymptotics) from Galvao, Juhl, Montes-Rojas, and Olmo (2017).
#'
#' @param formula A formula specifying the model (e.g., \code{y ~ x1 + x2}).
#' @param data A data frame containing panel data with variables specified in
#'   the formula.
#' @param id Character string specifying the panel (cross-section) identifier
#'   variable name.
#' @param time Character string specifying the time variable name.
#' @param tau Numeric vector of quantiles to test, each in (0,1).
#'   Default is \code{c(0.25, 0.50, 0.75)}.
#' @param bw Character string specifying the bandwidth selection method:
#'   \code{"hallsheather"} (default) or \code{"bofinger"}.
#' @param hac Integer specifying the HAC lag truncation parameter. If 0 (default),
#'   no HAC correction is applied.
#' @param marginal Logical. If \code{TRUE}, also compute marginal (per-variable)
#'   slope homogeneity tests. Default is \code{FALSE}.
#'
#' @return An object of class \code{"xtqsh"} containing:
#' \describe{
#'   \item{S}{Matrix of S-hat statistics (joint test) by quantile}
#'   \item{D}{Matrix of D-hat statistics (joint test) by quantile}
#'   \item{pval_S}{Matrix of p-values for S-hat statistics}
#'   \item{pval_D}{Matrix of p-values for D-hat statistics}
#'   \item{S_ols}{S-hat statistic for OLS (mean regression)}
#'   \item{D_ols}{D-hat statistic for OLS}
#'   \item{pval_S_ols}{P-value for S_ols}
#'   \item{pval_D_ols}{P-value for D_ols}
#'   \item{beta_md}{Matrix of minimum distance QR estimates by quantile}
#'   \item{beta_md_se}{Matrix of standard errors for beta_md}
#'   \item{beta_all}{Array of individual panel QR estimates}
#'   \item{beta_ols_all}{Matrix of individual panel OLS estimates}
#'   \item{S_marginal}{Matrix of marginal S-hat statistics (if \code{marginal = TRUE})}
#'   \item{D_marginal}{Matrix of marginal D-hat statistics (if \code{marginal = TRUE})}
#'   \item{pval_S_marginal}{P-values for marginal S-hat (if \code{marginal = TRUE})}
#'   \item{pval_D_marginal}{P-values for marginal D-hat (if \code{marginal = TRUE})}
#'   \item{tau}{Vector of quantiles tested}
#'   \item{n_panels}{Number of panels}
#'   \item{valid_panels}{Number of panels with valid estimates}
#'   \item{k}{Number of slope coefficients}
#'   \item{n_obs}{Total number of observations}
#'   \item{depvar}{Dependent variable name}
#'   \item{indepvars}{Independent variable names}
#'   \item{bw}{Bandwidth method used}
#'   \item{hac}{HAC lag truncation used}
#'   \item{call}{The matched call}
#' }
#'
#' @details
#' The null hypothesis is slope homogeneity across all panels:
#' \deqn{H_0: \beta_1(\tau) = \beta_2(\tau) = \cdots = \beta_n(\tau)}
#'
#' The alternative is that at least one panel has different slopes:
#' \deqn{H_1: \beta_i(\tau) \neq \beta_j(\tau) \text{ for some } i \neq j}
#'
#' Two test statistics are computed:
#'
#' \strong{S-hat statistic:} Under the null hypothesis with T tending to infinity
#' and n fixed, the S-hat statistic follows a chi-squared distribution with
#' k(n-1) degrees of freedom, where k is the number of slope coefficients and
#' n is the number of panels.
#'
#' \strong{D-hat statistic:} Under the null hypothesis with both n and T tending

#' to infinity, the D-hat statistic is asymptotically standard normal. This
#' statistic is standardized to have mean 0 and variance 1 under the null.
#'
#' The minimum distance (MD) estimator under homogeneity is computed as:
#' \deqn{\hat{\beta}_{MD}(\tau) = \left(\sum_{i=1}^{n} \hat{V}_i^{-1}\right)^{-1}
#' \sum_{i=1}^{n} \hat{V}_i^{-1} \hat{\beta}_i(\tau)}
#'
#' where \eqn{\hat{V}_i} is the estimated variance-covariance matrix of the
#' individual panel quantile regression estimator \eqn{\hat{\beta}_i(\tau)}.
#'
#' @references
#' Galvao AF, Juhl T, Montes-Rojas G, Olmo J (2017). "Testing Slope Homogeneity
#' in Quantile Regression Panel Data with an Application to the Cross-Section
#' of Stock Returns." \emph{Journal of Financial Econometrics}, 16(2), 211-243.
#' \doi{10.1093/jjfinec/nbx027}
#'
#' Koenker R, Bassett G (1978). "Regression Quantiles." \emph{Econometrica},
#' 46(1), 33-50. \doi{10.2307/1913643}
#'
#' Hall P, Sheather SJ (1988). "On the Distribution of a Studentized Quantile."
#' \emph{Journal of the Royal Statistical Society: Series B}, 50(3), 381-391.
#' \doi{10.1111/j.2517-6161.1988.tb01735.x}
#'
#' Bofinger E (1975). "Estimation of a Density Function Using Order Statistics."
#' \emph{Australian Journal of Statistics}, 17(1), 1-7.
#' \doi{10.1111/j.1467-842X.1975.tb01366.x}
#'
#' @examples
#' \donttest{
#' # Load example panel data
#' data(qsh_sample)
#'
#' # Test slope homogeneity at 25th, 50th, and 75th quantiles
#' fit <- xtqsh(
#'   formula = y ~ x1 + x2,
#'   data = qsh_sample,
#'   id = "id",
#'   time = "time",
#'   tau = c(0.25, 0.50, 0.75)
#' )
#'
#' # View results
#' summary(fit)
#'
#' # Include marginal tests
#' fit_marg <- xtqsh(
#'   formula = y ~ x1 + x2,
#'   data = qsh_sample,
#'   id = "id",
#'   time = "time",
#'   tau = c(0.25, 0.50, 0.75),
#'   marginal = TRUE
#' )
#' summary(fit_marg)
#' }
#'
#' @importFrom stats model.frame model.matrix model.response var cov na.omit
#'   pnorm pchisq dnorm qnorm formula terms complete.cases as.formula coef
#'   lm.fit residuals fitted sd
#' @importFrom quantreg rq
#' @export
xtqsh <- function(formula, data, id, time,
                  tau = c(0.25, 0.50, 0.75),
                  bw = c("hallsheather", "bofinger"),
                  hac = 0L,
                  marginal = FALSE) {

  # Match arguments
  bw <- match.arg(bw)
  call <- match.call()

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!id %in% names(data)) {
    stop("Panel identifier '", id, "' not found in data")
  }
  if (!time %in% names(data)) {
    stop("Time variable '", time, "' not found in data")
  }
  if (!all(tau > 0 & tau < 1)) {
    stop("All quantiles in 'tau' must be between 0 and 1 (exclusive)")
  }
  if (length(tau) < 1) {
    stop("At least one quantile required in 'tau'")
  }

  # Sort data by panel and time
  data <- data[order(data[[id]], data[[time]]), ]

  # Parse formula
  mf <- model.frame(formula, data = data, na.action = na.omit)
  y <- model.response(mf)
  X <- model.matrix(formula, data = mf)

  # Remove intercept from slope coefficient matrix
  has_intercept <- "(Intercept)" %in% colnames(X)
  if (has_intercept) {
    X_slopes <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  } else {
    X_slopes <- X
  }

  depvar <- all.vars(formula)[1]
  indepvars <- colnames(X_slopes)
  k <- ncol(X_slopes)

  if (k < 1) {
    stop("At least one independent variable required")
  }

  # Get panel information
  # Need to match rows in model frame to original data
  mf_rows <- as.numeric(rownames(mf))
  panel_ids <- data[[id]][mf_rows]
  time_vals <- data[[time]][mf_rows]

  panels <- unique(panel_ids)
  n_panels <- length(panels)
  ntau <- length(tau)

  # Initialize storage
  # Per-panel estimates: list of results for each panel
  panel_results <- vector("list", n_panels)
  names(panel_results) <- as.character(panels)

  valid_panels <- 0
  n_obs <- length(y)

  # Panel-by-panel estimation
  for (i in seq_along(panels)) {
    panel_id <- panels[i]
    idx <- which(panel_ids == panel_id)

    yi <- y[idx]
    Xi <- X[idx, , drop = FALSE]
    Xi_slopes <- X_slopes[idx, , drop = FALSE]
    Ti <- length(yi)

    # Need sufficient observations
    if (Ti < k + 5) {
      panel_results[[i]] <- list(valid = FALSE)
      next
    }

    # OLS estimation for this panel
    ols_fit <- tryCatch({
      lm.fit(Xi, yi)
    }, error = function(e) NULL)

    if (is.null(ols_fit)) {
      panel_results[[i]] <- list(valid = FALSE)
      next
    }

    ols_coef <- ols_fit$coefficients
    ols_resid <- ols_fit$residuals

    # OLS slope coefficients
    if (has_intercept) {
      beta_ols <- ols_coef[names(ols_coef) != "(Intercept)"]
    } else {
      beta_ols <- ols_coef
    }
    names(beta_ols) <- indepvars

    # OLS variance-covariance (robust to heteroskedasticity)
    V_ols <- .compute_panel_vcov_ols(Xi_slopes, ols_resid, hac)

    # QR estimation for each quantile
    beta_qr <- matrix(NA, nrow = ntau, ncol = k)
    V_qr <- array(NA, dim = c(k, k, ntau))
    colnames(beta_qr) <- indepvars

    for (ti in seq_along(tau)) {
      tauval <- tau[ti]

      qr_result <- tryCatch({
        fit <- quantreg::rq(yi ~ Xi - 1, tau = tauval)
        coefs <- coef(fit)

        # Extract slope coefficients
        if (has_intercept) {
          beta <- coefs[names(coefs) != "Xi(Intercept)"]
        } else {
          beta <- coefs
        }
        names(beta) <- indepvars

        # Compute variance-covariance using sparsity/bandwidth
        V <- .compute_panel_vcov_qr(Xi_slopes, yi, beta, tauval, bw, hac, Ti)

        list(beta = beta, V = V, success = TRUE)
      }, error = function(e) {
        list(beta = NULL, V = NULL, success = FALSE)
      })

      if (qr_result$success) {
        beta_qr[ti, ] <- qr_result$beta
        V_qr[, , ti] <- qr_result$V
      }
    }

    panel_results[[i]] <- list(
      valid = TRUE,
      Ti = Ti,
      beta_ols = beta_ols,
      V_ols = V_ols,
      beta_qr = beta_qr,
      V_qr = V_qr
    )

    valid_panels <- valid_panels + 1
  }

  if (valid_panels < 2) {
    stop("At least 2 valid panels required for homogeneity test")
  }

  # Compute test statistics
  # First, gather all valid panel results
  valid_idx <- which(sapply(panel_results, function(x) x$valid))
  n_valid <- length(valid_idx)

  # Degrees of freedom for chi-squared: k * (n - 1)
  df_S <- k * (n_valid - 1)

  # Initialize result matrices
  S_mat <- matrix(NA, nrow = 1, ncol = ntau)
  D_mat <- matrix(NA, nrow = 1, ncol = ntau)
  pS_mat <- matrix(NA, nrow = 1, ncol = ntau)
  pD_mat <- matrix(NA, nrow = 1, ncol = ntau)
  colnames(S_mat) <- colnames(D_mat) <- colnames(pS_mat) <- colnames(pD_mat) <- paste0("tau_", tau)

  beta_md <- matrix(NA, nrow = ntau, ncol = k)
  beta_md_se <- matrix(NA, nrow = ntau, ncol = k)
  colnames(beta_md) <- colnames(beta_md_se) <- indepvars
  rownames(beta_md) <- rownames(beta_md_se) <- paste0("tau_", tau)

  # Store all individual estimates
  beta_all <- array(NA, dim = c(n_valid, k, ntau))
  dimnames(beta_all) <- list(
    panel = as.character(panels[valid_idx]),
    var = indepvars,
    tau = paste0("tau_", tau)
  )

  beta_ols_all <- matrix(NA, nrow = n_valid, ncol = k)
  colnames(beta_ols_all) <- indepvars
  rownames(beta_ols_all) <- as.character(panels[valid_idx])

  # Fill individual estimates
  for (j in seq_along(valid_idx)) {
    pr <- panel_results[[valid_idx[j]]]
    beta_ols_all[j, ] <- pr$beta_ols
    for (ti in seq_along(tau)) {
      beta_all[j, , ti] <- pr$beta_qr[ti, ]
    }
  }

  # Compute joint test statistics for each quantile
  for (ti in seq_along(tau)) {
    result <- .compute_joint_test(panel_results, valid_idx, ti, k)

    if (!is.null(result)) {
      S_mat[1, ti] <- result$S
      D_mat[1, ti] <- result$D
      pS_mat[1, ti] <- pchisq(result$S, df = df_S, lower.tail = FALSE)
      pD_mat[1, ti] <- 2 * pnorm(abs(result$D), lower.tail = FALSE)
      beta_md[ti, ] <- result$beta_md
      beta_md_se[ti, ] <- result$beta_md_se
    }
  }

  # OLS test
  ols_result <- .compute_ols_test(panel_results, valid_idx, k)

  S_ols <- ols_result$S
  D_ols <- ols_result$D
  pS_ols <- pchisq(S_ols, df = df_S, lower.tail = FALSE)
  pD_ols <- 2 * pnorm(abs(D_ols), lower.tail = FALSE)

  # Marginal tests (if requested)
  S_marginal <- D_marginal <- pS_marginal <- pD_marginal <- NULL

  if (marginal) {
    df_S_marg <- n_valid - 1  # k = 1 for marginal tests

    S_marginal <- matrix(NA, nrow = k, ncol = ntau)
    D_marginal <- matrix(NA, nrow = k, ncol = ntau)
    pS_marginal <- matrix(NA, nrow = k, ncol = ntau)
    pD_marginal <- matrix(NA, nrow = k, ncol = ntau)

    rownames(S_marginal) <- rownames(D_marginal) <- indepvars
    rownames(pS_marginal) <- rownames(pD_marginal) <- indepvars
    colnames(S_marginal) <- colnames(D_marginal) <- paste0("tau_", tau)
    colnames(pS_marginal) <- colnames(pD_marginal) <- paste0("tau_", tau)

    for (j in seq_len(k)) {
      for (ti in seq_along(tau)) {
        marg_result <- .compute_marginal_test(panel_results, valid_idx, ti, j)

        if (!is.null(marg_result)) {
          S_marginal[j, ti] <- marg_result$S
          D_marginal[j, ti] <- marg_result$D
          pS_marginal[j, ti] <- pchisq(marg_result$S, df = df_S_marg, lower.tail = FALSE)
          pD_marginal[j, ti] <- 2 * pnorm(abs(marg_result$D), lower.tail = FALSE)
        }
      }
    }
  }

  # Construct result object
  result <- list(
    S = S_mat,
    D = D_mat,
    pval_S = pS_mat,
    pval_D = pD_mat,
    S_ols = S_ols,
    D_ols = D_ols,
    pval_S_ols = pS_ols,
    pval_D_ols = pD_ols,
    beta_md = beta_md,
    beta_md_se = beta_md_se,
    beta_all = beta_all,
    beta_ols_all = beta_ols_all,
    S_marginal = S_marginal,
    D_marginal = D_marginal,
    pval_S_marginal = pS_marginal,
    pval_D_marginal = pD_marginal,
    tau = tau,
    n_panels = n_panels,
    valid_panels = n_valid,
    k = k,
    df_S = df_S,
    n_obs = n_obs,
    depvar = depvar,
    indepvars = indepvars,
    bw = bw,
    hac = hac,
    marginal = marginal,
    call = call
  )

  class(result) <- "xtqsh"
  return(result)
}


#' @keywords internal
.compute_panel_vcov_ols <- function(X, resid, hac) {
  # Heteroskedasticity-robust variance-covariance matrix for OLS
  # V = (X'X)^{-1} X' diag(e^2) X (X'X)^{-1}

  n <- nrow(X)
  k <- ncol(X)

  XtX <- crossprod(X)
  XtX_inv <- tryCatch(solve(XtX), error = function(e) NULL)

  if (is.null(XtX_inv)) {
    return(diag(k) * NA)
  }

  if (hac > 0) {
    # HAC variance (Newey-West style)
    meat <- matrix(0, k, k)
    for (j in 0:hac) {
      weight <- if (j == 0) 1 else (1 - j / (hac + 1))
      for (t in (j + 1):n) {
        xt <- X[t, , drop = FALSE]
        xt_j <- X[t - j, , drop = FALSE]
        meat <- meat + weight * (resid[t] * resid[t - j]) * (t(xt) %*% xt_j + t(xt_j) %*% xt)
      }
    }
    # Remove double counting of j=0
    for (t in 1:n) {
      xt <- X[t, , drop = FALSE]
      meat <- meat - 0.5 * resid[t]^2 * (t(xt) %*% xt)
    }
    V <- XtX_inv %*% meat %*% XtX_inv
  } else {
    # HC0 variance
    e2 <- resid^2
    meat <- crossprod(X * sqrt(e2))
    V <- XtX_inv %*% meat %*% XtX_inv
  }

  return(V)
}


#' @keywords internal
.compute_panel_vcov_qr <- function(X, y, beta, tau, bw, hac, T_i) {
  # Variance-covariance matrix for quantile regression
  # V = tau(1-tau) * (f(0))^{-2} * (X'X)^{-1}
  # where f(0) is the density of errors at 0, estimated via bandwidth

  k <- ncol(X)

  XtX <- crossprod(X)
  XtX_inv <- tryCatch(solve(XtX), error = function(e) NULL)

  if (is.null(XtX_inv)) {
    return(diag(k) * NA)
  }

  # Compute residuals
  resid <- y - X %*% beta

  # Estimate sparsity (1/f(0)) using bandwidth
  h <- .compute_bandwidth(tau, T_i, bw)

  # Kernel density estimation at 0
  # Using indicator function (simple approach)
  f0 <- sum(abs(resid) <= h) / (2 * h * T_i)

  if (f0 < 1e-10) {
    # Fallback: use normal approximation
    sigma <- sd(resid)
    f0 <- dnorm(qnorm(tau)) / sigma
  }

  sparsity <- 1 / f0

  # Variance matrix
  V <- tau * (1 - tau) * sparsity^2 * XtX_inv / T_i

  if (hac > 0) {
    # HAC adjustment for serial correlation
    # Simplified Newey-West type correction
    psi <- ifelse(resid < 0, tau - 1, tau)

    meat <- matrix(0, k, k)
    for (j in 0:hac) {
      weight <- if (j == 0) 1 else (1 - j / (hac + 1))
      for (t in (j + 1):T_i) {
        xt <- X[t, , drop = FALSE]
        xt_j <- X[t - j, , drop = FALSE]
        meat <- meat + weight * psi[t] * psi[t - j] * (t(xt) %*% xt_j + t(xt_j) %*% xt)
      }
    }
    # Adjust for j=0 double counting
    for (t in 1:T_i) {
      xt <- X[t, , drop = FALSE]
      meat <- meat - 0.5 * psi[t]^2 * (t(xt) %*% xt)
    }

    V <- XtX_inv %*% meat %*% XtX_inv / T_i
  }

  return(V)
}


#' @keywords internal
.compute_bandwidth <- function(tau, n, method) {
  # Bandwidth selection for sparsity estimation
  # Following Hall-Sheather (1988) or Bofinger (1975)

  alpha <- 0.05  # nominal level
  z_alpha <- qnorm(1 - alpha / 2)

  if (method == "hallsheather") {
    # Hall-Sheather (1988) bandwidth
    f_tau <- dnorm(qnorm(tau))
    h <- n^(-1/3) * z_alpha^(2/3) * ((1.5 * f_tau^2) / (2 * qnorm(tau)^2 + 1))^(1/3)
  } else {
    # Bofinger (1975) bandwidth
    f_tau <- dnorm(qnorm(tau))
    h <- n^(-1/5) * ((4.5 * f_tau^4) / (2 * qnorm(tau)^2 + 1)^2)^(1/5)
  }

  # Ensure positive bandwidth
  h <- max(h, 1e-6)

  return(h)
}


#' @keywords internal
.compute_joint_test <- function(panel_results, valid_idx, tau_idx, k) {
  # Compute joint S and D statistics for a given quantile

  n <- length(valid_idx)

  # Collect beta estimates and variances
  betas <- matrix(NA, nrow = n, ncol = k)
  V_list <- vector("list", n)
  V_inv_list <- vector("list", n)

  for (j in seq_along(valid_idx)) {
    pr <- panel_results[[valid_idx[j]]]
    betas[j, ] <- pr$beta_qr[tau_idx, ]
    V_list[[j]] <- pr$V_qr[, , tau_idx]

    # Invert variance matrix
    V_inv <- tryCatch(solve(V_list[[j]]), error = function(e) NULL)
    if (is.null(V_inv)) {
      return(NULL)
    }
    V_inv_list[[j]] <- V_inv
  }

  # Check for any NA values
  if (any(is.na(betas))) {
    return(NULL)
  }

  # Minimum distance estimator: beta_md = (sum V_i^{-1})^{-1} sum V_i^{-1} beta_i
  sum_V_inv <- Reduce(`+`, V_inv_list)
  sum_V_inv_beta <- Reduce(`+`, lapply(seq_len(n), function(j) V_inv_list[[j]] %*% betas[j, ]))

  V_md <- tryCatch(solve(sum_V_inv), error = function(e) NULL)
  if (is.null(V_md)) {
    return(NULL)
  }

  beta_md <- as.vector(V_md %*% sum_V_inv_beta)
  beta_md_se <- sqrt(diag(V_md))

  # S statistic: sum_i (beta_i - beta_md)' V_i^{-1} (beta_i - beta_md)
  S <- 0
  for (j in seq_len(n)) {
    diff_j <- betas[j, ] - beta_md
    S <- S + as.numeric(t(diff_j) %*% V_inv_list[[j]] %*% diff_j)
  }

  # D statistic: standardized version for large n, T asymptotics
  # D = (S - k(n-1)) / sqrt(2 * k * (n-1))
  df <- k * (n - 1)
  D <- (S - df) / sqrt(2 * df)

  list(S = S, D = D, beta_md = beta_md, beta_md_se = beta_md_se)
}


#' @keywords internal
.compute_ols_test <- function(panel_results, valid_idx, k) {
  # Compute joint S and D statistics for OLS

  n <- length(valid_idx)

  # Collect beta estimates and variances
  betas <- matrix(NA, nrow = n, ncol = k)
  V_inv_list <- vector("list", n)

  for (j in seq_along(valid_idx)) {
    pr <- panel_results[[valid_idx[j]]]
    betas[j, ] <- pr$beta_ols

    V_inv <- tryCatch(solve(pr$V_ols), error = function(e) NULL)
    if (is.null(V_inv)) {
      # Use identity as fallback
      V_inv <- diag(k)
    }
    V_inv_list[[j]] <- V_inv
  }

  # Minimum distance estimator
  sum_V_inv <- Reduce(`+`, V_inv_list)
  sum_V_inv_beta <- Reduce(`+`, lapply(seq_len(n), function(j) V_inv_list[[j]] %*% betas[j, ]))

  V_md <- tryCatch(solve(sum_V_inv), error = function(e) diag(k))
  beta_md <- as.vector(V_md %*% sum_V_inv_beta)

  # S statistic
  S <- 0
  for (j in seq_len(n)) {
    diff_j <- betas[j, ] - beta_md
    S <- S + as.numeric(t(diff_j) %*% V_inv_list[[j]] %*% diff_j)
  }

  # D statistic
  df <- k * (n - 1)
  D <- (S - df) / sqrt(2 * df)

  list(S = S, D = D)
}


#' @keywords internal
.compute_marginal_test <- function(panel_results, valid_idx, tau_idx, var_idx) {
  # Compute marginal S and D statistics for a single variable

  n <- length(valid_idx)

  # Collect scalar beta estimates and variances for this variable
  betas <- numeric(n)
  vars <- numeric(n)

  for (j in seq_along(valid_idx)) {
    pr <- panel_results[[valid_idx[j]]]
    betas[j] <- pr$beta_qr[tau_idx, var_idx]
    vars[j] <- pr$V_qr[var_idx, var_idx, tau_idx]
  }

  # Check for valid values
  if (any(is.na(betas)) || any(is.na(vars)) || any(vars <= 0)) {
    return(NULL)
  }

  v_inv <- 1 / vars

  # MD estimator (scalar case)
  beta_md <- sum(v_inv * betas) / sum(v_inv)

  # S statistic
  S <- sum(v_inv * (betas - beta_md)^2)

  # D statistic (k = 1)
  df <- n - 1
  D <- (S - df) / sqrt(2 * df)

  list(S = S, D = D)
}
