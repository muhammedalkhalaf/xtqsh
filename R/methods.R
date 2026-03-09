#' Print Method for xtqsh Objects
#'
#' @param x An object of class \code{"xtqsh"}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.xtqsh <- function(x, ...) {
  cat("\n")
  cat("Quantile Regression Slope Homogeneity Test\n")
  cat("===========================================\n")
  cat("H0: beta_1(tau) = beta_2(tau) = ... = beta_n(tau)  (slope homogeneity)\n")
  cat("H1: beta_i(tau) != beta_j(tau) for some i != j    (heterogeneity)\n")
  cat("-------------------------------------------\n")
  cat("Dep. variable  :", x$depvar, "\n")
  cat("Covariates     :", paste(x$indepvars, collapse = ", "), "\n")
  cat("Observations   :", x$n_obs, "\n")
  cat("Panels         :", x$valid_panels, "of", x$n_panels, "\n")
  cat("Bandwidth      :", x$bw, "\n")
  if (x$hac > 0) {
    cat("HAC lags       :", x$hac, "\n")
  }
  cat("===========================================\n\n")

  invisible(x)
}


#' Summary Method for xtqsh Objects
#'
#' @param object An object of class \code{"xtqsh"}.
#' @param ... Additional arguments (ignored).
#'
#' @return An object of class \code{"summary.xtqsh"}.
#'
#' @export
summary.xtqsh <- function(object, ...) {
  result <- object
  class(result) <- "summary.xtqsh"
  result
}


#' Print Method for summary.xtqsh Objects
#'
#' @param x An object of class \code{"summary.xtqsh"}.
#' @param digits Number of significant digits. Default is 4.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.summary.xtqsh <- function(x, digits = 4, ...) {
  cat("\n")
  cat(rep("=", 78), "\n", sep = "")
  cat("  Quantile Regression Slope Homogeneity Test\n")
  cat(rep("=", 78), "\n", sep = "")
  cat("  H0: beta_1(tau) = beta_2(tau) = ... = beta_n(tau)  (slope homogeneity)\n")
  cat("  H1: beta_i(tau) != beta_j(tau) for some i != j    (heterogeneity)\n")
  cat(rep("-", 78), "\n", sep = "")
  cat("  Dep. variable  :", x$depvar, "\n")
  cat("  Covariates     :", paste(x$indepvars, collapse = ", "), "\n")
  cat("  Observations   :", x$n_obs, "\n")
  cat(rep("=", 78), "\n\n", sep = "")

  # Joint test results
  cat(rep("-", 78), "\n", sep = "")
  cat("  JOINT SLOPE HOMOGENEITY TEST (H0: all slopes equal across panels)\n")
  cat(rep("-", 78), "\n", sep = "")

  # Header
  cat(sprintf("  %8s  %14s  %10s  %14s  %10s  %6s\n",
              "Quantile", "S-hat stat", "p(S-hat)", "D-hat stat", "p(D-hat)", ""))
  cat(rep("-", 78), "\n", sep = "")

  # OLS/Mean row
  cat(sprintf("  %8s  %14.2f  ", "Mean", x$S_ols))
  .print_pval(x$pval_S_ols, digits)
  cat(sprintf("  %14.3f  ", x$D_ols))
  .print_pval(x$pval_D_ols, digits)
  .print_stars(x$pval_D_ols)
  cat("\n")

  cat("  ", rep("-", 74), "\n", sep = "")

  # QR rows
  for (ti in seq_along(x$tau)) {
    tauval <- x$tau[ti]
    S_val <- x$S[1, ti]
    D_val <- x$D[1, ti]
    pS_val <- x$pval_S[1, ti]
    pD_val <- x$pval_D[1, ti]

    cat(sprintf("  tau = %4.2f  ", tauval))

    if (is.na(S_val)) {
      cat(sprintf("%14s  %10s  %14s  %10s\n", "---", "---", "---", "---"))
    } else {
      cat(sprintf("%14.2f  ", S_val))
      .print_pval(pS_val, digits)
      cat(sprintf("  %14.3f  ", D_val))
      .print_pval(pD_val, digits)
      .print_stars(pD_val)
      cat("\n")
    }
  }

  cat(rep("-", 78), "\n", sep = "")
  cat("  *** p<0.01, ** p<0.05, * p<0.10\n")
  cat(sprintf("  S-hat ~ chi2(%d) under H0 (T -> inf, n fixed)\n", x$df_S))
  cat("  D-hat ~ N(0,1) under H0 (T, n -> inf)\n")
  cat(sprintf("  Bandwidth: %s  |  Panels: n = %d  |  k = %d\n",
              x$bw, x$valid_panels, x$k))

  # Marginal tests (if available)
  if (x$marginal && !is.null(x$S_marginal)) {
    cat("\n")
    cat(rep("-", 78), "\n", sep = "")
    cat("  MARGINAL SLOPE HOMOGENEITY TESTS (per-variable, H0: beta_j same for all i)\n")
    cat(rep("-", 78), "\n", sep = "")

    for (j in seq_len(x$k)) {
      xvar <- x$indepvars[j]
      cat("\n")
      cat("  -- Variable:", xvar, rep("-", max(0, 55 - nchar(xvar))), "\n", sep = " ")
      cat(sprintf("  %8s  %14s  %10s  %14s  %10s  %6s\n",
                  "Quantile", "S-hat stat", "p(S-hat)", "D-hat stat", "p(D-hat)", ""))
      cat("  ", rep("-", 74), "\n", sep = "")

      for (ti in seq_along(x$tau)) {
        tauval <- x$tau[ti]
        S_val <- x$S_marginal[j, ti]
        D_val <- x$D_marginal[j, ti]
        pS_val <- x$pval_S_marginal[j, ti]
        pD_val <- x$pval_D_marginal[j, ti]

        cat(sprintf("  tau = %4.2f  ", tauval))

        if (is.na(S_val)) {
          cat(sprintf("%14s  %10s  %14s  %10s\n", "---", "---", "---", "---"))
        } else {
          cat(sprintf("%14.2f  ", S_val))
          .print_pval(pS_val, digits)
          cat(sprintf("  %14.3f  ", D_val))
          .print_pval(pD_val, digits)
          .print_stars(pD_val)
          cat("\n")
        }
      }
    }

    cat(rep("-", 78), "\n", sep = "")
    cat("  *** p<0.01, ** p<0.05, * p<0.10\n")
    cat(sprintf("  Marginal test: k=1, S-hat ~ chi2(%d), D-hat ~ N(0,1)\n",
                x$valid_panels - 1))
  }

  # MD-QR coefficient table
  cat("\n")
  cat(rep("-", 78), "\n", sep = "")
  cat("  MINIMUM DISTANCE QR ESTIMATES  beta_MD(tau) = (sum V_i^-1)^-1 sum V_i^-1 beta_i\n")
  cat(rep("-", 78), "\n", sep = "")

  # Header for coefficients
  cat(sprintf("  %8s", "Quantile"))
  for (j in seq_len(x$k)) {
    xvar <- substr(x$indepvars[j], 1, 12)
    cat(sprintf("  %12s", xvar))
  }
  cat("\n")
  cat(rep("-", 78), "\n", sep = "")

  # Coefficient rows
  for (ti in seq_along(x$tau)) {
    tauval <- x$tau[ti]

    # Estimates
    cat(sprintf("  tau = %4.2f", tauval))
    for (j in seq_len(x$k)) {
      est <- x$beta_md[ti, j]
      if (is.na(est)) {
        cat(sprintf("  %12s", "---"))
      } else {
        cat(sprintf("  %12.4f", est))
      }
    }
    cat("\n")

    # Standard errors
    cat(sprintf("  %8s", ""))
    for (j in seq_len(x$k)) {
      se <- x$beta_md_se[ti, j]
      if (is.na(se)) {
        cat(sprintf("  %12s", ""))
      } else {
        cat(sprintf("  (%10.4f)", se))
      }
    }
    cat("\n")
  }

  cat(rep("-", 78), "\n", sep = "")
  cat("  Standard errors in parentheses\n")
  cat("\n")

  invisible(x)
}


#' @keywords internal
.print_pval <- function(pval, digits = 4) {
  if (is.na(pval)) {
    cat(sprintf("%10s", "---"))
  } else if (pval < 0.01) {
    cat(sprintf("%10.*f", digits, pval))
  } else if (pval < 0.05) {
    cat(sprintf("%10.*f", digits, pval))
  } else {
    cat(sprintf("%10.*f", digits, pval))
  }
}


#' @keywords internal
.print_stars <- function(pval) {
  if (is.na(pval)) {
    cat("")
  } else if (pval < 0.01) {
    cat(" ***")
  } else if (pval < 0.05) {
    cat(" **")
  } else if (pval < 0.10) {
    cat(" *")
  } else {
    cat("")
  }
}


#' Extract Coefficients from xtqsh Objects
#'
#' Returns the minimum distance quantile regression estimates.
#'
#' @param object An object of class \code{"xtqsh"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A matrix of coefficients with rows for each quantile and columns
#'   for each variable.
#'
#' @export
coef.xtqsh <- function(object, ...) {
  object$beta_md
}


#' Extract Variance-Covariance Matrix from xtqsh Objects
#'
#' Returns the variance-covariance matrix of the minimum distance estimator
#' for the median (or middle quantile).
#'
#' @param object An object of class \code{"xtqsh"}.
#' @param tau Quantile for which to return the variance-covariance matrix.
#'   Default is the median (0.50) if available, otherwise the middle quantile.
#' @param ... Additional arguments (ignored).
#'
#' @return A variance-covariance matrix.
#'
#' @export
vcov.xtqsh <- function(object, tau = NULL, ...) {
  if (is.null(tau)) {
    # Default to median or middle quantile
    if (0.5 %in% object$tau) {
      tau_idx <- which(object$tau == 0.5)
    } else {
      tau_idx <- ceiling(length(object$tau) / 2)
    }
  } else {
    tau_idx <- which(object$tau == tau)
    if (length(tau_idx) == 0) {
      stop("Specified tau not found in results")
    }
  }

  # Return diagonal variance matrix based on standard errors
  se <- object$beta_md_se[tau_idx, ]
  V <- diag(se^2)
  colnames(V) <- rownames(V) <- object$indepvars
  V
}
