#' xtqsh: Panel Quantile Regression Slope Homogeneity Test
#'
#' Tests for slope homogeneity in panel quantile regression models.
#' Implements the methodology of Galvao, Juhl, Montes-Rojas, and Olmo (2017)
#' for testing whether slope coefficients are equal across panels at
#' various quantiles of the conditional distribution.
#'
#' @section Main Function:
#' \describe{
#'   \item{\code{\link{xtqsh}}}{Panel quantile regression slope homogeneity test}
#' }
#'
#' @section Test Statistics:
#' The package implements two test statistics:
#'
#' \strong{S-hat statistic:} Under the null hypothesis of slope homogeneity
#' with T (time periods) tending to infinity and n (panels) fixed, the S-hat
#' statistic follows a chi-squared distribution with k(n-1) degrees of freedom,
#' where k is the number of slope coefficients.
#'
#' \strong{D-hat statistic:} Under the null with both n and T tending to
#' infinity, the D-hat statistic is asymptotically standard normal. This
#' provides a standardized test suitable for large panel settings.
#'
#' @section Minimum Distance Estimator:
#' Under slope homogeneity, the minimum distance (MD) quantile regression
#' estimator optimally combines individual panel estimates:
#' \deqn{\hat{\beta}_{MD}(\tau) = \left(\sum_{i=1}^{n} \hat{V}_i^{-1}\right)^{-1}
#' \sum_{i=1}^{n} \hat{V}_i^{-1} \hat{\beta}_i(\tau)}
#'
#' @section Bandwidth Selection:
#' Two methods are available for estimating the sparsity function (inverse of
#' the error density at zero):
#' \itemize{
#'   \item \code{hallsheather}: Hall and Sheather (1988) bandwidth
#'   \item \code{bofinger}: Bofinger (1975) bandwidth
#' }
#'
#' @section References:
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
#' @docType package
#' @name xtqsh-package
#' @aliases xtqsh-package
"_PACKAGE"
