#' Variance-Covariance Matrix for Penalized PLS Coefficients
#'
#' Returns the estimated variance-covariance matrix of the regression coefficients from a jackknife-based PPLS model.
#'
#' @name vcov.mypls
#' @param object An object of class \code{"mypls"}, typically returned by \code{\link{jack.ppls}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' The function retrieves the covariance matrix stored in the \code{object$covariance} field. If this field is \code{NULL}, a warning is issued and the function returns \code{NULL}.
#'
#' This method can be used in conjunction with \code{\link{coef.mypls}} and \code{\link{ttest.ppls}} to conduct inference on the coefficients of a penalized PLS model.
#'
#' @returns A numeric matrix representing the variance-covariance matrix of the regression coefficients, or \code{NULL} if unavailable.
#'
#' @seealso \code{\link{coef.mypls}}, \code{\link{jack.ppls}}, \code{\link{ttest.ppls}}
#'
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(30 * 100), ncol = 30)
#' y <- rnorm(100)
#'
#' ppls.cv <- penalized.pls.cv(X, y, lambda = c(1, 10), ncomp = 3)
#' myjack <- jack.ppls(ppls.cv)
#' Sigma <- vcov(myjack)
#' Sigma[1:5, 1:5]
#'
#' @export
vcov.mypls <- function(object, ...) {
  dummy <- object$covariance
  if (is.null(dummy)) {
    warning("Covariance of regression coefficients is
            not available. Returning NULL.")
  }
  return(dummy)
}
