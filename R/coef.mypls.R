#' Extract Regression Coefficients from a mypls Object
#'
#' Returns the regression coefficients (without intercept) from an object of class \code{mypls}, typically produced by the function \code{\link{jack.ppls}}.
#'
#' @name coef.mypls
#' @param object An object of class \code{mypls}, which must contain the elements \code{coefficients} and \code{covariance}.
#' @param ... Additional arguments passed to methods (currently unused).
#'
#' @returns A numeric vector containing the regression coefficients corresponding to the penalized PLS model.
#'
#' @details This method returns the vector of regression coefficients associated with the penalized PLS fit stored in the \code{mypls} object. These coefficients can be used together with the variance-covariance matrix returned by \code{\link{vcov.mypls}} to construct confidence intervals or hypothesis tests.
#'
#' @seealso \code{\link{vcov.mypls}}, \code{\link{jack.ppls}}
#'
#' @references
#' N. Kraemer, A.-L. Boulesteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, 94(1), 60–69. \doi{10.1016/j.chemolab.2008.06.009}
#'
#' @examples
#' n <- 50  # number of observations
#' p <- 5   # number of variables
#' X <- matrix(rnorm(n * p), ncol = p)
#' y <- rnorm(n)
#'
#' pls.object <- penalized.pls.cv(X, y)
#' my.jack <- jack.ppls(pls.object)
#' my.coef <- coef(my.jack)
#' print(my.coef)
#'
#' @export
coef.mypls <- function(object, ...){
    return(object$coefficients)
}
