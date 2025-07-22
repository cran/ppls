#' Jackknife Estimation for Penalized PLS Coefficients
#'
#' This function computes jackknife estimates (mean and covariance) of the regression coefficients obtained from a cross-validated Penalized Partial Least Squares (PPLS) model.
#'
#' @name jack.ppls
#' @param ppls.object An object returned by \code{\link{penalized.pls.cv}}. Must contain the array \code{coefficients.jackknife} as well as fields \code{lambda}, \code{ncomp.opt}, and \code{index.lambda}.
#' @param ncomp Integer. Number of PLS components to use. Default is \code{ppls.object$ncomp.opt}.
#' @param index.lambda Integer. Index of the penalization parameter \code{lambda}. Default is \code{ppls.object$index.lambda}.
#'
#' @details
#' The jackknife estimates are computed using the array of regression coefficients obtained in each cross-validation fold. The function returns both the mean coefficients and the associated variance-covariance matrix.
#'
#' If the requested number of components \code{ncomp} or the lambda index \code{index.lambda} exceeds the available dimensions of the \code{coefficients.jackknife} array, they are adjusted to their maximum allowable values, with a message.
#'
#' Note: This jackknife procedure is not discussed in Kraemer et al. (2008), but it is useful for statistical inference, such as confidence intervals or hypothesis tests.
#'
#' @returns An object of class \code{"mypls"}, which is a list containing:
#' \describe{
#'   \item{coefficients}{The mean regression coefficients across cross-validation splits.}
#'   \item{covariance}{The estimated covariance matrix of the coefficients.}
#'   \item{k}{Number of cross-validation folds used.}
#'   \item{ncomp}{Number of components used in estimation.}
#'   \item{index.lambda}{Index of the lambda value used.}
#' }
#'
#' @seealso \code{\link{penalized.pls.cv}}, \code{\link{coef.mypls}}, \code{\link{vcov.mypls}}, \code{\link{ttest.ppls}}
#'
#' @references
#' N. Kraemer, A.-L. Boulesteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, 94(1), 60–69. \doi{10.1016/j.chemolab.2008.06.009}
#'
#' @examples
#' data(cookie)  # load example data
#' X <- as.matrix(cookie$NIR)  # NIR spectra
#' y <- cookie$constituents$fat    # extract one constituent
#'
#' pls.object <- penalized.pls.cv(X, y, ncomp = 10, kernel = TRUE)
#' my.jack <- jack.ppls(pls.object)
#' coef(my.jack)
#' vcov(my.jack)
#' @importFrom stats cov
#' @export

jack.ppls <- function(ppls.object,
                      ncomp = ppls.object$ncomp.opt,
                      index.lambda = ppls.object$index.lambda) {
  mydims <- dim(ppls.object$coefficients.jackknife)
  ncomp.cv <- mydims[2]
  k <- mydims[4]
  if (ncomp > ncomp.cv) {
    ncomp = ncomp.cv
    message(paste("ncomp is too large and set to ", ncomp, ".\n"))
  }
  if (index.lambda > length(ppls.object$lambda)) {
    index.lambda <- length(ppls.object$lambda)
    message(paste("index of lambda is too large and set to ", index.lambda, " .\n"))
  }
  index.lambda <- ppls.object$index.lambda
  p <- length(ppls.object$coefficients)
  l <- length(ppls.object$lambda)
  coefficients <- ppls.object$coefficients.jackknife[, ncomp, index.lambda, ]
  mean.ppls <- apply(coefficients, 1, mean)
  vcov.ppls <- cov(t(coefficients)) * ((k - 1)^2) / k
  outlist <- list(
    coefficients = mean.ppls,
    covariance = vcov.ppls,
    k = k,
    ncomp = ncomp,
    index.lambda = index.lambda
  )
  class(outlist) <- 'mypls'
  return(outlist)
}
