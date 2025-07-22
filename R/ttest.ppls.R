#' t-Test for Penalized PLS Regression Coefficients
#'
#' Computes two-sided t-tests and p-values for the regression coefficients of a penalized PLS model based on jackknife estimation.
#'
#' @name ttest.ppls
#' @param ppls.object An object returned by \code{\link{penalized.pls.cv}}, containing the jackknife array \code{coefficients.jackknife}.
#' @param ncomp Integer. Number of PLS components to use. Default is \code{ppls.object$ncomp.opt}.
#' @param index.lambda Integer. Index of the penalty parameter \code{lambda} to use. Default is \code{ppls.object$index.lambda}.
#'
#' @details
#' This function calls \code{\link{jack.ppls}} to estimate:
#' \itemize{
#'   \item The mean of the jackknife coefficients (point estimates),
#'   \item The covariance matrix (for standard errors),
#'   \item The degrees of freedom, equal to \code{k - 1}, where \code{k} is the number of cross-validation folds.
#' }
#'
#' It then performs standard two-sided t-tests:
#' \deqn{t_j = \frac{\hat{\beta}_j}{\text{SE}_j}, \quad \text{df} = k - 1}
#' and computes associated p-values.
#'
#' These p-values can be used for variable selection or inference, although they are based on cross-validation folds and should be interpreted with caution.
#'
#' @returns A list with:
#' \describe{
#'   \item{tvalues}{Numeric vector of t-statistics.}
#'   \item{pvalues}{Numeric vector of two-sided p-values.}
#' }
#'
#' @seealso \code{\link{jack.ppls}}, \code{\link{coef.mypls}}, \code{\link{vcov.mypls}}
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(20 * 100), ncol = 20)
#' y <- rnorm(100)
#' result <- penalized.pls.cv(X, y, lambda = c(0, 1), ncomp = 3)
#' tstats <- ttest.ppls(result)
#' print(tstats$pvalues)
#'
#' @import stats
#' @export
ttest.ppls <- function(ppls.object,
                       ncomp = ppls.object$ncomp.opt,
                       index.lambda = ppls.object$index.lambda) {
  jack.object <- jack.ppls(ppls.object, ncomp = ncomp, index.lambda = index.lambda)
  my.mean <- coef(jack.object)
  my.sd <- sqrt(diag(vcov(jack.object)))
  my.df <- jack.object$k - 1
  tvalues <- my.mean / my.sd
  pvalues <- 2 * pt(abs(tvalues), df = my.df, lower.tail = FALSE)
  return(list(pvalues = pvalues, tvalues = tvalues))

}
