#' Predict New Data Using a Penalized PLS Model
#'
#' Given a fitted penalized PLS model and new test data, this function predicts the response for all components. If true response values are provided, it also returns the mean squared error (MSE) for each component.
#'
#' @name penalized.pls
#' @param ppls A fitted penalized PLS model, as returned by \code{\link{penalized.pls}}.
#' @param Xtest A numeric matrix of new input data for prediction.
#' @param ytest Optional. A numeric response vector corresponding to \code{Xtest}, for evaluating prediction error.
#'
#' @details
#' The fitted model \code{ppls} contains intercepts and regression coefficients for each number of components (from 1 to \code{ncomp}). The function computes:
#' \itemize{
#'   \item the matrix of predicted values for each component (as columns),
#'   \item and, if \code{ytest} is provided, a vector of mean squared errors for each component.
#' }
#'
#' The prediction is performed as:
#' \deqn{\hat{y}^{(i)} = X_\text{test} \cdot \beta^{(i)} + \text{intercept}^{(i)},}
#' for each number of components \eqn{i = 1, \ldots, ncomp}.
#'
#' @returns A list containing:
#' \describe{
#'   \item{ypred}{A numeric matrix of predicted responses. Each column corresponds to a different number of PLS components.}
#'   \item{mse}{A numeric vector of mean squared errors, if \code{ytest} is provided. Otherwise \code{NULL}.}
#' }
#'
#' @seealso \code{\link{penalized.pls}}, \code{\link{penalized.pls.cv}}, \code{\link{ppls.splines.cv}}
#'
#' @references
#' N. Kraemer, A.-L. Boulesteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, 94(1), 60–69. \doi{10.1016/j.chemolab.2008.06.009}
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(50 * 200), ncol = 50)
#' y <- rnorm(200)
#'
#' Xtrain <- X[1:100, ]
#' ytrain <- y[1:100]
#' Xtest <- X[101:200, ]
#' ytest <- y[101:200]
#'
#' pen.pls <- penalized.pls(Xtrain, ytrain, ncomp = 10)
#' pred <- new.penalized.pls(pen.pls, Xtest, ytest)
#' head(pred$ypred)
#' pred$mse
#'
#' @export
new.penalized.pls <- function(ppls, Xtest, ytest = NULL) {
  ncomp <- ncol(ppls$coefficients)

  #if (!is.matrix(Xtest)) Xtest=matrix(Xtest,ncol=1)

  ypred <- Xtest %*% ppls$coefficients +
    rep(1, nrow(Xtest)) %*% t(ppls$intercept)

  mse <- NULL
  if (!is.null(ytest)) {
    yres <- ypred - ytest
    mse <- apply(yres^2, 2, mean)
  }

  return(list(ypred = ypred, mse = mse))
}
