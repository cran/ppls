#' Cross-Validation for Penalized Partial Least Squares (PPLS)
#'
#' Performs k-fold cross-validation to evaluate and select the optimal penalization parameter \code{lambda} and the number of components \code{ncomp} in a PPLS model.
#'
#' @name penalized.pls
#' @param X A numeric matrix of input data.
#' @param y A numeric response vector.
#' @param P Optional penalty matrix. If \code{NULL}, ordinary PLS is computed (i.e., no penalization).
#' @param lambda A numeric vector of candidate penalty parameters. Default is \code{1}.
#' @param ncomp Integer. Number of PLS components to compute. Default is \code{min(nrow(X) - 1, ncol(X))}.
#' @param k Integer. Number of cross-validation folds. Default is \code{5}.
#' @param kernel Logical. If \code{TRUE}, uses the kernel representation of PPLS. Default is \code{FALSE}.
#' @param scale Logical. If \code{TRUE}, scales predictors in \code{X} to unit variance. Default is \code{FALSE}.
#'
#' @details
#' The function splits the data into \code{k} cross-validation folds, and for each value of \code{lambda} and number of components up to \code{ncomp}, computes the mean squared prediction error.
#'
#' The optimal parameters are selected as those minimizing the prediction error across all folds. Internally, for each fold and \code{lambda} value, the function calls \code{\link{penalized.pls}} to fit the model and \code{\link{new.penalized.pls}} to evaluate predictions.
#'
#' The returned object can be further used for statistical inference (e.g., via jackknife) or prediction.
#'
#' @returns An object of class \code{"mypls"}, a list with the following components:
#' \describe{
#'   \item{error.cv}{A matrix of mean squared errors. Rows correspond to different \code{lambda} values; columns to different numbers of components.}
#'   \item{lambda}{The vector of candidate lambda values.}
#'   \item{lambda.opt}{The lambda value giving the minimum cross-validated error.}
#'   \item{index.lambda}{The index of \code{lambda.opt} in \code{lambda}.}
#'   \item{ncomp.opt}{The optimal number of PLS components.}
#'   \item{min.ppls}{The minimum cross-validated error.}
#'   \item{intercept}{Intercept of the optimal model (fitted on the full dataset).}
#'   \item{coefficients}{Coefficient vector for the optimal model.}
#'   \item{coefficients.jackknife}{An array of shape \code{ncol(X) x ncomp x length(lambda) x k}, containing the coefficients from each CV split and parameter setting.}
#' }
#'
#' @seealso \code{\link{penalized.pls}}, \code{\link{new.penalized.pls}}, \code{\link{jack.ppls}}, \code{\link{ppls.splines.cv}}
#'
#' @references
#' N. Kraemer, A.-L. Boulesteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, 94(1), 60–69. \doi{10.1016/j.chemolab.2008.06.009}
#'
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(20 * 100), ncol = 20)
#' y <- rnorm(100)
#'
#' # Example with no penalty
#' result <- penalized.pls.cv(X, y, lambda = c(0, 1, 10), ncomp = 5)
#' result$lambda.opt
#' result$ncomp.opt
#' result$min.ppls
#'
#' # Using jackknife estimation after CV
#' jack <- jack.ppls(result)
#' coef(jack)
#'
#' @export
penalized.pls.cv <- function(X,
                             y,
                             P = NULL,
                             lambda = 1,
                             ncomp = NULL,
                             k = 5,
                             kernel = FALSE,
                             scale = FALSE) {
  p <- ncol(X)
  n <- nrow(X)
  if (is.null(ncomp)) ncomp <- min(n - 1, ncol(X))
  lambda <- as.vector(lambda)
  if (is.null(P)) lambda = 0
  all.folds <- split(sample(1:n), rep(1:k, length = n))
  ntrain = vector(length = k)
  for (i in 1:k) {
    ntrain[i] <- n - length(all.folds[[i]])
  }
  ntrain.min <- min(ntrain)
  ncomp <- min(ncomp, ntrain.min - 1)
  error.cv <- matrix(0, length(lambda), ncomp)
  coefficients.jackknife <- array(dim = c(p, ncomp, length(lambda), k))
  for (i in seq(k)) {
    omit <- all.folds[[i]]
    Xtrain <- X[-omit, , drop = FALSE]
    ytrain <- y[-omit]
    Xtest <- X[omit, , drop = FALSE]
    ytest <- y[omit]
    for (j in 1:length(lambda)) {
      if (is.null(P) == TRUE) {
        Pj <- NULL
      }
      if (is.null(P) == FALSE) {
        Pj <- lambda[j] * P
      }
      penpls <- penalized.pls(
        Xtrain,
        ytrain,
        P = Pj,
        ncomp,
        kernel = kernel,
        scale = scale
      )
      coefficients.jackknife[, , j, i] <- penpls$coefficients
      error.cv[j, ] <- error.cv[j, ] + length(ytest) * (new.penalized.pls(penpls, Xtest, ytest)$mse)
    }
  }
  error.cv <- error.cv / n
  value1 <- apply(error.cv, 1, min)
  index.lambda <- which.min(value1)
  lambda.opt <- lambda[index.lambda]
  ncomp.opt <- which.min(error.cv[lambda == lambda.opt, ])
  min.ppls <- min(value1)
  if (is.null(P)) {
    P.opt = NULL
  }
  if (!is.null(P)) {
    P.opt <- lambda.opt * P
  }
  ppls <- penalized.pls(X, y, P.opt, ncomp = ncomp.opt, kernel)
  intercept <- ppls$intercept[ncomp.opt]
  coefficients <- ppls$coefficients[, ncomp.opt]
  outlist <- list(
    error.cv = error.cv,
    lambda = lambda,
    ncomp = ncomp,
    lambda.opt = lambda.opt,
    index.lambda = index.lambda,
    ncomp.opt = ncomp.opt,
    min.ppls = min.ppls,
    intercept = intercept,
    coefficients = coefficients,
    coefficients.jackknife = coefficients.jackknife
  )
  class(outlist) <- "mypls"
  return(outlist)
}
