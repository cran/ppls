#' Cross-Validation for Penalized PLS with Spline-Transformed Predictors
#'
#' Performs cross-validation to select the optimal number of components and penalization parameter for a penalized partial least squares model (PPLS) fitted to spline-transformed predictors.
#'
#' @name ppls.splines.cv
#' @param X A numeric matrix of input predictors.
#' @param y A numeric response vector.
#' @param lambda A numeric vector of penalty parameters. Default is \code{1}.
#' @param ncomp Integer. Maximum number of PLS components. Default is \code{min(nrow(X) - 1, ncol(X))}.
#' @param degree Integer. Degree of B-splines (e.g., 3 for cubic splines). Default is 3.
#' @param order Integer. Order of the differences used in the penalty matrix. Default is 2.
#' @param nknot Integer or vector. Number of knots per variable (before adjustment). If \code{NULL}, defaults to \code{rep(20, ncol(X))}.
#' @param k Number of folds for cross-validation. Default is 5.
#' @param kernel Logical. Whether to use the kernel representation of PPLS. Default is \code{FALSE}.
#' @param scale Logical. Whether to standardize predictors to unit variance. Default is \code{FALSE}.
#' @param reduce.knots Logical. If \code{TRUE}, adaptively reduces the number of knots when overfitting is detected. Default is \code{FALSE}.
#' @param select Logical. If \code{TRUE}, applies block-wise variable selection. Default is \code{FALSE}.
#'
#' @details
#' This function performs the following steps for each cross-validation fold:
#' \enumerate{
#'   \item Transforms predictors using B-spline basis functions via \code{\link{X2s}}.
#'   \item Computes the penalty matrix using \code{\link{Penalty.matrix}}.
#'   \item Fits a penalized PLS model using \code{\link{penalized.pls}} with the given lambda and number of components.
#'   \item Evaluates prediction performance on the test fold using \code{\link{new.penalized.pls}}.
#' }
#'
#' The optimal parameters are those minimizing the average squared prediction error across all folds.
#'
#' @returns A list with the following components:
#' \describe{
#'   \item{error.cv}{Matrix of prediction errors: rows = lambda values, columns = components.}
#'   \item{min.ppls}{The minimum cross-validated error.}
#'   \item{lambda.opt}{Optimal lambda value.}
#'   \item{ncomp.opt}{Optimal number of components.}
#' }
#'
#' @seealso \code{\link{X2s}}, \code{\link{Penalty.matrix}}, \code{\link{penalized.pls}}, \code{\link{penalized.pls.cv}}
#'
#' @references
#' N. Kraemer, A.-L. Boulesteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, 94(1), 60–69. \doi{10.1016/j.chemolab.2008.06.009}
#'
#' @examples
#' # Simulated data
#' set.seed(123)
#' X <- matrix(rnorm(30 * 100), ncol = 30)
#' y <- rnorm(100)
#'
#' # Run CV with 3 lambdas and max 4 components
#' result <- ppls.splines.cv(X, y, lambda = c(1, 10, 100), ncomp = 4)
#' result$lambda.opt
#' result$ncomp.opt
#'
#' @importFrom splines splineDesign
#' @export
ppls.splines.cv <- function(X,
                            y,
                            lambda = 1,
                            ncomp = NULL,
                            degree = 3,
                            order = 2,
                            nknot = NULL,
                            k = 5,
                            kernel = FALSE,
                            scale = FALSE,
                            reduce.knots = FALSE,
                            select = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(ncomp)) ncomp <- min(n - 1, p)
  lambda <- as.vector(lambda)
  all.folds <- split(sample(1:n), rep(1:k, length = n))

  # ensure that ncomp does not exceed the sample size on the cv splits
  ntrain <- vector(length = k)
  for (i in 1:k) {
    ntrain[i] <- n - length(all.folds[[i]])
  }
  ntrain.min <- min(ntrain)
  ncomp <- min(ncomp, ntrain.min - 1)
  #
  error.cv <- matrix(0, length(lambda), ncomp)

  for (i in seq(k)) {
    omit <- all.folds[[i]]
    Xtrain <- X[-omit, , drop = FALSE]
    ytrain <- y[-omit]
    Xtest <- X[omit, , drop = FALSE]
    ytest <- y[omit]

    Z <- X2s(Xtrain, Xtest, degree, nknot, reduce.knots = reduce.knots)

    Ztrain <- Z$Z

    Ztest <- Z$Ztest

    P <- Penalty.matrix(m = Z$sizeZ, order = order)
    blocks <- c()
    for (b in 1:length(Z$sizeZ)) {
      blocks <- c(blocks, rep(b, Z$sizeZ[b]))
    }

    for (j in 1:length(lambda)) {
      penpls <- penalized.pls(
        Ztrain,
        ytrain,
        lambda[j] * P,
        ncomp,
        kernel,
        blocks = blocks,
        select = select,
        scale = scale
      )
      error.cv[j, ] <- error.cv[j, ] + length(ytest) *
        (new.penalized.pls(penpls, Ztest, ytest)$mse)
    }

  }
  #cat(paste("cv completed \n"))
  error.cv <- error.cv / n
  value1 <- apply(error.cv, 1, min)
  lambda.opt <- lambda[which.min(value1)]
  ncomp.opt <- which.min(error.cv[lambda == lambda.opt, ])
  min.ppls <- min(value1)
  return(
    list(
      error.cv = error.cv,
      min.ppls = min.ppls,
      lambda.opt = lambda.opt,
      ncomp.opt = ncomp.opt
    )
  )
}
