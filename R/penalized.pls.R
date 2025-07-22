#' Penalized Partial Least Squares (PPLS)
#'
#' Computes the regression coefficients for a Penalized Partial Least Squares (PPLS) model, using either a classical NIPALS algorithm or a kernel-based version. Optionally allows block-wise variable selection.
#'
#' @name penalized.pls
#' @param X A numeric matrix of input data (predictors).
#' @param y A numeric response vector.
#' @param P Optional penalty matrix. If \code{NULL}, no penalization is applied.
#' @param ncomp Integer. Number of PLS components to compute. Default is \code{min(ncol(X), nrow(X) - 1)}.
#' @param kernel Logical. If \code{TRUE}, uses the kernel representation of PPLS (faster if \code{ncol(X)} >> \code{nrow(X)}). Default is \code{FALSE}.
#' @param scale Logical. If \code{TRUE}, each column of \code{X} is standardized to have unit variance. Default is \code{FALSE}.
#' @param blocks Integer vector of length \code{ncol(X)}. Encodes the block structure of the predictors (e.g., \code{c(1,1,2,2,...)}). Default is \code{1:ncol(X)} (each variable is its own block).
#' @param select Logical. If \code{TRUE}, block-wise variable selection is applied in each iteration. Only one block contributes to the latent direction per component. Default is \code{FALSE}.
#'
#' @details
#' This function centers \code{X} and \code{y}, and optionally scales \code{X}, then computes PPLS components using one of:
#' \itemize{
#'   \item the classical NIPALS algorithm (\code{kernel = FALSE}), or
#'   \item the kernel representation (\code{kernel = TRUE}), often faster when \code{p > n} (high-dimensional case).
#' }
#'
#' When a penalty matrix \code{P} is supplied, a transformation \eqn{M = (I + P)^{-1}} is computed internally. The algorithm then maximizes the penalized covariance between \code{Xw} and \code{y}:
#' \deqn{\text{argmax}_w \; \text{Cov}(Xw, y)^2 - \lambda \cdot w^\top P w}
#'
#' The block-wise selection strategy (when \code{select = TRUE}) restricts the weight vector \code{w} at each iteration to be non-zero in a single block, selected greedily.
#'
#' @returns A list with components:
#' \describe{
#'   \item{intercept}{A numeric vector of intercepts for 1 to \code{ncomp} components.}
#'   \item{coefficients}{A numeric matrix of size \code{ncol(X)} x \code{ncomp}, each column being the coefficient vector for the corresponding number of components.}
#' }
#'
#' @seealso \code{\link{penalized.pls.cv}}, \code{\link{new.penalized.pls}}, \code{\link{ppls.splines.cv}}, \code{\link{Penalty.matrix}}
#'
#' @references
#' N. Kraemer, A.-L. Boulesteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, 94(1), 60–69. \doi{10.1016/j.chemolab.2008.06.009}
#'
#' @examples
#' ## Example from Kraemer et al. (2008)
#' data(BOD)
#' X <- BOD[, 1]
#' y <- BOD[, 2]
#'
#' Xtest <- seq(min(X), max(X), length = 200)
#' dummy <- X2s(X, Xtest, deg = 3, nknot = 20)  # Spline transformation
#' Z <- dummy$Z
#' Ztest <- dummy$Ztest
#' size <- dummy$sizeZ
#' P <- Penalty.matrix(size, order = 2)
#' lambda <- 200
#' number.comp <- 3
#'
#' ppls <- penalized.pls(Z, y, P = lambda * P, ncomp = number.comp)
#' new.ppls <- new.penalized.pls(ppls, Ztest)$ypred
#'
#' # Plot fitted values for 2 components
#' plot(X, y, lwd = 3, xlim = range(Xtest))
#' lines(Xtest, new.ppls[, 2], col = "blue")
#'
#' @importFrom stats var
#' @export
penalized.pls <- function(X,
                          y,
                          P = NULL,
                          ncomp = NULL,
                          kernel = FALSE,
                          scale = FALSE,
                          blocks = 1:ncol(X),
                          select = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  y <- as.vector(y)

  meanx <- apply(X, 2, mean)
  meany <- mean(y)

  if (!is.logical(scale) | length(scale) > 1) {
    stop("The parameter 'scale' should be a boolean of length 1,
         either TRUE or FALSE.")
  }
  if (!is.logical(kernel) | length(kernel) > 1) {
    stop("The parameter 'kernel' should be a boolean of length 1,
         either TRUE or FALSE.")
  }

  if (scale) {
    sdx <- sqrt(apply(X, 2, var))
    sdx[sdx == 0] <- 1 # take care of columns with zero variance
  } else {
    sdx <- rep(1, ncol(X))
  }

  if (is.null(ncomp)) ncomp <- min(p, nrow(X) - 1)
  X <- (X - rep(1, n) %*% t(meanx)) / (rep(1, n) %*% t(sdx))
  y <- scale(y, center = TRUE, scale = FALSE)

  M <- NULL
  if (!is.null(P)) {
    M <- solve(P + diag(p))
  }

  if (select) {
    ppls <- penalized.pls.select(X, y, M, ncomp, blocks)
  } else if (kernel) {
    ppls <- penalized.pls.kernel(X, y, M, ncomp)
  } else {
    ppls <- penalized.pls.default(X, y, M, ncomp)
  }

  coefficients <- ppls$coefficients /
    (sdx %*% t(rep(1, ncol(ppls$coefficients))))
  intercept <- rep(meany, ncomp) - t(coefficients) %*% meanx

  return(list(intercept = intercept, coefficients = coefficients))
}
