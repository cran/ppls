#' Kernel Algorithm for Penalized Partial Least Squares (PPLS)
#'
#' Computes the regression coefficients using the kernel-based version of Penalized PLS, especially useful when the number of predictors exceeds the number of observations (\code{p >> n}).
#'
#' @name penalized.pls
#' @param X A numeric matrix of centered (and optionally scaled) predictors.
#' @param y A centered numeric response vector.
#' @param M Optional matrix transformation \eqn{M = (I + P)^{-1}} derived from a penalty matrix \code{P}. If \code{NULL}, no penalization is applied.
#' @param ncomp Integer. Number of PLS components to compute.
#'
#' @details
#' The kernel PPLS algorithm is based on representing the model in terms of the Gram matrix \eqn{K = X M X^\top} (or simply \eqn{K = X X^\top} if \code{M = NULL}). The algorithm iteratively computes orthogonal latent components \eqn{t_i} in sample space.
#'
#' Steps:
#' \enumerate{
#'   \item Initialize residual vector \eqn{u = y}, then normalize \eqn{t = Ku}.
#'   \item Orthogonalize \eqn{t} with respect to previous components (if needed).
#'   \item Repeat for \code{ncomp} components.
#' }
#'
#' The regression coefficients are recovered as:
#' \deqn{\beta = X^\top A, \quad \text{where } A = UU L (T^\top y),}
#' with \eqn{UU} and \eqn{TT} the matrices of latent vectors and components, and \eqn{L = R^{-1}} the back-solved triangular system.
#'
#' @returns A list with:
#' \describe{
#'   \item{coefficients}{A matrix of size \code{ncol(X) x ncomp}, containing the estimated regression coefficients for each number of components.}
#' }
#'
#' @seealso \code{\link{penalized.pls}}, \code{\link{penalized.pls.default}}, \code{\link{normalize.vector}}
#'
#' @keywords internal
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(100 * 10), nrow = 100)
#' y <- rnorm(100)
#' K <- X %*% t(X)
#' coef <- penalized.pls.kernel(X, y, M = NULL, ncomp = 2)$coefficients
#' head(coef[, 1])  # coefficients for 1st component
#'
#' @importFrom stats var
#' @export
penalized.pls.kernel <- function(X, y, M = NULL, ncomp) {
  if (var(y) == 0) {
    coefficients <- matrix(0, ncol(X), ncomp) # if y is constant, the coefficients are zero
    return(list(coefficients = coefficients))
  }
  n <- nrow(X)
  yhat <- rep(0, n)
  if (is.null(M)) {
    K <- X %*% t(X)
  } else {
    K <- X %*% M %*% t(X)
  }
  UU <- matrix(NA, n, ncomp)
  TT <- matrix(NA, n, ncomp)
  for (i in 1:ncomp) {
    uu <- y - yhat
    uu <- uu / sqrt(sum((K %*% uu) * uu))
    UU[, i] <- uu
    if (i == 1) {
      tt <- K %*% uu
    }
    if (i > 1) {
      TTi <- TT[, 1:(i - 1), drop = FALSE]
      tt <- K %*% uu -  TTi %*% (t(TTi) %*% (K %*% uu))
      if (floor(i / 5) == i / 5) {
        tt <- tt - TTi %*% (t(TTi) %*% tt)
      }
    }
    tt <- normalize.vector(tt)
    TT[, i] <- tt
    yhat <- yhat + sum(tt * y) * tt
  }
  RR <- t(TT) %*% K %*% UU
  RR[row(RR) > col(RR)] <- 0 # inserted for numerical stability
  RR[row(RR) < col(RR) - 1] <- 0 # inserted for numerical stability
  LL <- backsolve(RR, diag(ncomp))
  AA <- matrix(, nrow(X), ncomp)
  for (i in 1:ncomp) {
    Li <- LL[1:i, 1:i, drop = FALSE]
    AA[, i] <- UU[, 1:i, drop = FALSE] %*%
      (Li %*% (t(TT[, 1:i, drop = FALSE]) %*% y))
  }
  coefficients <- t(X) %*% AA
  if (!is.null(M)) {
    coefficients <- M %*% coefficients
  }
  return(list(coefficients = coefficients))
}
