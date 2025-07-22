#' Default Algorithm for Penalized Partial Least Squares (PPLS)
#'
#' Computes the regression coefficients using the standard (NIPALS-based) version of Penalized PLS. This function is typically called internally by \code{\link{penalized.pls}}.
#' @name penalized.pls
#' @param X A numeric matrix of centered (and optionally scaled) predictor variables.
#' @param y A numeric vector of centered response values.
#' @param M Optional transformation matrix derived from the penalty matrix: \eqn{M = (I + P)^{-1}}. If \code{NULL}, no penalization is applied.
#' @param ncomp Integer. Number of PLS components to compute.
#'
#' @details
#' The method is based on iteratively computing latent directions that maximize the covariance with the response \code{y}. At each step:
#' \itemize{
#'   \item A weight vector \eqn{w} is computed as \eqn{w = M X^\top y} (if penalization is used).
#'   \item The latent component \eqn{t = X w} is extracted and normalized.
#'   \item The matrix \code{X} is deflated orthogonally with respect to \code{t}.
#' }
#'
#' The final regression coefficients are computed via a triangular system using the bidiagonal matrix \eqn{R = T^\top X W}, and backsolving:
#' \deqn{\beta = W L (T^\top y),}
#' where \eqn{L = R^{-1}}.
#'
#' @returns A list with:
#' \describe{
#'   \item{coefficients}{A matrix of size \code{ncol(X) x ncomp}, each column containing the regression coefficients for the first \eqn{i} components.}
#' }
#'
#' @seealso \code{\link{penalized.pls}}, \code{\link{penalized.pls.kernel}}, \code{\link{normalize.vector}}
#'
#' @keywords internal
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(20 * 50), nrow = 50)
#' y <- rnorm(50)
#' M <- diag(ncol(X))  # No penalty
#' coef <- penalized.pls.default(scale(X, TRUE, FALSE), scale(y, TRUE, FALSE),
#'   M, ncomp = 3)$coefficients
#' coef[, 1]  # coefficients for 1st component
#'
#' @importFrom stats var
#' @export
penalized.pls.default <- function(X, y, M = NULL, ncomp) {
  if (var(y) == 0) {
    coefficients <- matrix(0, ncol(X), ncomp) # if y is constant, all coefficients are zero
    return(list(coefficients = coefficients))
  }
  X0 <- X
  WW <- matrix(NA, ncol(X), ncomp)
  TT <- matrix(NA, nrow(X), ncomp)
  for (i in 1:ncomp) {
    ww <- t(X) %*% y
    if (!is.null(M)) ww <- M %*% ww
    # ww = normalize.vector(ww)
    tt <- X %*% ww
    tt <- normalize.vector(tt)
    WW[, i] <- ww
    TT[, i] <- tt
    X <- X - tt %*% (t(tt) %*% X)
  }
  B <- matrix(, ncol(X), ncomp)
  RR <- matrix(t(TT) %*% X0 %*% WW, nrow = ncomp) # the bidiagonal matrix
  RR[row(RR) > col(RR)] <- 0 # inserted for numerical stability
  RR[row(RR) < col(RR) - 1] <- 0 # inserted for numerical stability
  LL <- backsolve(RR, diag(ncomp))
  B <- matrix(, ncol(X), ncomp)
  for (i in 1:ncomp) {
    Li <- LL[1:i, 1:i, drop = FALSE]
    B[, i] <- WW[, 1:i, drop = FALSE] %*%
      (Li %*% (t(TT[, 1:i, drop = FALSE]) %*% y))
  }
  coefficients <- B
  return(list(coefficients = coefficients))
}
