#' Penalized PLS with Block-wise Variable Selection
#'
#' Computes the regression coefficients of a Penalized Partial Least Squares (PPLS) model using block-wise selection, where each component is restricted to use variables from only one block.
#'
#' @name penalized.pls
#' @param X A numeric matrix of centered (and optionally scaled) predictor variables.
#' @param y A centered numeric response vector.
#' @param M Optional penalty transformation matrix \eqn{M = (I + P)^{-1}}. If \code{NULL}, no penalization is applied.
#' @param ncomp Integer. Number of PLS components to compute.
#' @param blocks An integer vector of length \code{ncol(X)} that defines the block structure of the variables. All variables sharing the same value in \code{blocks} belong to the same block.
#'
#' @details
#' This function implements a sparse selection strategy inspired by sparse or group PLS. At each component iteration, it computes the penalized covariance between \code{X} and \code{y}, and selects the block \code{k} for which the mean squared weight of its variables is maximal:
#' \deqn{\text{score}_k = \frac{1}{|B_k|} \sum_{j \in B_k} w_j^2}
#'
#' Only the weights corresponding to the selected block are retained, and all others are set to zero. The rest of the algorithm follows the classical NIPALS-like PLS with orthogonal deflation.
#'
#' This procedure enhances interpretability by selecting only one block per component, making it suitable for structured variable selection (e.g., grouped predictors).
#'
#' @returns A list with:
#' \describe{
#'   \item{coefficients}{A matrix of size \code{ncol(X) x ncomp}, containing the regression coefficients after block-wise selection.}
#' }
#'
#' @seealso \code{\link{penalized.pls}}, \code{\link{penalized.pls.cv}}, \code{\link{normalize.vector}}
#'
#' @keywords internal
#'
#' @examples
#' set.seed(321)
#' X <- matrix(rnorm(40 * 30), ncol = 40)
#' y <- rnorm(30)
#'
#' # Define 4 blocks of 10 variables each
#' blocks <- rep(1:4, each = 10)
#' result <- penalized.pls.select(X, y, M = NULL, ncomp = 2, blocks = blocks)
#' result$coefficients[, 1]  # Coefficients for first component
#'
#' @importFrom stats var
#' @export

penalized.pls.select <- function(X, y, M = NULL, ncomp, blocks) {
  if (var(y) == 0) {
    coefficients <- matrix(0, ncol(X), ncomp) # if y is constant, all coefficients are zero
    return(list(coefficients = coefficients))
  }

  X0 <- X
  WW <- matrix(, ncol(X), ncomp)
  TT <- matrix(, nrow(X), ncomp)

  for (i in 1:ncomp) {
    ww <- t(X) %*% y
    if (is.null(M) == FALSE) {
      ww <- M %*% ww
    }
    # select optimal blocks
    score <- vector(length = max(blocks))
    for (k in 1:length(score)) {
      score[k] <- mean(ww[blocks == k]^2)
    }
    k.max <- which.max(score)
    w <- ww * 0
    w[blocks == k.max] <- ww[blocks == k.max]
    ww <- normalize.vector(w)
    #
    tt <- X %*% ww
    tt <- normalize.vector(tt)
    WW[, i] <- ww
    TT[, i] <- tt
    X <- X - tt %*% (t(tt) %*% X)
  }
  B <- matrix(, ncol(X), ncomp)
  RR <- matrix(t(TT) %*% X0 %*% WW, nrow = ncomp) # the triangular matrix
  RR[row(RR) > col(RR)] <- 0 # inserted for numerical stability
  LL <- backsolve(RR, diag(ncomp))
  B <- matrix(, ncol(X), ncomp)
  for (i in 1:ncomp) {
    Li <- LL[1:i, 1:i, drop = FALSE]
    B[, i] <- WW[, 1:i, drop = FALSE] %*% (Li %*% (t(TT[, 1:i, drop =
                                                         FALSE]) %*% y))
  }
  coefficients <- B
  return(list(coefficients = coefficients))

}
