#' Penalty Matrix for Higher Order Differences
#'
#' Computes the block-diagonal penalty matrix penalizing higher-order differences.
#'
#' @name Penalty.matrix
#' @param m Numeric vector indicating sizes of blocks.
#' @param order Integer indicating the order of differences (default is 2).
#'
#' @details
#' For each block of size \code{m[j]}, and default \code{order = 2}, computes:
#' \deqn{v^\top P_j v = \sum_{i=3}^{m[j]} (v_i - 2v_{i-1} + v_{i-2})^2.}
#' The final penalty matrix is block-diagonal composed of blocks \eqn{P_j}.
#'
#' @returns Penalty matrix (numeric matrix) of dimension \code{sum(m)} x \code{sum(m)}.
#'
#' @examples
#' P <- Penalty.matrix(c(6, 4), order = 2)
#'
#' @references
#' Kraemer, N., Boulesteix, A.-L., & Tutz, G. (2008). Penalized Partial Least Squares with Applications to B-Spline Transformations and Functional Data. \emph{Chemometrics and Intelligent Laboratory Systems}, 94, 60-69. https://doi.org/10.1016/j.chemolab.2008.06.009
#'
#' @export
Penalty.matrix <- function(m, order = 2){
  ######################################################
  ## internal function that computes the differences
  Diff.matrix <- function(m, order = 2) {
    # internal function that computes the first order differences
    d.matrix <- function(m) {
      A <- cbind(diag(m - 1), rep(0, m - 1))
      B <- cbind(rep(0, m - 1), -1 * diag(m - 1))
      d <- A + B
      return(d)
    }

    D <- d.matrix(m)
    if (order > 1) {
      for (k in 2:order) {
        D <- d.matrix(m - k + 1) %*% D
      }
    }
    return(D)
  }

  ########################################################################
  # m is a vector that contains the size of the blocks for each penalty term
  p <- length(m)

  start.block = cumsum(m) - m + 1
  end.block = cumsum(m)
  P <- matrix(0, sum(m), sum(m))

  for (i in 1:p) {
    D <- Diff.matrix(m[i], order = order)
    K <- t(D) %*% D
    P[start.block[i]:end.block[i], start.block[i]:end.block[i]] = K
  }
  return(P)

}
