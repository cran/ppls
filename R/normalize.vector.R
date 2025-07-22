#' Normalize a Numeric Vector to Unit Length
#'
#' Returns the input vector normalized to have unit Euclidean norm (i.e., length equal to 1).
#'
#' @name normalize.vector
#' @param v A numeric vector.
#'
#' @returns A numeric vector of the same length as \code{v}, with unit norm.
#'
#' @details
#' This function performs:
#' \deqn{v_\text{normalized} = \frac{v}{\sqrt{\sum v_i^2}}}
#' It is primarily used to normalize weight vectors or component directions in Partial Least Squares algorithms.
#'
#' Note: If the input vector has zero norm, the function returns \code{NaN} due to division by zero.
#'
#' @seealso \code{\link{penalized.pls}}, \code{\link{penalized.pls.default}}, \code{\link{penalized.pls.kernel}}
#'
#' @examples
#' v <- c(3, 4)
#' normalize.vector(v)  # returns c(0.6, 0.8)
#'
#' v2 <- rnorm(10)
#' sqrt(sum(normalize.vector(v2)^2))  # should be 1
#'
#' @export
normalize.vector <- function(v) {
  v / sqrt(sum(v^2))
}
