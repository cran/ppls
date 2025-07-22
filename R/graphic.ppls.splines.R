#' Plot Penalized PLS Components for Spline-Transformed Data
#'
#' This function applies a nonlinear regression model using penalized Partial Least Squares (PLS) on B-spline transformed variables, then visualizes each additive component.
#'
#' @name graphic.ppls.splines
#' @param X A numeric matrix of input data.
#' @param y A numeric response vector.
#' @param lambda A numeric value for the penalization parameter. Default is \code{NULL}.
#' @param add.data Logical. If \code{TRUE}, the original data points \code{X} and \code{y} are added to the plots. Default is \code{FALSE}.
#' @param select Logical. If \code{TRUE}, the function fits only one block (variable) per iteration (block-wise selection). Default is \code{FALSE}.
#' @param ncomp Integer. Number of PLS components to use. Default is 1.
#' @param deg Integer. Degree of the B-spline basis. Default is 3.
#' @param order Integer. Order of the differences to penalize. Default is 2.
#' @param nknot A numeric vector specifying the number of knots for each variable. Default is \code{NULL}, which uses \code{rep(20, ncol(X))}.
#' @param reduce.knots Logical. If \code{TRUE}, automatically reduces the number of knots for variables leading to constant basis functions. Default is \code{FALSE}.
#' @param kernel Logical. If \code{TRUE}, uses the kernelized version of PPLS. Default is \code{TRUE}.
#' @param window.size A numeric vector of length 2 indicating the number of plots per row and column. Default is \code{c(3, 3)} (3 rows and 3 columns).
#'
#' @details
#' This function first transforms the input data \code{X} and a test grid \code{Xtest} using B-spline basis functions, then fits a penalized PLS model using these transformed variables.
#' Each additive component (i.e., variable effect) is then plotted individually.
#'
#' If \code{add.data = TRUE}, the actual observations are plotted on top of the corresponding fitted component functions. While this can help visualize the fit, note that only the sum of all fitted components approximates \code{y}, and not each component individually.
#'
#' The function is intended for exploratory visualization and should be used after appropriate model selection using, e.g., \code{\link{ppls.splines.cv}}.
#'
#' @returns A numeric vector of regression coefficients for the final penalized PLS model.
#'
#' @seealso \code{\link{ppls.splines.cv}}, \code{\link{X2s}}, \code{\link{penalized.pls}}, \code{\link{Penalty.matrix}}
#'
#' @references
#' N. Kraemer, A.-L. Boulesteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, 94(1), 60–69. \doi{10.1016/j.chemolab.2008.06.009}
#'
#' @examples
#' # Load Boston housing data
#' library(MASS)
#' data(Boston)
#' y <- Boston[, 14]
#' X <- Boston[, -14]
#' X <- X[, -4]  # remove categorical variable
#' X <- as.matrix(X)
#'
#' # Plot with variable selection and original data
#' graphic.ppls.splines(
#'   X, y, lambda = 100, ncomp = 5,
#'   add.data = TRUE, select = TRUE, window.size = c(3, 4)
#' )
#'
#' # Plot without variable selection and without data
#' graphic.ppls.splines(
#'   X, y, lambda = 100, ncomp = 5,
#'   add.data = FALSE, select = FALSE, window.size = c(3, 4)
#' )
#'
#' @import graphics
#' @import MASS
#' @export
graphic.ppls.splines <- function(X,
                                 y,
                                 lambda = NULL,
                                 add.data = FALSE,
                                 select = FALSE,
                                 ncomp = 1,
                                 deg = 3,
                                 order = 2,
                                 nknot = NULL,
                                 reduce.knots = FALSE,
                                 kernel = TRUE,
                                 window.size = c(3, 3)) {
  p <- ncol(X)
  ntest <- 300 # number of test examples
  Xtest <- matrix(, ntest, p)
  for (i in 1:p) {
    Xtest[, i] = seq(min(X[, i]), max(X[, i]), length = 300)
  }
  # transform training and test data
  Z <- X2s(X,
           Xtest,
           deg = deg,
           nknot = nknot,
           reduce.knots = reduce.knots)
  Ztrain <- Z$Z
  Ztest <- Z$Ztest
  sizeZ <- Z$sizeZ
  P <- lambda * Penalty.matrix(m = sizeZ, order = order)
  blocks = c()
  for (b in 1:length(sizeZ)) {
    blocks = c(blocks, rep(b, sizeZ[b]))
  }

  ppls.object <- penalized.pls(
    Ztrain,
    y,
    P = P,
    ncomp = ncomp,
    select = select,
    kernel = kernel,
    blocks = blocks
  )
  ppls.coefficients <- ppls.object$coefficients[, ncomp]
  Ytest <- matrix(, ntest, ncol(X)) # prediction for each additive component
  for (i in 1:ncol(X)) {
    start <- cumsum(c(0, sizeZ))[i] + 1 # start of the ith block
    end <- cumsum(sizeZ)[i] # end of the ith block
    Ytest[, i] <- Ztest[, start:end] %*% ppls.coefficients[start:end]
  }
  # plot the predicted functions
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = window.size)
  for (i in 1:p) {
    plot(
      Xtest[, i],
      Ytest[, i],
      type = "l",
      lwd = 3,
      xlab = "x",
      ylab = "y",
      main = i,
      col = "blue"
    )
    if (add.data == TRUE) {
      lines(X[, i],
            scale(y, scale = FALSE),
            type = "p",
            lwd = 2)
      lines(
        Xtest[, i],
        Ytest[, i],
        type = "l",
        lwd = 3,
        col = "blue"
      )
    }
  }
  return(ppls.coefficients)

}
