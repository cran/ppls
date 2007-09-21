\name{ppls.splines.cv}
\alias{ppls.splines.cv}
\title{Cross-validation for penalized PLS based on Splines Transformations}
\description{
  This function computes the nonlinear-regression model for penalized PLS based on
  B-Splines Transformations.
}
\usage{
ppls.splines.cv(X, y, lambda, ncomp, degree, order, nknot, k, kernel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A n x p matrix of input data.}
  \item{y}{A numeric n-vector of responses.}
  \item{lambda}{A numeric vector giving the (positive) candidate values for the parameter lambda 
  determining the amount of penalization. Default value is 1.}
  \item{ncomp}{The number of PLS components, default value is
    \code{min(nrow(X)-1,ncol(Z))}, where \code{Z} denotes the transformed data
    obtained from the function \code{\link{X2s}}.}
  \item{degree}{The degree of the splines. Default value is 3.}
  \item{order}{The order of the differences to be computed for the penalty
    term. Default value is 2.}
  \item{nknot}{The number of knots. Default value is 20 for all variables.}
   \item{k}{The number of splits in  \code{k}-fold cross-validation. Default value is \code{k}=5. }
  \item{kernel}{If \code{kernel=TRUE}, the kernelized version of
  penalized PLS is computed. Default value is \code{kernel=FALSE}.}
}
\details{This function computes the cv-optimal nonlinear regression
  model with Penalized Partial Least Squares. In a nutshell, the
  algorithm works as follows. Starting with a generalized additive
model for the columns of \code{X}, each additive component is expanded in terms of a generous
amount of B-Splines basis functions. The basis functions are determined
via their \code{degree} and  \code{nknot}, the number of knots. In order to prevent
overfitting,  the additive  model is estimated via penalized PLS, where
the penalty term penalizes the differences of a specified \code{order} within each original variable. 
See Kraemer, Boulesteix, and Tutz (2007) for details.
}
\value{
  \item{lambda.opt }{Optimal value of lambda}
  \item{ncomp.opt}{Optimal number of penalized PLS components}
  \item{min.ppls}{Cross-validated error for the optimal penalized PLS
    solution}
}
\references{N. Kraemer, A.-L. Boulesteix, G. Tutz (2007) "Penalized
  Partial Least Squares with Applications to B-Splines Transformations
  and Functional Data", preprint

  available at \url{http://ml.cs.tu-berlin.de/~nkraemer/publications.html}}
\author{Nicole Kraemer}

\seealso{\code{\link{penalized.pls}},\code{\link{penalized.pls.cv}}}
\examples{
# this example does not make much sense, it only illustrates
# how to use the functions properly

X<-matrix(rnorm(100*5),ncol=5)
y<-sin(X[,1]) +X[,2]^2 + rnorm(100)
lambda<-c(0,1,10,100,1000)
cv.result<-ppls.splines.cv(X,y,ncomp=10,k=10,lambda=lambda)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
