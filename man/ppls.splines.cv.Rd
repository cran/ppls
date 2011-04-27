\name{ppls.splines.cv}
\alias{ppls.splines.cv}
\title{Cross-validation for penalized PLS based on Spline Transformations}
\description{
  Computes the nonlinear-regression model for penalized PLS based on
  B-Spline transformations.
}
\usage{
ppls.splines.cv(X, y, lambda, ncomp, degree, order, nknot, k, kernel,
               scale, reduce.knots, select)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of input data}
  \item{y}{vector of response data}
  \item{lambda}{vector of candidate parameters lambda for the penalty
    term. Default value is 1}
  \item{ncomp}{Number of PLS components, default value is
    min(nrow(X)-1,ncol(Z)), where Z denotes the transformed data
    obtained from the function \code{X2s} }
  \item{degree}{Degree of the splines. Default value is 3}
  \item{order}{Order of the differences to be computed for the penalty
    term. Default value is 2.}
  \item{nknot}{number of knots. Default value is 20 for all variables.}
   \item{k}{the number of splits in  \code{k}-fold cross-validation. Default value is \code{k}=5. }
  \item{kernel}{Logical value. If kernel=TRUE, the kernelized version of
    penalized PLS is computed. Default value is kernel=FALSE}
   \item{scale}{logical value. If scale=TRUE, the X variables are
    standardized to have unit variance. Default value is FALSE}
  \item{reduce.knots}{Logical variable. If \code{TRUE}, the function
    assures that there the transformed data does not contain a constant
    column.  Default value is \code{FALSE}. }
  \item{select}{Logical value. If \code{select=TRUE}, the function
    fits only one variable per iteration. }
}
\details{This function computes the cv-optimal nonlinear regression
  model with Penalized Partial Least Squares. In a nutshell, the
  algorithm works as follows. Starting with a generalized additive
model for the columns of \code{X}, each additive component is expanded in terms of a generous
amount of B-Splines basis functions. The basis functions are determined
via their \code{degree} and  \code{nknot}, the number of knots. In order to prevent
overfitting,  the additive  model is estimated via penalized PLS, where
the penalty term penalizes the differences of a specified \code{order}. Consult Kraemer, Boulesteix, and Tutz (2008) for details.

A graphical tool for penalized PLS on splines-transformed data is provided by \code{\link{graphic.ppls.splines}}.
}
\value{
  \item{error.cv}{matrix of cross-validated errors. The rows correspond to the values of lambda, the columns correspond to the number of components.}
  \item{lambda.opt }{Optimal value of lambda}
  \item{ncomp.opt}{Optimal number of penalized PLS components}
  \item{min.ppls}{Cross-validated error for the optimal penalized PLS
    solution}
}
\references{N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications
to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, 94, 60 - 69. \url{http://dx.doi.org/10.1016/j.chemolab.2008.06.009}}
\author{Nicole Kraemer}

\seealso{\code{\link{penalized.pls}},\code{\link{penalized.pls.cv}}, \code{\link{graphic.ppls.splines}}}
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
