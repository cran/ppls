\name{new.penalized.pls}
\alias{new.penalized.pls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction for Penalized Partial Least Squares}
\description{
  Given an object returned by the function \code{\link{penalized.pls}}, this function
  predicts the response for new data, with different numbers of components (as specified
when running \code{\link{penalized.pls}}). 
}
\usage{
new.penalized.pls(ppls, Xtest, ytest = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppls}{An object returned from \code{\link{penalized.pls}}}
  \item{Xtest}{A matrix of new input data with the same number of columns as the matrix \code{X}
given as an argument of the function \code{\link{penalized.pls}} to construct \code{ppls}.}
  \item{ytest}{The numeric vector of new response data, optional. Its length must equal the number of
rows of \code{Xtest}.}
}
\details{
  The function \code{\link{penalized.pls}} simply returns the intercepts and regression
  coefficients for all penalized PLS components up to \code{ncomp} which is
  specified as an argument. The function \code{new.penalized.pls} then computes the estimated response
  using these regression vectors. If \code{ytest} is given, the mean squared
  errors for all number of components are computed as well.
}
\value{
  \item{ypred}{The matrix of predicted responses. The number of rows equals the number of rows of \code{Xtest}, whereas
the number of columns equals the \code{ncomp} argument given to the function \code{penalized.pls} while
constructing the object \code{ppls}.
}
  \item{mse}{vector of mean squared errors, if \code{ytest} (giving the true response) is provided.}

}
\references{N. Kraemer, A.-L. Boulesteix, G. Tutz (2007) "Penalized
  Partial Least Squares with Applications to B-Splines Transformations
  and Functional Data", preprint

  available at \url{http://ml.cs.tu-berlin.de/~nkraemer/publications.html}}
\author{Nicole Kraemer}

\seealso{\code{\link{penalized.pls}}, \code{\link{penalized.pls.cv}}, \code{\link{ppls.splines.cv}}}
\examples{
# see also the example for penalised.pls 
X<-matrix(rnorm(50*200),ncol=50)
y<-rnorm(200)
Xtrain<-X[1:100,]
Xtest<-X[101:200,]
ytrain<-y[1:100]
ytest<-X[101:200]
pen.pls<-penalized.pls(Xtrain,ytrain,ncomp=10)
test.error<-new.penalized.pls(pen.pls,Xtest,ytest)$mse
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}