\name{new.penalized.pls}
\alias{new.penalized.pls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction for Penalized Partial Least Squares}
\description{
  Given a penalized.pls. object, and  new data, this function
  predicts the response for all components. 
}
\usage{
new.penalized.pls(ppls, Xtest, ytest = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppls}{Object returned from penalized.pls}
  \item{Xtest}{matrix of new input data}
  \item{ytest}{vector of new response data, optional}
}
\details{
  Penalized.pls simply returns the intercepts and regression
  coefficients for all penalized PLS components up to \code{ncomp} as
  specified in the function \code{penalized.pls}. \code{new.penalized.pls} then computes the estimated response
  based on these regression vectors. If \code{ytest} is given, the mean squared
  error for all components are computed as well.
}
\value{
  \item{ypres}{matrix of responses}
  \item{mse}{vector of mean squared errors, if ytest is provided.}

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