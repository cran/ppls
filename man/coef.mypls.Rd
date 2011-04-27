\name{coef.mypls}
\Rdversion{1.1}
\alias{coef.mypls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Regression coefficients
}
\description{
This function returns the regression coefficients of a mypls-object.
}
\usage{
\method{coef}{mypls}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{mypls} that is returned by the function \code{jack.ppls}. Objects of the class \code{mypls} require a slot \code{coefficients} and a slot \code{covariance}.}
\item{...}{additional parameters}
}
\details{
The function returns the regression coefficients (without intercept)  for the model parameters assigned to \code{jack.ppls}. Together with the covariance matrix returned by  \code{\link{vcov.mypls}}, it is possible to construct confidence intercals or tests.
}
\value{
regression coefficients.
}

\author{
Nicole Kraemer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{vcov.mypls}},\code{\link{jack.ppls}}
}
\examples{

n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)


pls.object<-penalized.pls.cv(X,y)
my.jack<-jack.ppls(pls.object)
mycoed<-coef(my.jack)
}

\keyword{models}
