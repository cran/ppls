\name{vcov.mypls}
\Rdversion{1.1}
\alias{vcov.mypls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variance-covariance matrix of the regression coefficients
}
\description{
This function returns the variance-covariance matrix of PLS regression coefficients. 
}
\usage{
\method{vcov}{mypls}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{mypls} that is returned by the function \code{jack.ppls}. Objects of the class \code{mypls} require a slot \code{coefficients} and a slot \code{covariance}.}
\item{...}{additional parameters}
}
\details{
The function returns the variance-covariance matrix for the model parameters assigned to \code{jack.ppls}. Together with the regression coefficients returned by  \code{\link{coef.mypls}}, it is possible to construct confidence intercals or tests.}
\value{
variance-covariance matrix
}

\author{
Nicole Kraemer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{jack.ppls}}, \code{\link{penalized.pls.cv}}
}
\examples{

n<-50 # number of observations
p<-5 # number of variables
X<-matrix(rnorm(n*p),ncol=p)
y<-rnorm(n)


pls.object<-penalized.pls.cv(X,y)
my.jack<-jack.ppls(pls.object)
myvcov<-vcov(my.jack)
}

\keyword{models}
