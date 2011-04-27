\name{jack.ppls}
\alias{jack.ppls}
\title{Jackknife estimation for PPLS-coefficients}
\description{
  This function computes the mean and the covariance of the regression coefficients of Penalized Partial Least Squares. 
}
\usage{
jack.ppls(ppls.object,ncomp,index.lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppls.object}{an object returned by \code{penalized.pls.cv}}
  \item{ncomp}{integer. The number of components that are used. The default value is the cross-validation optimal number of components.}
  \item{index.lambda}{integer. The index of the penalization intensity, given by the vector \code{lambda} that was provided to \code{penalized.pls.cv}. The default value is the cross-validation optimal index.}
  
}

\value{The function returns an object of class ''ppls''.
   \item{mean.ppls}{The mean of the regression coefficients over all cross-validation splits. This is a vector of length \code{ncol(X)}. Note that in general, this differs from the regression coefficients computed on the whole data set, but if the number of observations is fairly large, the difference should be small.}
   \item{vcov.ppls}{The covariance matrix of the regression coefficients. This is a symmetric matrix of size \code{ncol(X) x ncol(X)}. }
  \item{index.lambda}{Index for the value of lambda that determines the regression coefficient.}
  \item{ncomp}{Number of components that determines the regression coefficients.}
  \item{k}{The number of cross-validation splits. These can be used to construct a t-test for the coefficients.}

}

\details{The function needs an object returned by \code{penalized.pls.cv}. It estimates the mean and the covariance of the regression coefficient (with \code{ncomp} components and penalization intensity indexed by \code{index.lambda}). This is done via a jackknife estimate over the \code{k} cross-validation splits. We remark that this estimation step is not discussed in Kraemer, Boulesteix and Tutz (2008).}
\references{N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications
to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems 94, 60 - 69. \url{http://dx.doi.org/10.1016/j.chemolab.2008.06.009}}
\author{Nicole Kraemer}

\seealso{\code{\link{penalized.pls.cv}},\code{\link{ttest.ppls}}}
\examples{
data(cookie) # load data
X<-as.matrix(cookie[,1:700]) # extract NIR spectra
y<-as.vector(cookie[,701]) # extract one constituent

pls.object<-penalized.pls.cv(X,y,ncomp=10,kernel=TRUE) # PLS without penalization
my.jack<-jack.ppls(pls.object)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
