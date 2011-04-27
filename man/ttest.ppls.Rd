\name{ttest.ppls}
\alias{ttest.ppls}
\title{T-Test for  regression coefficients}
\description{
  This function computes test statistics and p-values for the regression coefficients of Penalized Partial Least Squares. 
}
\usage{
ttest.ppls(ppls.object,ncomp,index.lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppls.object}{an object returned by \code{penalized.pls.cv}}
  \item{ncomp}{integer. The number of components that are used. The default value is the cross-validation optimal number of components.}
  \item{index.lambda}{integer. The index of the penalization intensity, given by the vector \code{lambda} that was provided to \code{penalized.pls.cv}. The default value is the cross-validation optimal index.}
  
}

\value{
   \item{tvalues}{vector of test statistics}
   \item{pvalues}{vector of p-values}

}

\details{We note that neither the distribution of the regression coefficients nor the correct degrees of freedom are known. Hence, the assumptions of the T-Test might not be fulfilled.  We remark that this testing procedure is not discussed in Kraemer, Boulesteix and Tutz (2008).
In general, the p-values need to be corrected in order to account for the multiple testing problem.

}
\references{N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications
to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems 94, 60 - 69. \url{http://dx.doi.org/10.1016/j.chemolab.2008.06.009}}
\author{Nicole Kraemer}

\seealso{\code{\link{penalized.pls.cv}},\code{\link{jack.ppls}}}
\examples{
data(cookie) # load data
X<-as.matrix(cookie[,1:700]) # extract NIR spectra
y<-as.vector(cookie[,701]) # extract one constituent

pls.object<-penalized.pls.cv(X,y,ncomp=10,kernel=TRUE) # PLS without penalization
my.ttest<-ttest.ppls(pls.object) # test for the cv-optimal model

plot(sort(my.ttest$pvalues),type="l",ylab="sorted pvalues") # plot sorted p-values
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{model}
