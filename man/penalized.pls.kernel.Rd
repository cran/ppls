\name{penalized.pls.kernel}
\alias{penalized.pls.kernel}
\title{Kernel Penalized PLS}
\description{
  Internal function that computes the penalized PLS solutions based on a
  kernel matrix.
}
\usage{
penalized.pls.kernel(X, y, M, ncomp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of centered and (possibly) scaled input data }
  \item{y}{vector of centered and (possibly) scaled response data}
  \item{M}{matrix that is a transformation of  the penalty term P}
  \item{ncomp}{number of PLS components}
}
\details{This function assumes that the columns of \code{X} and \code{y}
  are centered. The matrix \code{M} is defined as the inverse of 
   \eqn{(I + P)}. 
  The
  computation of the regression coefficients is based on a Kernel
  representation of penalized PLS. If the number of observations
  is large with respect to the number of variables, it is
  computationally more efficient to use the function
  \code{penalized.pls.default}.  For more details, see Kraemer, Boulesteix,  and Tutz (2007).
}
\value{
  \item{coefficients}{Penalized PLS coefficients for all 1,2,...,ncomp compoents }
}
\references{N. Kraemer, A.-L. Boulesteix, G. Tutz (2007) "Penalized
  Partial Least Squares with Applications to B-Splines Transformations
  and Functional Data", preprint

  available at \url{http://ml.cs.tu-berlin.de/~nkraemer/publications.html}}
\author{Nicole Kraemer}
\note{This is an internal function that is called by penalized.pls

}
\seealso{\code{\link{penalized.pls}},\code{\link{penalized.pls.default}}}
\examples{
# this is an internal function
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
