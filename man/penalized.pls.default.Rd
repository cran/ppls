\name{penalized.pls.default}
\alias{penalized.pls.default}
\title{Penalized PLS based on NIPALS Algorithm}
\description{
  Internal function that computes the penalized PLS solutions. 
}
\usage{
penalized.pls.default(X, y, M, ncomp = min(nrow(X) - 1, ncol(X)))
}



\arguments{
  \item{X}{A n x p matrix of centered and (possibly) scaled input data.}
  \item{y}{A numeric n-vector of centered and (possibly) scaled response data.}
  \item{M}{A matrix that is a transformation of  the penalty term P.}
  \item{ncomp}{The number of PLS components, default value is the the rank of the centered matrix \code{X}.}
}

\details{This function assumes that the columns of \code{X} and \code{y}
  are centered and - optionally - scaled. The matrix \code{M} is defined
  as the inverse of \eqn{(I + P)} .  The
  computation of the regression coefficients is based on an extension of
  the classical NIPALS algorithm for PLS. If the number of observations
  is small with respect to the number of variables, it is
  computationally more efficient to use the function
  \code{penalized.pls.kernel}. For more details, see Kraemer,
  Boulesteix,  and Tutz (2007). 
}

\value{
  \item{coefficients}{Penalized PLS coefficients for all 1,2,...,\code{ncomp} components }
}
\references{N. Kraemer, A.-L. Boulesteix, G. Tutz (2007) "Penalized
  Partial Least Squares with Applications to B-Splines Transformations
  and Functional Data", preprint

  available at \url{http://ml.cs.tu-berlin.de/~nkraemer/publications.html}}
\author{Nicole Kraemer}
\note{This is an internal function that is called by penalized.pls. 

}
\seealso{\code{\link{penalized.pls}}, \code{\link{penalized.pls.kernel}}}
\examples{
# this is an internal function}
\keyword{multivariate}
