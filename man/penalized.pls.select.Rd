\name{penalized.pls.select}
\alias{penalized.pls.select}
\title{Penalized PLS based on NIPALS Algorithm and blockswise variable selection}
\description{
  Internal function that computes the penalized PLS solutions with
  included block-wise variable selection. 
}
\usage{
penalized.pls.select(X, y, M, ncomp,blocks)
}



\arguments{
  \item{X}{matrix of centered and (possibly) scaled input data }
  \item{y}{vector of centered and (possibly) scaled response data}
  \item{M}{matrix that is a transformation of  the penalty term P}
  \item{ncomp}{number of PLS components}
  \item{blocks}{vector of length \code{ncol(X)} that encodes the block
    structure of \code{X}.}
}

\details{This function assumes that the columns of \code{X} and \code{y}
  are centered and - optionally - scaled. The matrix \code{M} is defined
  as the inverse of \eqn{(I + P)} .  The
  computation of the regression coefficients is based on an extension of
  the classical NIPALS algorithm for PLS. Moreover, in each iteration,
  the weight vector is only defined by one block of variables.  For more details, see Kraemer,
  Boulesteix,  and Tutz (2007). 
}

\value{
  \item{coefficients}{Penalized PLS coefficients for all 1,2,...,ncomp components }

}
\references{N. Kraemer, A.-L. Boulesteix, G. Tutz (2007) "Penalized
  Partial Least Squares with Applications to B-Splines Transformations
  and Functional Data", preprint

  available at \url{http://ml.cs.tu-berlin.de/~nkraemer/publications.html}}
\author{Nicole Kraemer}
\note{This is an internal function that is called by penalized.pls.}
\seealso{\code{\link{penalized.pls}}, \code{\link{ppls.splines.cv}}}
\examples{
# this is an internal function}
\keyword{multivariate}
