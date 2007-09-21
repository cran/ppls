\name{Penalty.matrix}
\alias{Penalty.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Penalty matrix for higher order differences}
\description{
  This function computes a penalty matrix penalizing the higher order differences.
}
\usage{
Penalty.matrix(m,order = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{A numeric vector. The \eqn{j}th entry gives the size of the
    \eqn{j}th block in the penalty term.}
  \item{order}{The order of the differences. Default value is \code{order}=2.}
}
\details{For the \eqn{j}th entry of the vector \code{m}, and for the default values \code{order}=2, the penalty 
matrix \eqn{P_j} of size \code{m[j]} x \code{m[j]} penalizes the second order
  differences of a vector \eqn{v} of length m[j], that is
  \deqn{v^T P_j v=  \sum_{i=3} ^{m[j]} (\Delta v_i) ^2,}
  where
  \deqn{\Delta v_i= v_i -2 v_{i-1} + v_{i-2}}
  is the second order difference. This definition is easily extended to
  other values of \code{order}. The final penalty matrix \code{P} is a
  block-diagonal matrix with the \eqn{j}th block equal to \eqn{P_j}.  More
  details can be found in Kraemer, Boulesteix and Tutz (2007). 
}
\value{The penalty matrix of size \code{sum(m)} x \code{sum(m)}
 
}
\references{N. Kraemer, A.-L. Boulesteix, G. Tutz (2007) "Penalized
  Partial Least Squares with Applications to B-Splines Transformations
  and Functional Data", preprint

  available at
  \url{http://ml.cs.tu-berlin.de/~nkraemer/publications.html}

  C. de Boor (1978) "A practical guide to splines", Springer.
}
\author{Nicole Kraemer}
\section{Warning}{ All entries of the vector \code{m} must be larger than
  \code{order}, as the notion of kth order differences does not make
  sense for vectors of length \eqn{\le} k. } 
}

\seealso{\code{\link{penalized.pls}}}
\examples{
P<-Penalty.matrix(c(6,4),2)
# a more detailed example can be found under penalized.pls()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}
