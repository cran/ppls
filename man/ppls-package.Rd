\name{ppls-package}
\alias{ppls-package}
\alias{ppls}
\docType{package}
\title{ppls - Penalized Partial Least Squares
}
\description{
Partial Least Squares in combination with a  penalization term.
}
\details{
This package contains functions to estimate linear and nonlinear
regression methods with Penalized Partial Least Squares.

Partial Leasts Squares (PLS) is a regression method that constructs latent
components \eqn{Xw} from the data \eqn{X} with maximal covariance to a
response \eqn{y}. The
components are then used in a least-squares fit instead of
\eqn{X}. For a quadratic penalty term on \eqn{w}, Penalized Partial Least Squares constructs latent components
that maximize the penalized covariance. Applications include the
estimation of generalized additive models and functional data. More
details can be found in Kraemer, Boulesteix, and Tutz (2007).  
}
\author{Nicole Kraemer <nkraemer@cs.tu-berlin.de>
}
\references{
N. Kraemer, A.-L. Boulesteix, G. Tutz (2007) "Penalized
  Partial Least Squares with Applications to B-Splines Transformations
  and Functional Data", preprint

  available at \url{http://ml.cs.tu-berlin.de/~nkraemer/publications.html}
}

\keyword{ package }

