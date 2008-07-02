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
details can be found in Kr\"amer, Boulesteix, and Tutz (2008).  

The package also contains a data set from Near-Infrared Spectroscopy (Osborne et.al., 1984).
}
\author{Nicole Kr\"amer <nkraemer@cs.tu-berlin.de>
}
\references{N. Kr\"amer, A.-L. Boulsteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications
to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, doi: 10.1016/j.chemolab.2008.06.009


B.G. Osborne, T. Fearn, A.R. Miller, and S. Douglas (1984) \emph{Application of Near-Infrared Reflectance Spectroscopy to Compositional Analysis of Biscuits and Biscuit Dough}. Journal of the Science of Food and Agriculture, 35, pp. 99 - 105.
}

\keyword{ package }

