\name{X2s}
\alias{X2s}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Nolinear Transformation via B-splines}
\description{This function transforms each column of a matrix using a
  set of B-spline functions.
}
\usage{
X2s(X, Xtest = NULL, deg = 3, nknot = NULL,reduce.knots=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{data matrix}
  \item{Xtest}{optional matrix of test data}
  \item{deg}{degree of the splines. Default value is 3}
  \item{nknot}{vector of length \code{ncol(X)}. The \eqn{j}th entry
    determines the number of knots to be used for the \eqn{j}th column
    of \code{X}. Default value is \code{rep(20,ncol(X))}. }
  \item{reduce.knots}{Logical variable. If \code{TRUE}, the function
    assures that there the transformed data does not contain a constant
    column. See below for more details. Default value is \code{FALSE}. }
}
\details{Each column of the matrix \code{X} represents one variable. For
  each variable, we consider the set of B-splines functions
  \eqn{\phi_1,...,\phi_K} that are determined by the degree  \code{deg}
  of the splines and the number \code{nknot} of knots. The knots are
  equidistantly based on the range of the variable. The data  and -- if
  available -- the test data is the transformed nonlinearly using the
  B-splines function.
  For a large amount of knots, it is possible that some columns of the
transformed matrix \code{Z} only contain zeroes. If this is the case for
one variable and if \code{reduce.knots=TRUE}, the amount of knots is reduced until this phenomenon does
not occur anymore. Note that the
  penalized PLS algorithm runs correctly for constant columns in
  \code{Z}, unless you scale the
  columns of the data.
}
\value{
  \item{Z}{matrix of transformed data}
  \item{Ztest}{matrix of test data, if provided. Otherwise, the
    transformed training data is returned. }
  \item{sizeZ}{vector of length ncol(X). Each component contains the
    number of basis functions for each column of X.}
  }
  \note{Depending on the degrees of the splines - there
must be minimum number of knots. If \code{nknot} contains too few knots,
the function automatically increases the number. 
}
\references{C. de Boor (1978) \emph{A Practical Guide to Splines}, Springer.

  N. Kr\"amer, A.-L. Boulsteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications
to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, 94, 60 - 69.
}
\author{Nicole Kr\"amer}
\seealso{\code{\link{ppls.splines.cv}}}
\examples{
X<-matrix(rnorm(100),ncol=5)
Xtest<-matrix(rnorm(300),ncol=5)
dummy<-X2s(X,Xtest)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
