\name{sim.data.ppls}
\alias{sim.data.ppls}
\title{Simulated Data}
\description{
  generates data that can be used for simulations
}
\usage{
sim.data.ppls(ntrain,ntest,stnr,p,a=NULL,b=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ntrain}{number of training observations}
  \item{ntest}{number of test observations}
  \item{stnr}{signal to noise ratio}
  \item{p}{number of predictor variables}
  \item{a}{vector of length 5 that determines the regression problem to
    be simulated}
  \item{b}{vector of length 5 that determines the regression problem to
    be simulated}
}
\details{The matrix of training and test data is drawn from a uniform
  distribution over [-1,1] for each of the \code{p} variables. The response is
  generated via a nonlinear regression model of the form
  \deqn{Y=\sum _{j=1}  ^5 f_j(X_j) + \varepsilon}
    where \eqn{f_j(x)=a_j x + sin(6 b_jx)}. The values of \eqn{a_j} and
    \eqn{b_j} can be specified via \code{a} or \code{b}. If no values
    for \code{a} or \code{b} is given, they are drawn randomly from
    [-1,1]. The variance of the noise term is chosen such that the
    signal-to-noise-ratio equals \code{stnr} on the training data.
}
\value{
  \item{Xtrain}{matrix of size \code{ntrain} x \code{p}}
  \item{ytrain}{vector of lengt \code{ntrain}}
  \item{Xtest}{matrix of size \code{ntest} x \code{p}}
  \item{ytest}{vector of lengt \code{ntest}}
  \item{sigma}{standard deviation of the noise term}
  \item{a}{vector that determines the nonlinear function}
  \item{b}{vector that determines the nonlinear function}
}
\references{N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications
to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, 94, 60 - 69. \url{http://dx.doi.org/10.1016/j.chemolab.2008.06.009}}
\author{Nicole Kraemer}
\seealso{\code{\link{ppls.splines.cv}}}
\examples{
dummy<-sim.data.ppls(ntrain=50,ntest=200,p=16,stnr=16)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
