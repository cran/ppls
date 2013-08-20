\name{penalized.pls.cv}
\alias{penalized.pls.cv}
\title{Cross-validation for Penalized PLS}
\description{
  Computes the cross-validated error of penalized PLS for different values of
  lambda and  components, and returns the parameter values and
  coefficients for the optimal model.
}
\usage{
penalized.pls.cv(X, y, P, lambda, ncomp, k, kernel,scale)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of input data}
  \item{y}{vector of responses}
  \item{P}{Penalty matrix. For the default value \code{P=NULL}, no penalty term is used,
    i.e. ordinary PLS is computed.}
  \item{lambda}{vector of candidate parameters lambda for the amount of
    penalization. Default value is 1 }
  \item{ncomp}{Number of penalized PLS components to be
    computed. Default value is min(nrow(X)-1,ncol(X))}
  \item{k}{the number of splits in  \code{k}-fold cross-validation. Default value is \code{k}=5. }
  \item{kernel}{Logical value. If \code{kernel=TRUE}, the kernelized version of
    penalized PLS is computed. Default value is \code{kernel=FALSE}}
  \item{scale}{logical value. If scale=TRUE, the X variables are
    standardized to have unit variance. Default value is FALSE}
}

\value{
   \item{error.cv}{matrix of cross-validated errors. The rows correspond to the values of lambda, the columns correspond to the number of components.}
   \item{lambda}{vector of candidate parameters lambda for the amount of
    penalization}
  \item{lambda.opt }{Optimal value of lambda}
  \item{index.lambda}{Index for the optimal value of lambda}
  \item{ncomp.opt}{Optimal number of penalized PLS components}
  \item{min.ppls}{Cross-validated error for the optimal penalized PLS
    solution}
  \item{intercept}{Intercept for the optimal model, computed on the
    whole data set}
  \item{coefficients}{Regression coefficients for the optimal model,
    computed on the whole data set}
  \item{coefficients.jackknife}{array of regression coefficients for each cross-validation run and each parameter setting. The dimension is \code{ncol(X) x ncomp x length(lambda) x k}. This result can be used to estimate the variance of the regression coefficients.}

}
\references{N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications
to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems 94, 60 - 69. \url{http://dx.doi.org/10.1016/j.chemolab.2008.06.009}}
\author{Nicole Kraemer}

\seealso{\code{\link{ppls.splines.cv}},\code{\link{penalized.pls}},\code{\link{new.penalized.pls}},\code{\link{jack.ppls}}}
\examples{
# the penalty term in this example does not make much
# sense
X<-matrix(rnorm(20*100),ncol=20)
y<-rnorm(rnorm(100))
P<-Penalty.matrix(m=20)
pen.pls<-penalized.pls.cv(X,y,lambda=c(0,1,10),P=P,ncomp=10,kernel=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
