\name{penalized.pls.cv}
\alias{penalized.pls.cv}
\title{Cross-validation for penalized PLS}
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
  \item{scale}{logical value. I scale=TRUE, the X variables are
    standardized to have unit variance. Default value is FALSE}
}
\details{
}
\value{
  \item{lambda.opt }{Optimal value of lambda}
  \item{ncomp.opt}{Optimal number of penalized PLS components}
  \item{min.ppls}{Cross-validated error for the optimal penalized PLS
    solution}
  \item{intercept}{Intercept for the optimal model, computed on the
    whole data set}
  \item{coefficients}{Regression coefficients for the optimal model,
    computed on the whole data set}

}
\references{N. Kraemer, A.-L. Boulesteix, G. Tutz (2007) "Penalized
  Partial Least Squares with Applications to B-Splines Transformations
  and Functional Data", preprint

  available at \url{http://ml.cs.tu-berlin.de/~nkraemer/publications.html}}
\author{Nicole Kraemer}

\seealso{\code{\link{ppls.splines.cv}},\code{\link{penalized.pls}},\code{\link{new.penalized.pls}}}
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