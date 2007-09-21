\name{penalized.pls}
\alias{penalized.pls}
\title{Penalized Partial Least Squares}
\description{
  This function computes the regression coefficients using the Penalized Partial Least Squares methodology
(see Kraemer, Boulesteix and Tutz, 2007).
}
\usage{
penalized.pls(X, y, P = NULL, ncomp = NULL, kernel = FALSE, scale = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A n x p matrix of input data.}
  \item{y}{A numeric n-vector of responses.}
  \item{P}{The penalty matrix. Default value is P=NULL, i.e. no penalization is
  used.}
  \item{ncomp}{The number of penalized PLS components. Default value is \code{min(ncol(p),nrow(X)-1)}.}
  \item{kernel}{If \code{kernel=TRUE}, penalized PLS is computed
  based on the kernel algorithm. Default value is \code{kernel=FALSE}.}
  \item{scale}{If \code{scale=TRUE}, the \code{X} variables are
    standardized to have unit variance. Default value is \code{scale=FALSE}.}
}
\details{The regression coefficients can be computed in two different
  but equivalent ways. The first one is the extension of the classical
  NIPALS algorithm for PLS (which corresponds to \code{kernel=FALSE}),
  and the second one is based on a kernel representation. The latter
  method is in general faster if the number of observations is small
  compared to the number of variables.
  Note that \code{P=NULL} corresponds to Partial Least Squares without penalization.
}
\value{
  \item{intercept}{A vector of length \code{ncomp}. The \eqn{i}th entry corresponds to
    the intercept for penalized PLS with \eqn{i} components.}
  \item{coefficients}{A \code{ncol(X)} x \code{ncomp} matrix. The \eqn{i}th
    column corresponds to the regressions coefficients for penalized
    PLS with \eqn{i} components.}
}
\references{N. Kraemer, A.-L. Boulesteix, G. Tutz (2007) "Penalized
  Partial Least Squares with Applications to B-Splines Transformations
  and Functional Data", preprint

  available at \url{http://ml.cs.tu-berlin.de/~nkraemer/publications.html}}
\author{Nicole Kraemer}
\seealso{\code{\link{new.penalized.pls}},
  \code{\link{penalized.pls.cv}}, \code{\link{ppls.splines.cv}},
  \code{\link{Penalty.matrix}}}
\examples{
## example from the paper ##
# load BOD data
data(BOD)
X<-BOD[,1]
y<-BOD[,2]

Xtest=seq(min(X),max(X),length=200) # generate test data for plot
dummy<-X2s(X,Xtest,deg=3,nknot=20) # transformation of the data
Z<-dummy$Z # transformed X data
Ztest<-dummy$Ztest # transformed Xtest data
size<-dummy$sizeZ # size of the transformed data
P<-Penalty.matrix(size,order=2)  # Penalty matrix
lambda<-200 # amount of penalization
number.comp<-3 # number of components

ppls<-penalized.pls(Z,y,P=lambda*P,ncomp=number.comp) # fit
new.ppls<-new.penalized.pls(ppls,Ztest)$ypred # prediction for test data
## plot fitted values for 2 components
plot(X,y,lwd=3,xlim=range(Xtest))
lines(Xtest,new.ppls[,2])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
