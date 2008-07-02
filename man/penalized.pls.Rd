\name{penalized.pls}
\alias{penalized.pls}
\title{Penalized Partial Least Squares}
\description{
  computes the regression coefficients for Penalized Partial Least Squares.
}
\usage{
penalized.pls(X, y, P, ncomp, kernel, scale ,blocks,select)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of input data}
  \item{y}{vector of response data}
  \item{P}{penalty matrix. Default value is P=NULL, i.e. no penalization is
  used}
  \item{ncomp}{number of components, default value is min(ncol(p),nrow(X)-1)}
  \item{kernel}{logical value. If kernel=TRUE, penalized PLS is computed
  based on the kernel algorithm. Default value is kernel=FALSE}
  \item{scale}{logical value. If scale=TRUE, the X variables are
    standardized to have unit variance. Default value is FALSE}
  \item{blocks}{vector of length \code{ncol(X)} that encodes a block
    structure of the data. Default value is \code{1:ncol(X)}. See below
    for more details.}
  \item{select}{logical variable. If \code{logical=TRUE}, block-wise
    variable selection is applied. Default value is \code{FALSE}. See
    below for more details.}
}
\details{The regression coefficients can be computed in two different
  but equivalent ways. The first one is the extension of the classical
  NIPALS algorithm for PLS (which corresponds to \code{kernel=FALSE}),
  and the second one is based on a kernel representation. The latter
  method is in general faster if the number of observations is small
  compared to the number of variables.
  Note that \code{P=NULL} corresponds to Partial Least Squares without
  penalization.
  In addition, it is possible to select blocks of variables in each
  iteration step of penalized PLS. The block structure is encoded in the
  vector \code{blocks} of length \code{ncol(X)} that has the form
  1,..,1,2,..,2,3,...,3,... . If \code{select=TRUE}, the algorithm
  select the weight vector with maximal penalized covariance under the
  constraint that only a single block in the weight vector is
  non-zero. This strategy is used for the combination of penalized PLS
  and B-splines transformations.
}
\value{
  \item{intercept}{vector of length \code{ncomp}. The ith entry corresponds to
    the intercept for penalized PLS with i components}
  \item{coefficients}{matrix of dimension ncol(X) x \code{ncomp}. The ith
    column corresponds to the regressions coefficients for penalized
    PLS with i components}
}
\references{N. Kr\"amer, A.-L. Boulsteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications
to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, doi: 10.1016/j.chemolab.2008.06.009}
\author{Nicole Kr\"amer}
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
Z=dummy$Z # transformed X data
Ztest=dummy$Ztest # transformed Xtest data
size=dummy$sizeZ # size of the transformed data
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
