\name{graphic.ppls.splines}
\alias{graphic.ppls.splines}
\title{Plots for penalized PLS based on Spline Transformations}
\description{plotting device for penalized PLS on splines transformed variables
}
\usage{
graphic.ppls.splines(X, y, lambda, add.data, select, ncomp, deg, order,
                     nknot, reduce.knots, kernel, window.size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{matrix of input data}
  \item{y}{vector of response data}
\item{add.data}{logical value. If TRUE, the data \code{X} and \code{y} are also plotted. Default is FALSE. See warning below!}
\item{select}{Logical value. If \code{select=TRUE}, the function
    fits only one variable per iteration. Default is FALSE.}
  \item{lambda}{vector of candidate parameters lambda for the penalty
    term. Default value is NULL}
  \item{ncomp}{Number of PLS components, default value is 1 }
  \item{deg}{Degree of the splines. Default value is 3}
  \item{order}{Order of the differences to be computed for the penalty
    term. Default value is 2.}
  \item{nknot}{number of knots. Default value is 20 for all variables.}
  \item{kernel}{Logical value. If kernel=TRUE, the kernelized version of
    penalized PLS is computed. Default value is kernel=TRUE}
  \item{reduce.knots}{Logical variable. If \code{TRUE}, the function
    assures that there the transformed data does not contain a constant
    column.  Default value is \code{FALSE}. }
\item{window.size}{vector of length size 2. Determines the number of plots on one page. Default is c(3,3), that is 3 rows and 3 columns.}
  
}
\details{This function computes a nonlinear regression
  model with Penalized Partial Least Squares using penalized PLS on B-spline transformed variables. The model parameters have to be provided - for proper model selection, we recommend to determine the optiaml parameters with \code{\link{ppls.splines.cv}}.   Consult Kraemer, Boulesteix, and Tutz (2008) for details.

The function plots the additive components for each variable. 

WARNING: If add.data=TRUE, the function also plots the data X and y. While it seems convenient to compare the data (x_j,y) and the fitted functions (x_j,f_j(x_j)), one should keep in mind that only the sum of the fitted functions f_j(x) are an approximation of y.

}
\value{
  \item{ppls.coefficients}{The regression coefficients for the transformed variables.}
}

\references{N. Kraemer, A.-L. Boulsteix, and G. Tutz (2008). \emph{Penalized Partial Least Squares with Applications
to B-Spline Transformations and Functional Data}. Chemometrics and Intelligent Laboratory Systems, 94, 60 - 69. \url{http://dx.doi.org/10.1016/j.chemolab.2008.06.009}}
\author{Nicole Kraemer}

\seealso{\code{\link{ppls.splines.cv}},\code{\link{X2s}}}
\examples{
# ------------------------------------------------------
# load boston housing data

library(MASS)
data(Boston)
y<-Boston[,14]
X<-Boston[,-14]
X<-X[,-4] # remove categorical variable
X<-as.matrix(X)

# ----------------------------------------------------------------------
# plot ppls results for some random parameters 

# with variable selection , and with data (add.data=TRUE)
dummy<-graphic.ppls.splines(X,y,lambda=100,ncomp=5,add.data=TRUE,select=TRUE,window.size=c(3,4))

# without variable selection and without data

dev.new()
dummy<-graphic.ppls.splines(X,y,lambda=100,ncomp=5,add.data=FALSE,select=FALSE,window.size=c(3,4))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
