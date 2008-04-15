`penalized.pls` <-
function(X,y,P=NULL,ncomp=NULL,kernel=FALSE,scale=FALSE,blocks=1:ncol(X),select=FALSE){

  p<-ncol(X)

  y<-as.vector(y)
  
  if (is.null(P)) P=matrix(0,p,p)

  meanx=apply(X,2,mean)

  meany=mean(y)
  if (scale==TRUE) {
    sdx=sqrt(apply(X,2,var))
    }
    else {
    sdx=rep(1,ncol(X))
  }

  if (is.null(ncomp)) ncomp=min(p,nrow(X)-1)

  X<-scale(X,center=TRUE,scale=scale)

  y<-scale(y,center=TRUE,scale=FALSE)

  Minv<-diag(p)+P

  M<-solve(Minv)

  if (select==FALSE){  

	if (kernel==TRUE) ppls=penalized.pls.kernel(X,y,M,ncomp); 
	if (kernel==FALSE)  ppls=penalized.pls.default(X,y,M,ncomp);
}
if (select==TRUE) ppls=penalized.pls.select(X,y,M,ncomp,blocks)
  
  
  coefficients=ppls$coefficients /(sdx %*%t(rep(1,ncol(ppls$coefficients))))
  
  
   intercept=rep(meany,ncomp) - t(coefficients)%*%meanx
  
  return(list(intercept=intercept,coefficients=coefficients))

}
