`penalized.pls.cv` <-
function(X,y,P=NULL,lambda=1,ncomp=NULL,k=5,kernel=FALSE,scale=FALSE){

  n<-nrow(X)
  
  if (is.null(ncomp)) ncomp=min(n-1,ncol(X))
  
  lambda=as.vector(lambda)
  if (is.null(P)==TRUE) P=matrix(0,ncol(X),ncol(X))
  error.cv=matrix(0,length(lambda),ncomp)
  all.folds <- split(1:n, rep(1:k,length=n))
  
  for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain=X[-omit,,drop=FALSE]
        ytrain=y[-omit]
        Xtest=X[omit,,drop=FALSE]
        ytest=y[omit]
        for (j in 1:length(lambda)){
            penpls=penalized.pls(Xtrain,ytrain,P=lambda[j]*P,ncomp,kernel=kernel,scale=scale)
            error.cv[j,]=error.cv[j,]+ length(ytest)*(new.penalized.pls(penpls,Xtest,ytest)$mse)
        }
        
    }

  error.cv=error.cv/n
  value1=apply(error.cv,1,min)
  lambda.opt=lambda[which.min(value1)]
  ncomp.opt=which.min(error.cv[lambda==lambda.opt,])
  min.ppls=min(value1)
  ppls=penalized.pls(X,y,P*lambda.opt,ncomp=ncomp.opt,kernel=kernel,scale=scale)
  intercept=ppls$intercept[ncomp.opt]
  coefficients=ppls$coefficients[,ncomp.opt]
  return(list(lambda.opt=lambda.opt,ncomp.opt=ncomp.opt,min.ppls=min.ppls,intercept=intercept,coefficients=coefficients))

}
