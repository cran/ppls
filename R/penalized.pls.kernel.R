`penalized.pls.kernel` <-
function(X,y,M=NULL,ncomp){
  n=nrow(X)
  yhat=rep(0,n)
  if (is.null(M)==TRUE){
	K=X%*%t(X)
  }
  if (is.null(M)==FALSE){
  	K=X%*%M%*%t(X)
  }
  UU=matrix(,n,ncomp)
  TT=matrix(,n,ncomp)
  for (i in 1:ncomp){
        uu=y-yhat
        uu=uu/sqrt(sum((K%*%uu)*uu))
        UU[,i]=uu
        if (i==1) tt=K%*%uu else {
            TTi=TT[,1:(i-1),drop=FALSE]
            Kuu=K%*%uu
            tt=Kuu - TTi %*%(t(TT)%*%Kuu)
            # Theoretically, we do not need the projection TT%*%t(TT) on ALL of the components but
            # only the projection tt%*%t(tt) on the last component (see the equation above algorithm
            # 3 in our paper) 
            # but to ensure pairwise orthogonality, we use the projection TT%*%t(TT)
            if (floor(i/5)==i/5){
                tt=tt-TTi%*%(t(TTi)%*%tt) # after a while, the components are not orthogonal anymore ....
            }
        }
        tt=normalize.vector(tt)
        TT[,i]=tt
        yhat=yhat + sum(tt*y) * tt
    }
    RR=t(TT)%*%K%*%UU
    RR[row(RR)>col(RR)]<-0 # inserted for numerical stability
    RR[row(RR)< col(RR) -1]<-0 # inserted for numerical stability
    LL=backsolve(RR,diag(ncomp))
    AA=matrix(,nrow(X),ncomp)
    for (i in 1:ncomp){
        Li=LL[1:i,1:i,drop=FALSE]
        AA[,i]= UU[,1:i,drop=FALSE]%*%(Li%*%(t(TT[,1:i,drop=FALSE])%*%y))
    }
    coefficients=t(X)%*%AA
	if (is.null(M)==FALSE){
	coefficients=M%*%coefficients	
}
  return(list(coefficients=coefficients))

}
