`penalized.pls.default` <-
function(X,y,M,ncomp=min(nrow(X)-1,ncol(X))){

  X0<-X

  WW<-c()

  TT<-c()
  
  length.tt=vector(length=ncomp)

  for (i in 1:ncomp){

      ww<-M%*%t(X)%*%y
      # ww=normalize.vector(ww)
      tt<-X%*%ww
      length.tt[i]=sqrt(sum(tt^2))
      tt=normalize.vector(tt)
      WW=cbind(WW,ww)
      TT=cbind(TT,tt)
      X=X- tt%*%t(tt)%*%X
    }
    B=matrix(,ncol(X),ncomp)
        RR=matrix(t(TT)%*%X0%*%WW,nrow=ncomp) # the bidiagonal matrix
        RR[row(RR)>col(RR)]<-0 # inserted for numerical stability
        RR[row(RR)< col(RR) -1]<-0 # inserted for numerical stability
        LL=backsolve(RR,diag(ncomp))
      Wtilde=WW%*%LL
        B=matrix(,ncol(X),ncomp)
        for (i in 1:ncomp) {
            Li=LL[1:i,1:i,drop=FALSE]
        B[,i]=WW[,1:i,drop=FALSE]%*%Li%*%t(TT[,1:i,drop=FALSE])%*%y
    }
  return(list(coefficients=B,Wtilde=Wtilde))

}
