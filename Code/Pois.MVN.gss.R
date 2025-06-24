# Random Poisson-Multivariate SSM data generator not including NA generation

Pois.MVN.gss <- function(pars=list(mu=rep(5,3),theta=rep(0.23,3),sigsq=0.09,mig=0.15,
                                   dmat=matrix(c(0,5,10,5,0,5,10,5,0), nrow=3,ncol=3)),
                         len=30){
  
  mu <- pars$mu
  p <- length(mu)
  theta <- pars$theta
  sigsq <- pars$sigsq
  mig <- pars$mig
  dmat <- pars$dmat
  
  
  cvec <- exp(-theta)  
  B <- (dmat/dmat^2)*mig
  diag(B) <- cvec
  
  A  <- (diag(p) - B)%*%mu
  
  Sigma <- diag(rep(sigsq,p))
  Vec.V <- ginv(diag(p*p) - kronecker(B,B))%*%as.vector(Sigma) # Eq. 17th Ives et al 2003
  V <- matrix(Vec.V, nrow=p,ncol=p,byrow=FALSE)
  
  Xmat <- matrix(0,ncol=len,nrow=p)
  Nmat <- Xmat	
  Xo <- randmvn(n=1,mu.vec= ginv(diag(p)-B)%*%A, cov.mat=V)
  Xmat[,1] <- Xo
  Nmat[,1] <- rpois(n=p, lambda=exp(Xo))
  
  for( i in 2:len){
    
    im1 <- i-1;
    Xim1  <- matrix(Xmat[,im1], nrow=3,ncol=1)
    mui <- A + B%*%Xim1
    rand.trans <- randmvn(n=1,mu.vec= mui, cov.mat=Sigma)
    Xmat[,i] <- rand.trans
    Nmat[,i] <- rpois(n=3,lambda=exp(rand.trans))
  }
  
  return(list(Xmat=Xmat, Nmat=Nmat, tvec=0:(len-1)))	
}
