#    rand.MVN:  Multivariate Normal random number generator
#    n = number of random samples of a MVN vector
#    mu = mean vector of the MVN distribution to sample from
#    cov.mat = Variance-covariance matrix of the MVN distribution to sample from
randmvn <- function(n,mu.vec, cov.mat){
  
  p <- length(mu.vec);
  Tau <- chol(cov.mat, pivot=TRUE);
  Zmat <- matrix(rnorm(n=p*n,mean=0,sd=1),nrow=p,ncol=n); #generate normal deviates outside loop
  out <- matrix(0,nrow=p,ncol=n);
  for(i in 1:n){
    
    Z <- Zmat[,i];
    out[,i] <- t(Tau)%*%Z + mu.vec
    
  }
  
  return(out)
  
}
