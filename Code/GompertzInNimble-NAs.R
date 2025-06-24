lib.load <- function(){
  library(nimble)
  library(dclone)
  library(bbmle)
  library("MASS")
}

# load
lib.load()


# 1 sp sim Gompertz:
Pois.gssm <- function(pars=c(mu,theta,sigsq),len=30,T.out=TRUE){
  
  cc <- exp(-theta);
  a  <- mu*(1-cc);
  Vinf <- sigsq/(1-exp(-2*theta));
  
  X <- rep(0,len)
  N <- rep(0,len)
  X[1] <- rnorm(n=1, mean=mu, sd=sqrt(Vinf))
  N[1] <- rpois(n=1, lambda=exp(X[1]))
  for(i in 2:len){
    X[i] <- a+ cc*X[(i-1)] + rnorm(n=1, mean=0, sd=sqrt(sigsq))
    N[i] <-  rpois(n=1,lambda=exp(X[i]))
  }
  
  out <- cbind(X,N)
  colnames(out) <- c("X", "N")
  
  if(T.out==TRUE){
      out <- cbind(0:(len-1),X,N)
      colnames(out) <- c("Tvec", "X", "N")
  }
  
  return(out)  
}




##### --------------------------------------------------------------------------
##### Nimble Model and functions
cumsumNimble <- nimbleFunction(
  run = function(x = double(1)) {
    n <- length(x)
    result <- numeric(n)
    result[1] <- x[1]
    for (i in 2:n) {
      result[i] <- result[i-1] + x[i]
    }
    return(result)
    returnType(double(1))
  }
)
cumsumNimble <- compileNimble(cumsumNimble)


GSSM <- nimbleCode({
  
  ##### ------------------------------------------------------------------------
  ##### Priors
  lmu    ~ dnorm(mean=6, sd=3);
  ltheta ~ dnorm(mean=0.5, sd=0.01);
  lsigmasq ~ dnorm(mean=0, sd=1);
  
  mu <- exp(lmu);
  theta <- exp(ltheta);
  sigmasq <- exp(lsigmasq);
  
  
  #tausq ~ dchisq(df = 1);
  
  car.cap <- exp(mu);
  cc <- exp(-theta);
  a  <- mu*(1-cc);
  Vinf <- sigmasq/(1-exp(-2*theta));
  
  ##### ----------------------------------------- -------------------------------
  ##### Likelihood
  for(k in 1:ncl){
        
      X[1,k] ~ dnorm(mean=mu, sd=sqrt(Vinf))
      N[1,k] ~ dpois(lambda=exp(X[1,k]))
      for( i in 2:qp1){ 
        X[i, k]  ~ dnorm( mean = (a +cc*X[(i-1),k]), sd=sqrt(sigmasq) ); ## Process error
        N[i, k]  ~ dpois(lambda = exp(X[i,k])) ## Observation error 
      }
  }
})


GSSM.kalman <- nimbleCode({
  
  ##### ----------------------------------------- -------------------------------
  ##### Likelihood
  X[1] ~ dnorm(mean=mu, sd=sqrt(Vinf))
  Y[1] ~ dpois(lambda=exp(X[1]))
  Nhid[1] <- exp(X[1])
    for( i in 2:qp1){
      X[i]  ~ dnorm( mean = (a +cc*X[(i-1)]), sd=sqrt(sigmasq) ); ## Process error
      N[i]  ~ dpois(lambda = exp(X[i])) ## Observation error 
      Nhid[i] <- exp(X[i])
    }
})



guess.calc <- function(Yobs,Tvec){
  
  T.t <-Tvec-Tvec[1]; #  For calculations, time starts at zero.
  q <- length(Yobs)-1;      #  Number of time series transitions, q.
  qp1 <- q+1;              #  q+1 gets used a lot, too.
  S.t <- T.t[2:qp1]-T.t[1:q];  #  Time intervals.
  Ybar <- mean(Yobs);
  Yvar <- sum((Yobs-Ybar)*(Yobs-Ybar))/q;
  mu1 <- Ybar;
  
  # Kludge an initial value for theta based on mean of Y(t+s) given Y(t).
  th1<- -mean(log(abs((Yobs[2:qp1]-mu1)/(Yobs[1:q]-mu1)))/S.t);            
  bsq1<- 2*th1*Yvar/(1+2*th1);         # Moment estimate using stationary
  tsq1<- bsq1;                         #   variance, with betasq=tausq.
  
  #three 0's 
  three0s <- sum(c(th1,bsq1,tsq1))
  if(three0s==0|is.na(three0s)){th1 <- 0.5;bsq1 <- 0.09; tsq1 <- 0.23;}
  
  
  out1 <- c(th1,bsq1,tsq1);
  if(sum(out1<1e-7)>=1){out1 <- c(0.5,0.09,0.23)}
  out <- c(mu1,out1);
  return(abs(out))
  
}

#Let's wrap this initial guess estimator function into another function 
#that takes into account possible NAs

guess.calc2.0<- function(TimeAndNs){
  
  newmat <- TimeAndNs # to be replaced if is.na ==TRUE
  isnas <- sum(is.na(TimeAndNs))
  
  if(isnas >= 1){
    
    isnaind <- which(is.na(TimeAndNs[,2]), arr.ind=TRUE)
    newmat <- TimeAndNs[-isnaind,]
    newmat[,1] <- newmat[,1] - newmat[1,1]
    
  }
  
  init.guess <- guess.calc(Yobs = log(newmat[,2]), Tvec=newmat[,1])
  
  mu1  <- init.guess[1]
  th1  <- init.guess[2]
  bsq1 <- init.guess[3]
  sigsq1<- ((1-exp(-2*th1))*bsq1)/(2*th1)
  
  out <- c(mu=mu1, theta=th1, sigmasq = sigsq1)
  return(out)
}



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

####### What follows is code for NIMBLE #########
kroneckerProd <- nimbleFunction(
  run = function(A = double(2), B = double(2)) {
    # Get dimensions
    p <- dim(A)[1]
    # Output will be (p*p) x (p*p)
    out <- matrix(0, nrow = p*p, ncol = p*p)
    for(i in 1:p) {
      for(j in 1:p) {
        for(k in 1:p) {
          for(l in 1:p) {
            # Row index in output
            rowIdx <- (i - 1) * p + k
            # Col index in output
            colIdx <- (j - 1) * p + l
            out[rowIdx, colIdx] <- A[i, j] * B[k, l]
          }
        }
      }
    }
    returnType(double(2))
    return(out)
  }
)
compileNimble(kroneckerProd)


GSSM_MV <- nimbleCode({
  ##### Priors
  lmuvec[1:p]    ~ dmnorm(mupmean[1:p], mupvar[1:p,1:p]);
  lthetavec[1:p] ~ dmnorm(lthetamean[1:p], lthetavar[1:p,1:p]);
  lsigmasq ~ dnorm(mean=0, sd=1);
  mig ~ dbeta(2,2)

  muvec[1:p] <- exp(lmuvec[1:p]);
  thetavec[1:p] <- exp(lthetavec[1:p]);
  sigmasq <- exp(lsigmasq);
  
  cvec[1:p] <- exp(-thetavec[1:p]);
  preB[1:p,1:p] <- dinvmat[1:p,1:p]*mig
  B[1:p,1:p] <- put.diag(preB[1:p,1:p],cvec[1:p])   #3
  A[1:p]  <- (create.I(p) - B[1:p,1:p])%*%muvec[1:p] #4

  Sigma[1:p,1:p] <- sigmasq*create.I(p) #5

  ##### ----------------------------------------- -------------------------------
  ##### Likelihood
    Xarr[1:p,1] ~ dmnorm(muvec[1:p], Sigma[1:p,1:p])#dmnorm(muvec[1:p],V[1:p,1:p])
    for(h in 1:p){
      Narr[h,1] ~ dpois(exp(Xarr[h,1]))
    }
    for( i in 2:qp1){ 
      Xarr[1:p,i]  ~ dmnorm(Xarr[1:p,i-1], Sigma[1:p,1:p]); ## Process error
      for(j in 1:p){
        Narr[j,i]  ~ dpois(exp(Xarr[j,i])) ## Observation error
      }
    }
}
)

get.diag<- nimbleFunction(
  run=function(A=double(2)){
    
    out <- diag(A)
    returnType(double(1))
    return(out)
  }
)
compileNimble(get.diag)
#A <- matrix(1:9, 3, 3)
#get.diag(A)

put.diag <- nimbleFunction(
  
  run=function(A=double(2),d=double(1)){
    
    diag(A) <- d
    newA <- A
    returnType(double(2))
    return(newA)
  }
)
compileNimble(put.diag)
#A <- matrix(1:9, 3, 3)
#put.diag(A, c(0,0,0))


create.I <- nimbleFunction(
  run=function(p=double(1)){
    
    out <- diag(p)
    returnType(double(2))
    return(out)
  }
)
compileNimble(create.I)
#create.I(3)

asVector_nimble <- nimbleFunction(
  run = function(mat = double(2)) {
    nrow <- dim(mat)[1]
    ncol <- dim(mat)[2]
    out <- numeric(nrow * ncol, init = FALSE)
    idx <- 1
    for(j in 1:ncol) {
      for(i in 1:nrow) {
        out[idx] <- mat[i, j]
        idx <- idx + 1
      }
    }
    returnType(double(1))
    return(out)
  }
)

compileNimble(asVector_nimble)
A <- matrix(1:4, 2, 2)
as.vector(A)
asVector_nimble(A)


sliceArray <- nimbleFunction(
  run = function(arr = double(3), col = integer(0), slice = integer(0)) {
    p <- dim(arr)[1]
    out <- numeric(p, init = FALSE)
    for(i in 1:p) {
      out[i] <- arr[i, col, slice]
    }
    returnType(double(1))
    return(out)
  }
)
compileNimble(sliceArray)
sliceArray(data.mv$Narr, col=2, slice=1)





# A <- matrix(1:4, 2, 2)
# B <- matrix(5:8, 2, 2)
# kroneckerProd(A,B)
# kronecker(A,B)


### Old code:

# GSSM_MV <- nimbleCode({
#   
#   # data should be: 
#   # countsmat = matrix of counts, rows being demes, columns being time steps  
#   # p = scalar, number of demes
#   # qp1 = number of time steps, including first. qp1 == nrow(countsmat)
#   # ncl = number of clones
#   # Tvec = vector of time steps, starting at 0, then 1,2,3,.. 
#   
#   ##### Priors
#   mupmean[1:p] <- rep(6,p)
#   mupvar[1:p,1:p] <- p*create.I(p) #1
#   lmuvec[1:p]    ~ dmnorm(mupmean[1:p], mupvar[1:p,1:p]);
#   
#   lthetamean[1:p] <- rep(0.5,p)
#   lthetavar[1:p,1:p]  <- 0.1*create.I(p) #2
#   lthetavec[1:p] ~ dmnorm(lthetamean[1:p], lthetavar[1:p,1:p]);
#   lsigmasq ~ dnorm(mean=0, sd=1);
#   mig ~ dbeta(2,2)
#   
#   
#   muvec[1:p] <- exp(lmuvec[1:p]);
#   thetavec[1:p] <- exp(lthetavec[1:p]);
#   sigmasq <- exp(lsigmasq);
#   
#   cvec[1:p] <- exp(-thetavec[1:p]);
#   preB[1:p,1:p] <- dinvmat[1:p,1:p]*mig
#   B[1:p,1:p] <- put.diag(preB[1:p,1:p],cvec[1:p])   #3
#   A[1:p]  <- (create.I(p) - B[1:p,1:p])%*%muvec[1:p] #4
#   
#   Sigma[1:p,1:p] <- sigmasq*create.I(p) #5
#   #myI[1:psq,1:psq] <- create.I(psq)
#   #kronB[1:psq,1:psq] <- kroneckerProd(B[1:p,1:p],B[1:p,1:p])
#   #ImkronBInv[1:psq,1:psq] <- inverse(myI[1:psq,1:psq] - kronB[1:psq,1:psq])
#   #Vec.V[1:psq] <- ImkronBInv[1:psq,1:psq]%*%asVector_nimble(Sigma[1:p,1:p]) #6 # Eq. 17th Ives et al 2003
#   #Vpre[1:p,1:p] <- nimMatrix(Vec.V[1:psq], nrow=p,ncol=p)
#   #Vdiag1[1:p] <- get.diag(Vpre[1:p,1:p])
#   #Vdiag2[1:p] <- Vdiag11[1:p] + rep(1e-6,p)
#   #V[1:p,1:p] <- put.diag(Vpre[1:p,1:p],Vdiag2[1:p])
#   
#   ##### ----------------------------------------- -------------------------------
#   ##### Likelihood
#   for(k in 1:ncl){  
#     Xarr[1:p,1,k] ~ dmnorm(muvec[1:p], Sigma[1:p,1:p])#dmnorm(muvec[1:p],V[1:p,1:p])
#     for(h in 1:p){
#       Narr[h,1,k] ~ dpois(exp(Xarr[h,1,k]))
#     }
#     for( i in 2:qp1){ 
#       
#       meanvec[1:p,i,k] <- Xarr[1:p,i-1,k] #A[1:p] #+ Xarr[1:p,i-1,k]# sliceArray(Xarr[1:p,1:qp1,1:ncl], col=im1s[i], slice=k) + B[1:p,1:p]%*%sliceArray(Xarr[1:p,1:qp1,1:ncl], col=im1s[i], slice=k)
#       Xarr[1:p,i, k]  ~ dmnorm(meanvec[1:p,i,k], Sigma[1:p,1:p]); ## Process error
#       for(j in 1:p){
#         Narr[j,i, k]  ~ dpois(exp(Xarr[j,i,k])) ## Observation error
#       }
#     }
#   }
# }
# )