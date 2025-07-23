######################. FUNCTIONS AND LIBRARIES TO RUN THE CODE #############


lib.load <- function(){
  library(nimble)
  library(dclone)
  library(bbmle)
  library("MASS")
}

# load
lib.load()


# 1 sp sim Gompertz:
Pois.gssm <- function(pars=c(mu,theta,sigsq),len=30, pvec, T.out=TRUE){
  
  cc <- exp(-theta);
  a  <- mu*(1-cc);
  Vinf <- sigsq/(1-exp(-2*theta));
  
  X <- rep(0,len)
  N <- rep(0,len)
  X[1] <- rnorm(n=1, mean=mu, sd=sqrt(Vinf))
  N[1] <- rpois(n=1, lambda=pvec[1]*exp(X[1]))
  for(i in 2:len){
    X[i] <- a+ cc*X[(i-1)] + rnorm(n=1, mean=0, sd=sqrt(sigsq))
    N[i] <-  rpois(n=1,lambda=pvec[i]*exp(X[i]))
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
    N[1,k] ~ dpois(lambda=pvec[1]*exp(X[1,k]))
    for( i in 2:qp1){ 
      X[i, k]  ~ dnorm( mean = (a +cc*X[(i-1),k]), sd=sqrt(sigmasq) ); ## Process error
      N[i, k]  ~ dpois(lambda = pvec[i]*exp(X[i,k])) ## Observation error 
    }
  }
})


GSSM.kalman <- nimbleCode({
  
  ##### ----------------------------------------- -------------------------------
  ##### Likelihood
  X[1] ~ dnorm(mean=mu, sd=sqrt(Vinf))
  Y[1] ~ dpois(lambda=pvec[1]*exp(X[1]))
  Nhid[1] <- exp(X[1])
  for( i in 2:qp1){
    X[i]  ~ dnorm( mean = (a +cc*X[(i-1)]), sd=sqrt(sigmasq) ); ## Process error
    N[i]  ~ dpois(lambda = pvec[i]*exp(X[i])) ## Observation error 
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




############# PROCEDURAL SECTION ##################################

mu <- 5
theta <- 0.23
sigsq <- 0.09

cc <- exp(-theta)
a <- mu*(1-cc)
Vinf <- sigsq/(1-cc^2)

len <- 30
pvec <- rbeta(n=len, shape1=3,shape2=1)
Onesim <- Pois.gssm(pars=c(mu=mu, theta=theta, sigsq=sigsq), pvec=pvec,len=len)
plot(1:len, Onesim[,3], pch=16, type="b", xlab="Time", ylab="Population size",
     ylim=c(0,max(exp(Onesim[,2]))),bty="l")
points(1:len,exp(Onesim[,2]), pch=1, type="b")

# Data Cloning in Nimble
ncl <- 30
data <- list(N=nimMatrix(Onesim[,3], nrow= nrow(Onesim), ncol=ncl))
constants <- list(qp1 = nrow(Onesim), ncl=ncl, pvec=pvec)
monitors <- c("mu", "theta", "sigmasq", "a", "cc", "car.cap", "Vinf")
inits <- list(lmu = log(4), ltheta=log(0.5), lsigmasq=log(0.1))

tmp.model <- nimbleMCMC( code=GSSM, data = data, 
                         constants=constants, inits = inits,
                         monitors = monitors, thin = 10,
                         niter = 50000, nburnin = 10000, nchains = 3,
                         summary = TRUE, samplesAsCodaMCMC = TRUE)

cbind(tmp.model$summary$all.chains[,1], c(
  Vinf,a,car.cap=exp(mu),cc,mu,sigsq, theta))

# Kalman Estimates and plot:

data <- list(N=Onesim[,3])
constants <- append( as.list(tmp.model$summary$all.chains[,1]),  constants)

gss.kalman <- nimbleMCMC( code=GSSM.kalman, data = data, 
                          constants=constants, 
                          monitors = c("Nhid"), thin = 10,
                          niter = 50000, nburnin = 10000, nchains = 3,
                          summary = TRUE,
                          dimensions = list(Nhid = length(Onesim[,2])))

N.hat <-  gss.kalman$summary$all.chains[,1]
N.low <-  gss.kalman$summary$all.chains[,4]
N.hi <-  gss.kalman$summary$all.chains[,5]

par(oma=c(1,1,2,1), mar=c(4,5,1,1))
plot(1:len, exp(Onesim[,2]), pch=16, type="b", xlab="Time", ylab="Population size", bty="l", ylim=c(0,350), col="blue", cex.lab=1.15)
polygon(c(1:len,rev(1:len)), c(N.low,rev(N.hi)), col="pink", border=NA)
points(1:len, exp(Onesim[,2]), pch=16, type="b", col="blue")
points(1:len, N.hat, pch=16, type="b", col="red")
points(1:len, data$N, pch=16, type="b", col="black")
legend(0,373, legend=c("True abundance", "Estimated abundance", "Observed counts"),
       col=c("blue", "red","black"), pch=c(16,16,16), lty=c(1,1,1), cex=1.15, bty="n", horiz=TRUE)

