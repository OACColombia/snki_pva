# source the library with the main functions:
source("GompertzinNimble-NAs.R")

library(dclone); library(MASS); library(mcmcplots)

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
  #Xo <- randmvn(n=1,mu.vec= ginv(diag(p)-B)%*%A, cov.mat=V)
  Xo <- randmvn(n=1,mu.vec= ginv(diag(p)-B)%*%A, cov.mat=Sigma)
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



trial.sim <- Pois.MVN.gss(pars=list(mu=rep(5,3),theta=rep(0.23,3),sigsq=0.09,mig=0.15,
                                    dmat=matrix(c(0,5,10,5,0,5,10,5,0), nrow=3,ncol=3)),
                          len=30)
#true.c <- exp(-0.23)
#true.c

# Bits and pieces (constants, mostly) to run the MCMC
p <- nrow(trial.sim$Nmat)
dmat <- matrix(c(0,5,10,5,0,5,10,5,0), nrow=p,ncol=p)
dinvmat <- (dmat/dmat^2)  
diag(dinvmat) <- rep(0,p)
onets4guess <- cbind(trial.sim$tvec, log(trial.sim$Nmat[1,]))
naive.guess <- guess.calc2.0(TimeAndNs = onets4guess)
myinits <- list(lmuvec = log(rep(naive.guess[1],p)), lthetavec=log(rep(naive.guess[2],p)), 
              lsigmasq=log(naive.guess[3]), mig=0.20)

MGSS_JAGS <- function(){
  ##### Priors
  lmuvec    ~ dmnorm(mupmean, muprec)
  lthetavec ~ dmnorm(lthetamean, lthetaprec)
  lsigmasq ~ dnorm(0, 1)
  mig ~ dbeta(2,2)

  sigmasq <- exp(lsigmasq)  
  for(pp in 1:p){
    muvec[pp] <- exp(lmuvec[pp])
    thetavec[pp] <- exp(lthetavec[pp])
    cvec[pp] <- exp(-thetavec[pp])
  }

  for(m in 1:p){
    B[m,m] <- cvec[m]
    for(n in 1:(m-1)){
      B[m,n] <- dinvmat[m,n]*mig
    }
    for(n in (m+1):p){
      B[m,n] <- dinvmat[m,n]*mig
    }
  }
  
  
  A  <- (Ip - B)%*%muvec
  
  Sigma <- (1/sigmasq)*Ip

  ##### ----------------------------------------- -------------------------------
  ##### Likelihood
  X[1:p,1] ~ dmnorm(muvec, Sigma)
  for(h in 1:p){
    N[h,1] ~ dpois(exp(X[h,1]))
  }
  for( i in 2:qp1){ 
    X[1:p,i]  ~ dmnorm(A + B%*%X[1:p,(i-1)], Sigma); ## Process error
    for(j in 1:p){
      N[j,i]  ~ dpois(exp(X[j,i])) ## Observation error
    }
  }
}



data4jags <- list(N=trial.sim$Nmat, qp1=ncol(trial.sim$Nmat), p=p,
                  mupmean=rep(6,p), muprec=(1/3)*diag(p),
                  lthetamean = rep(0.5,p),
                  lthetaprec =  10*diag(p),
                  Ip =  diag(p),
                  dinvmat = dinvmat)




n.iter<-100000
n.adapt<-1000
n.update<-100
thin<-10
n.chains<-3

out.parms <- c("muvec", "cvec", "mig", "sigmasq", "A")
mgssjags <- jags.fit(data4jags, params=out.parms, model=MGSS_JAGS, n.chains = n.chains, 
                    n.adapt=n.adapt, n.update=n.update, n.iter = n.iter, thin=thin,
                    inits=myinits) 

all.chains <- mgssjags[[1]]
for(i in 2:n.chains){all.chains <- rbind(all.chains,mgssjags[[i]])}
dim(all.chains)
posterior.quantiles <- quantile.mcmc.list(mgssjags)

################ Next step: Data Cloning #####################

MGSS_DC <- function(){
  ##### Priors
  lmuvec    ~ dmnorm(mupmean, muprec)
  lthetavec ~ dmnorm(lthetamean, lthetaprec)
  lsigmasq ~ dnorm(0, 1)
  mig ~ dbeta(2,2)
  
  sigmasq <- exp(lsigmasq)  
  for(pp in 1:p){
    muvec[pp] <- exp(lmuvec[pp])
    thetavec[pp] <- exp(lthetavec[pp])
    cvec[pp] <- exp(-thetavec[pp])
  }
  
  Sigma <- (1/sigmasq)*Ip
    

  for(m in 1:p){
    B[m,m] <- cvec[m]
    for(n in 1:(m-1)){
      B[m,n] <- dinvmat[m,n]*mig
    }
    for(n in (m+1):p){
      B[m,n] <- dinvmat[m,n]*mig
    }
  }

  A[1:p]  <- (Ip - B[1:p,1:p])%*%muvec
  
  for(k in 1:K){
      
    ##### ------------------------------------------------------------------
    ##### Likelihood
    X[1:p,1,k] ~ dmnorm(muvec, Sigma)
    for(h in 1:p){
      N[h,1,k] ~ dpois(exp(X[h,1,k]))
    }
    for( i in 2:qp1){ 
      X[1:p,i,k]  ~ dmnorm(A[1:p] + B[1:p,1:p]%*%X[1:p,(i-1),k], Sigma); ## Process error
      for(j in 1:p){
        N[j,i,k]  ~ dpois(exp(X[j,i,k])) ## Observation error
      }
    }
    
  }  
}





############### Using dclone functions to speed the process ##########


data4dclone <- list(K=1, N=dcdim(array(trial.sim$Nmat,dim=c(p,ncol(trial.sim$Nmat),1))), 
                qp1=ncol(trial.sim$Nmat), p=p, mupmean=rep(6,p), muprec=(1/3)*diag(p),
                lthetamean = rep(0.5,p), lthetaprec =  10*diag(p),
                Ip =  diag(p),dinvmat = dinvmat)

cl.seq <- c(1,4,8,10);
n.iter<-100000;n.adapt<-50000;n.update<-100;thin<-10;n.chains<-3;

out.parms <- c("muvec", "cvec", "mig", "sigmasq", "A")
mgssdclone <- dc.fit(data4dclone, params=out.parms, model=MGSS_DC, n.clones=cl.seq,
                     multiply="K", unchanged = c("qp1","p","mupmean","muprec",
                                                 "lthetamean","lthetaprec",
                                                 "Ip","dinvmat"),
                     n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                     n.iter = n.iter, thin=thin, inits=myinits) 
dcdiag(mgssdclone)

############### Using straight up JAGS code #########
K <- 10
Narray <- array(NA,dim=c(p,ncol(trial.sim$Nmat),K))
for(i in 1:K){Narray[,,i] <- trial.sim$Nmat}

data4dc <- list(N=Narray,K=K,qp1=ncol(trial.sim$Nmat), p=p,
                  mupmean=rep(6,p), muprec=(1/3)*diag(p),
                  lthetamean = rep(0.5,p),
                  lthetaprec =  10*diag(p),
                  Ip =  diag(p),
                  dinvmat = dinvmat)


n.iter<-100000;n.adapt<-1000;n.update<-1000;thin<-10;n.chains<-3;
out.parms <- c("muvec", "cvec", "mig", "sigmasq", "A")
mgssjags <- jags.fit(data4dc, params=out.parms, model=MGSS_DC, n.chains = n.chains, 
                     n.adapt=n.adapt, n.update=n.update, n.iter = n.iter, thin=thin,
                     inits=myinits) 

#for Data Cloning:

full.chain <- do.call(rbind, mgssjags)
vcovmat <- var(full.chain)
lambdas.vcov <- eigen(vcovmat)
lambda.max <- max(lambdas.vcov$values)
lambda.max

Fish.Inv <- (1/K)*vcovmat
Wald.halfcis <- 1.96*sqrt(diag(Fish.Inv))

MLES <- apply(full.chain, 2, mean) 


##### Kalman estimates #########

MGSS_Kalman <- function(){

  Sigma <- (1/sigmasq)*Ip
  
  
  for(m in 1:p){
    B[m,m] <- cvec[m]
    for(n in 1:(m-1)){
      B[m,n] <- dinvmat[m,n]*mig
    }
    for(n in (m+1):p){
      B[m,n] <- dinvmat[m,n]*mig
    }
  }
  
  #A[1:p]  <- (Ip - B[1:p,1:p])%*%muvec
  
    ##### ------------------------------------------------------------------
    ##### Likelihood
    X[1:p,1] ~ dmnorm(muvec, Sigma)
    for(h in 1:p){
      N[h,1] ~ dpois(exp(X[h,1]))
    }
    for( i in 2:qp1){ 
      X[1:p,i]  ~ dmnorm(A[1:p] + B[1:p,1:p]%*%X[1:p,(i-1)], Sigma); ## Process error
      for(j in 1:p){
        N[j,i]  ~ dpois(exp(X[j,i])) ## Observation error
      }
    }
    
}


data4kalman <- list(N=Narray[,,1],qp1=ncol(trial.sim$Nmat), p=p,
                Ip =  diag(p), dinvmat = dinvmat,
                A = MLES[grep("A",names(MLES))], cvec= MLES[grep("cvec",names(MLES))],
                mig=MLES[grep("mig",names(MLES))], muvec=MLES[grep("muvec", names(MLES))],
                sigmasq=MLES[grep("sigmasq",names(MLES))])


n.iter<-100000;n.adapt<-1000;n.update<-1000;thin<-10;n.chains<-3;

mgsskalman <- jags.fit(data4kalman, params=c("X"), model=MGSS_Kalman, n.chains = n.chains, 
                     n.adapt=n.adapt, n.update=n.update, n.iter = n.iter, thin=thin) 

Xkalman <-  do.call(rbind, mgsskalman)

Kalman.quants <- apply(Xkalman, 2, FUN=function(x){quantile(x, probs=c(0.025,0.5,0.975))}) 

len <- ncol(trial.sim$Nmat)
Nmat <- trial.sim$Nmat
Xmat <- trial.sim$Xmat

names.quants <- substr(colnames(Kalman.quants), start=1, stop=3)
X1s <- Kalman.quants[,which(names.quants=="X[1", arr.ind=TRUE)]
X2s <- Kalman.quants[,which(names.quants=="X[2", arr.ind=TRUE)]
X3s <- Kalman.quants[,which(names.quants=="X[3", arr.ind=TRUE)]



par(oma=c(1,1,1,1), mar=c(4,5,2,1))
plot(1:len, Xmat[1,], pch=16, type="b", col="red", ylim=c(4,max(Xmat)), main ="3-Demes log-abundances", 
     xlab="Time", ylab="3-Demes true vs. estimated log-population sizes", bty="l")
polygon(c(1:len,rev(1:len)), c(X1s[1,],rev(X1s[3,])), col=scales::alpha("red",0.1), border=NA)
points(1:len, Xmat[1,], pch=16, type="b", col="red")
points(1:len, X1s[2,],pch=15, type="b", col="red")

polygon(c(1:len,rev(1:len)), c(X2s[1,],rev(X2s[3,])), col=scales::alpha("blue",0.1), border=NA)
points(1:len, Xmat[2,], pch=16, type="b", col="blue")
points(1:len, X2s[2,],pch=15, type="b", col="blue")

polygon(c(1:len,rev(1:len)), c(X3s[1,],rev(X3s[3,])), col=scales::alpha("darkgreen",0.1), border=NA)
points(1:len, Xmat[3,], pch=16, type="b", col="darkgreen")
points(1:len, X3s[2,],pch=15, type="b", col="darkgreen")





