############## Loading the actual data as list of matrices and running 
############## JAGS on the full model 
library(dclone); library(MASS); library(mcmcplots)

source("guess.calc.R")
source("guess.calc2.0.R")
source("randmvn.R")
source("Pois.MVN.gss.R")

load("SNKIdata.RData") # from DataProcessing.R

S4mat <- Counts.list$S4mat[1:6,1:28] # Cut the last time step: no data anywhere
S5mat <- Counts.list$S5mat[1:6,1:28]
S6mat <- Counts.list$S6mat[1:6,1:28]

NA4mat <- NAS.list$NA4mat[1:6,1:28]
NA5mat <- NAS.list$NA5mat[1:6,1:28]
NA6mat <- NAS.list$NA6mat[1:6,1:28]

# Now we pick a single time series of counts, remove NAs, and "cbind it" to a 
# column of time, thus getting a matrix with time as a first column, counts 
# as a second column.  We use this matrix to compute a guess of the model params

ts.4guess  <- S4mat[3,]
tvec4guess  <- as.numeric(colnames(S4mat))
onets4guess <- cbind(tvec4guess, log(ts.4guess)) # No NAS for KRV
naive.guess <- guess.calc2.0(TimeAndNs = onets4guess)

# Bits and pieces (constants, mostly) to run the MCMC

# First get the distance matrix and re-order it according to the order in which
# I got the data (in terms of demes' order)
gen.dists <- readRDS("SNKI_dist_deems_general.rds")
min.dists <- readRDS("SNKI_dist_deems_minimal.rds")
demes.order0 <- colnames(gen.dists)

p <- nrow(S4mat)
demes <- rownames(S4mat) 
gen.dists1 <- matrix(0,nrow=p, ncol=p)
rownames(gen.dists1) <- demes
colnames(gen.dists1) <- demes
min.dists1 <- gen.dists1

for(i in 1:p){
  
  ith.d   <- demes[i]
  other.i <- which(demes.order0==ith.d,arr.ind=TRUE)
  
  for(j in 1:p){
    
    jth.d <- demes[j]
    other.j <- which(demes.order0==jth.d, arr.ind=TRUE)
    
    gen.dists1[i,j] <- gen.dists[other.i,other.j]
    min.dists1[i,j] <- min.dists[other.i,other.j]    
    
  }
  
}


myinits <- list(lmuvec = log(rep(naive.guess[1],p)), lthetavec=log(rep(naive.guess[2],p)), 
                lsigmasq=log(naive.guess[3]), mig=0.20)

dmat <- gen.dists1 # Try a generalized distance first
dinvmat <- (dmat/dmat^2)  
diag(dinvmat) <- rep(0,p)

# data4jags <- list(N=S4mat,
#                   qp1=ncol(S4mat), 
#                   p=p,
#                   mupmean=rep(6,p), 
#                   muprec=(1/3)*diag(p),
#                   lthetamean = rep(0.5,p),
#                   lthetaprec =  10*diag(p),
#                   Ip =  diag(p),
#                   dinvmat = dinvmat
# )


data4jags <- list(N4=S4mat,
                  N5=S5mat,
                  N6=S6mat,
                  qp1=ncol(S4mat), 
                  p=p,
                  mupmean=rep(6,p), 
                  muprec=(1/3)*diag(p),
                  lthetamean = rep(0.5,p),
                  lthetaprec =  10*diag(p),
                  Ip =  diag(p),
                  dinvmat = dinvmat
                  )




MGSS_JAGS2.0 <- function(){
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
    N4[h,1] ~ dpois(exp(X[h,1]))
    N5[h,1] ~ dpois(exp(X[h,1]))    
    N6[h,1] ~ dpois(exp(X[h,1]))    
  }
  for( i in 2:qp1){ 
    X[1:p,i]  ~ dmnorm(A + B%*%X[1:p,(i-1)], Sigma); ## Process error
    for(j in 1:p){
      N4[j,i]  ~ dpois(exp(X[j,i])) ## Observation error
      N5[j,i]  ~ dpois(exp(X[j,i])) ## Observation error      
      N6[j,i]  ~ dpois(exp(X[j,i])) ## Observation error
    }
  }
}

n.iter<-100000
n.adapt<-1000
n.update<-100
thin<-10
n.chains<-3

out.parms <- c("muvec", "cvec", "mig", "sigmasq", "A")
mgssjags <- jags.fit(data4jags, params=out.parms, model=MGSS_JAGS2.0, n.chains = n.chains, 
                     n.adapt=n.adapt, n.update=n.update, n.iter = n.iter, thin=thin,
                     inits=myinits) 

all.chains <- mgssjags[[1]]
for(i in 2:n.chains){all.chains <- rbind(all.chains,mgssjags[[i]])}
dim(all.chains)
posterior.quantiles <- quantile.mcmc.list(mgssjags)

################ Next step: Data Cloning #####################

MGSS_DC2.0 <- function(){
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
    for(n in 1:p){
      B[m,n] <- ifelse(m==n, cvec[m], dinvmat[m,n]*mig) #
    }
  }

  # for(m in 1:p){
  #   B[m,m] <- cvec[m]
  #   for(n in 1:(m-1)){
  #     B[m,n] <- dinvmat[m,n]*mig
  #   }
  #   for(n in (m+1):p){
  #     B[m,n] <- dinvmat[m,n]*mig
  #   }
  # }
  
  A[1:p]  <- (Ip - B[1:p,1:p])%*%muvec
  
  for(k in 1:K){
    
    ##### ------------------------------------------------------------------
    ##### Likelihood
    X[1:p,1,k] ~ dmnorm(muvec, Sigma)
    for(h in 1:p){
      N4[h,1,k] ~ dpois(exp(X[h,1,k]))
      N5[h,1,k] ~ dpois(exp(X[h,1,k]))      
      N6[h,1,k] ~ dpois(exp(X[h,1,k]))
    }
    for( i in 2:qp1){ 
      X[1:p,i,k]  ~ dmnorm(A[1:p] + B[1:p,1:p]%*%X[1:p,(i-1),k], Sigma); ## Process error
      for(j in 1:p){
        N4[j,i,k]  ~ dpois(exp(X[j,i,k])) ## Observation error
        N5[j,i,k]  ~ dpois(exp(X[j,i,k])) ## Observation error        
        N6[j,i,k]  ~ dpois(exp(X[j,i,k])) ## Observation error        
      }
    }
    
  }  
}

data4dclone <- list(K=1, 
                    N4=dcdim(array(S4mat,dim=c(p,ncol(S4mat),1))), 
                    N5=dcdim(array(S5mat,dim=c(p,ncol(S5mat),1))), 
                    N6=dcdim(array(S6mat,dim=c(p,ncol(S6mat),1))), 
                    qp1=ncol(S4mat), 
                    p=p, 
                    mupmean=rep(6,p), 
                    muprec=(1/3)*diag(p),
                    lthetamean = rep(0.5,p), 
                    lthetaprec =  10*diag(p),
                    Ip =  diag(p),
                    dinvmat = dinvmat
                    )

cl.seq <- c(1,4,8,10);
n.iter<-100000;n.adapt<-50000;n.update<-100;thin<-10;n.chains<-3;

out.parms <- c("muvec", "cvec", "mig", "sigmasq", "A")
mgssdclone <- dc.fit(data4dclone, params=out.parms, model=MGSS_DC2.0, n.clones=cl.seq,
                     multiply="K", unchanged = c("qp1","p","mupmean","muprec",
                                                 "lthetamean","lthetaprec",
                                                 "Ip","dinvmat"),
                     n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                     n.iter = n.iter, thin=thin, inits=myinits) 
dcdiag(mgssdclone)

summary(mgssdclone)

confint(mgssdclone)

mlesgendist <- summary(mgssdclone)[[1]]
confintgendist <- confint(mgssdclone)

# Migration proabilities MLES and Confidence Intervals
mig.row <- which(rownames(mlesgendist)=="mig", arr.ind=TRUE)

m.mle <- mlesgendist[mig.row,1]
m.var <- mlesgendist[mig.row,3]^2

var.migmat <- (dinvmat^2)*m.var
sd.migmat <- sqrt(var.migmat)

mle.migmat <- m.mle*dinvmat

low.ci.migmat <- mle.migmat - 1.96*sd.migmat
up.ci.migmat <- mle.migmat + 1.96*sd.migmat


##### Now re-running the estimation but using the minimal distance

dmat <- min.dists1 # Change the distance matrix!
dinvmat <- (dmat/dmat^2)  
diag(dinvmat) <- rep(0,p)

data4dclone <- list(K=1, 
                    N4=dcdim(array(S4mat,dim=c(p,ncol(S4mat),1))), 
                    N5=dcdim(array(S5mat,dim=c(p,ncol(S5mat),1))), 
                    N6=dcdim(array(S6mat,dim=c(p,ncol(S6mat),1))), 
                    qp1=ncol(S4mat), 
                    p=p, 
                    mupmean=rep(6,p), 
                    muprec=(1/3)*diag(p),
                    lthetamean = rep(0.5,p), 
                    lthetaprec =  10*diag(p),
                    Ip =  diag(p),
                    dinvmat = dinvmat
)

cl.seq <- c(1,4,8,10);
n.iter<-100000;n.adapt<-50000;n.update<-100;thin<-10;n.chains<-3;

out.parms <- c("muvec", "cvec", "mig", "sigmasq", "A")
# Btest <- dc.fit(data4dclone, params=out.parms, model=MGSS_DC2.0, n.clones=cl.seq,
#                       multiply="K", unchanged = c("qp1","p","mupmean","muprec",
#                                                   "lthetamean","lthetaprec",
#                                                   "Ip","dinvmat"),
#                       n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
#                       n.iter = n.iter, thin=thin, inits=myinits) 

mgssdclone2 <- dc.fit(data4dclone, params=out.parms, model=MGSS_DC2.0, n.clones=cl.seq,
                     multiply="K", unchanged = c("qp1","p","mupmean","muprec",
                                                 "lthetamean","lthetaprec",
                                                 "Ip","dinvmat"),
                     n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                     n.iter = n.iter, thin=thin, inits=myinits) 
dcdiag(mgssdclone2)

summary(mgssdclone2)

confint(mgssdclone2)

mlesmindist <- summary(mgssdclone2)[[1]]
confintmindist <- confint(mgssdclone2)

# Migration proabilities MLES and Confidence Intervals

mig.row <- which(rownames(mlesgendist)=="mig", arr.ind=TRUE)

m.mle2 <- mlesmindist[mig.row,1]
m.var2 <- mlesmindist[mig.row,3]^2

var.migmat2 <- (dinvmat^2)*m.var2
sd.migmat2 <- sqrt(var.migmat2)

mle.migmat2 <- m.mle2*dinvmat

low.ci.migmat2 <- mle.migmat2 - 1.96*sd.migmat2
up.ci.migmat2 <- mle.migmat2 + 1.96*sd.migmat2




save.image("SNKI-Results.RData")

################ Now, getting the Kalman Estimates ############################
################. Will PP counts be predicted to be LOW before 2019?????#######


MGSS_KALMAN2.0 <- function(){

  A  <- (Ip - B)%*%muvec
  
  Sigma <- (1/sigmasq)*Ip
  
  ##### ----------------------------------------- -------------------------------
  ##### Likelihood
  X[1:p,1] ~ dmnorm(muvec, Sigma)
  for(h in 1:p){
    N4[h,1] ~ dpois(exp(X[h,1]))
    N5[h,1] ~ dpois(exp(X[h,1]))    
    N6[h,1] ~ dpois(exp(X[h,1]))    
  }
  for( i in 2:qp1){ 
    X[1:p,i]  ~ dmnorm(A + B%*%X[1:p,(i-1)], Sigma); ## Process error
    for(j in 1:p){
      N4[j,i]  ~ dpois(exp(X[j,i])) ## Observation error
      N5[j,i]  ~ dpois(exp(X[j,i])) ## Observation error      
      N6[j,i]  ~ dpois(exp(X[j,i])) ## Observation error
    }
  }
}


# Try using the minimum distance first:

dmat <- min.dists1 # Change the distance matrix!
dinvmat <- (dmat/dmat^2)  
diag(dinvmat) <- rep(0,p)

# Retrieving the MLES obtained by using the minimum distance first:
mles4kalman <- mlesmindist[,1]
avec4kalman <- mles4kalman[1:6]
cvec4kalman <- mles4kalman[7:12]
mig4kalman  <- mles4kalman[13]
muvec4kalman <- mles4kalman[14:19]
sigsq4kalman <- mles4kalman[20]


B <- matrix(0,p,p)

for(m in 1:p){
  for(n in 1:p){
    B[m,n] <- ifelse(m==n, cvec4kalman[m], dinvmat[m,n]*mig4kalman) #
  }
}



data4kalman <- list(N4=S4mat,
                  N5=S5mat,
                  N6=S6mat,
                  qp1=ncol(S4mat), 
                  p=p,
                  Ip =  diag(p),
                  muvec = muvec4kalman,
                  sigmasq = sigsq4kalman, 
                  B=B
)

n.iter<-100000
n.adapt<-1000
n.update<-100
thin<-10
n.chains<-1

out.parms <- c("X")
kalman.mcmcw <- jags.fit(data4kalman, params=out.parms, model=MGSS_KALMAN2.0, n.chains = n.chains, 
                     n.adapt=n.adapt, n.update=n.update, n.iter = n.iter, thin=thin) 

all.chains <- mgssjags[[1]]
for(i in 2:n.chains){all.chains <- rbind(all.chains,mgssjags[[i]])}
dim(all.chains)
posterior.quantiles <- quantile.mcmc.list(mgssjags)
