# OK, tryingto see if I can code an assymetric effects model:

library(dclone); library(MASS); library(mcmcplots)

source("guess.calc.R")
source("guess.calc2.0.R")
source("randmvn.R")
source("Pois.MVN.gss.R")

load("SNKIdata.RData") # Load raw data

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

centr.dists <- read.csv("Demes_centroid distance in meters.csv")[,-1]
centr.dists <- as.matrix(centr.dists)/1000
rownames(centr.dists) <- rownames(S4mat)

min.condists <- read.csv("Demes_minimum distance in meters.csv")[,-1]
min.condists <- as.matrix(min.condists)/1000
rownames(min.condists) <- rownames(S4mat)


p <- nrow(S4mat)
demes <- rownames(S4mat) 
myinits <- list(lmuvec = log(rep(naive.guess[1],p)), lthetavec=log(rep(naive.guess[2],p)), 
                lsigmasq=log(naive.guess[3])) # remove 'mig' from inits


# Defining the "dinvmat" choosing one of the distance matrices

dmat <- centr.dists 
dinvmat <- (dmat/dmat^2)  
diag(dinvmat) <- rep(0,p)


# Doing Data Cloning

################ Next step: Data Cloning #####################

MGSS_DC3.0 <- function(){
  ##### Priors
  lmuvec    ~ dmnorm(mupmean, muprec)
  lthetavec ~ dmnorm(lthetamean, lthetaprec)
  lsigmasq ~ dnorm(0, 1)
  #mig ~ dbeta(2,2)
  mig.near ~ dnorm(0,1)
  mig.far ~ dnorm(0,1)
  
  sigmasq <- exp(lsigmasq)  
  for(pp in 1:p){
    muvec[pp] <- exp(lmuvec[pp])
    thetavec[pp] <- exp(lthetavec[pp])
    cvec[pp] <- exp(-thetavec[pp])
  }
  
  Sigma <- (1/sigmasq)*Ip
  
  for(m in 1:p){
    for(n in 1:p){
      B[m,n] <- ifelse((m==n),cvec[m],ifelse(
        (m==step.stone[m,1])&&(n==step.stone[m,2]), 
        dinvmat[m,n]*mig.near,dinvmat[m,n]*mig.far)
      ) #
      
    }

  }
  
  
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
                    dinvmat = dinvmat,
                    step.stone = matrix(c(1,4,2,4,3,5,4,1,5,3,6,3),
                                        nrow=6,ncol=2,byrow=TRUE)
)

cl.seq <- c(1,8,16);
n.iter<-100000;n.adapt<-50000;n.update<-100;thin<-10;n.chains<-3;

out.parms <- c("muvec", "cvec", "mig.near", "mig.far", "sigmasq", "A")
mgssdclone <- dc.fit(data4dclone, params=out.parms, model=MGSS_DC3.0, n.clones=cl.seq,
                      multiply="K", unchanged = c("qp1","p","mupmean","muprec", 
                      "step.stone", "lthetamean","lthetaprec", "Ip","dinvmat"),
                      n.chains = n.chains, n.adapt=n.adapt, n.update=n.update, 
                      n.iter = n.iter, thin=thin, inits=myinits) 
dcdiag(mgssdclone)

save.image("AssymetricB-trial.RData")
 
summary(mgssdclone)

confint(mgssdclone)

mlescentrdist <- summary(mgssdclone)[[1]]
confintcentrdist <- confint(mgssdclone)


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



# Retrieving the MLES obtained by using the minimum distance first:
mles4kalman <- mlescentrdist[,1]
avec4kalman <- mles4kalman[1:6]
cvec4kalman <- mles4kalman[7:12]
mig.far4kalman  <- mles4kalman[13]
mig.near4kalman <- mles4kalman[14]
muvec4kalman <- mles4kalman[15:20]
sigsq4kalman <- mles4kalman[21]

step.stone <- matrix(c(1,4,2,4,3,5,4,1,5,3,6,3),nrow=6,ncol=2,byrow=TRUE)

B <- matrix(0,p,p)

for(m in 1:p){
  for(n in 1:p){
    B[m,n] <- ifelse((m==n),cvec4kalman[m],ifelse(
      (m==step.stone[m,1])&&(n==step.stone[m,2]), 
      dinvmat[m,n]*mig.near4kalman,dinvmat[m,n]*mig.far4kalman)
    ) #
    
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
n.chains<-3

out.parms <- c("X")
kalman.mcmcw <- jags.fit(data4kalman, params=out.parms, model=MGSS_KALMAN2.0, n.chains = n.chains, 
                         n.adapt=n.adapt, n.update=n.update, n.iter = n.iter, thin=thin) 


Xkalman <-  do.call(rbind, kalman.mcmcw)

Kalman.quants <- apply(Xkalman, 2, FUN=function(x){quantile(x, probs=c(0.025,0.5,0.975))}) 

len <- ncol(S4mat)

names.quants <- substr(colnames(Kalman.quants), start=1, stop=3)
X1s <- Kalman.quants[,which(names.quants=="X[1", arr.ind=TRUE)]
X2s <- Kalman.quants[,which(names.quants=="X[2", arr.ind=TRUE)]
X3s <- Kalman.quants[,which(names.quants=="X[3", arr.ind=TRUE)]
X4s <- Kalman.quants[,which(names.quants=="X[4", arr.ind=TRUE)]
X5s <- Kalman.quants[,which(names.quants=="X[5", arr.ind=TRUE)]
X6s <- Kalman.quants[,which(names.quants=="X[6", arr.ind=TRUE)]



par(oma=c(1,1,1,1), mar=c(4,5,2,1))
plot(1:len, X1s[2,], pch=16, type="b", col="red", ylim=c(min(Kalman.quants),max(Kalman.quants)), main ="6-Demes log-abundances", 
     xlab="Time", ylab="6-Demes log-population sizes", bty="l")
polygon(c(1:len,rev(1:len)), c(X1s[1,],rev(X1s[3,])), col=scales::alpha("red",0.1), border=NA)
points(1:len, X1s[2,],pch=16, type="b", col="red")

polygon(c(1:len,rev(1:len)), c(X2s[1,],rev(X2s[3,])), col=scales::alpha("blue",0.1), border=NA)
points(1:len, X2s[2,],pch=16, type="b", col="blue")

polygon(c(1:len,rev(1:len)), c(X3s[1,],rev(X3s[3,])), col=scales::alpha("darkgreen",0.1), border=NA)
points(1:len, X3s[2,],pch=16, type="b", col="darkgreen")

polygon(c(1:len,rev(1:len)), c(X4s[1,],rev(X4s[3,])), col=scales::alpha("yellow",0.1), border=NA)
points(1:len, X4s[2,],pch=16, type="b", col="yellow")

polygon(c(1:len,rev(1:len)), c(X5s[1,],rev(X5s[3,])), col=scales::alpha("purple",0.1), border=NA)
points(1:len, X5s[2,],pch=16, type="b", col="purple")

polygon(c(1:len,rev(1:len)), c(X6s[1,],rev(X6s[3,])), col=scales::alpha("orange",0.1), border=NA)
points(1:len, X6s[2,],pch=16, type="b", col="orange")


