# R functions and datasets to support "Modern Applied Statistics with S", 
#a book from W.N. Venables and B.D. Ripley
library(MASS); 
# R functions for Data Cloning (maximum likelihood estimation using Bayesian MCMC)
library(dclone); 
# Create plots for MCMC output
library(mcmcplots)
# Data manipulation and elegant figures
library(tidyverse)
# Palettes inspired by works at the Metropolitan Museum of Art in New York
library(MetBrewer)

rawdat <- read.csv("Code/snailkite counts 1997_2025.csv")
#head(rawdat, 12) # first 12 rows
#tail(rawdat, 12) # last 12 rows

names.rawdat <- names(rawdat)
deme.names <- c("EAST","EVER","KRV","OKEE","PP","SJM","OTHER" )
cols.with.counts <- match(deme.names, names.rawdat)

years <- unique(rawdat$year)
nyears <- length(years)
year.col <- which(names.rawdat=="year",arr.in=TRUE)


tt <- years - years[1]

# List the survey numbers:
surv.num <- unique(rawdat$survey_num)

S4mat <- S5mat <- S6mat <- matrix(0, nrow=7, ncol=nyears)

Counts.list <- list(S4mat = S4mat, S5mat=S5mat, S6mat=S6mat)
NAS.list  <- list(NA4mat = S4mat, NA5mat=S5mat, NA6mat=S6mat)
short.inds.list <- list()

for(i in 1:3){
  
  SN <- surv.num[(i+3)]
  
  ith.rows   <- which(rawdat$survey_num==SN, arr.ind=TRUE)
  
  #ith.datmat <- rawdat[ith.rows, c(year.col,cols.with.counts)]
  ith.datmat <- t(rawdat[ith.rows, cols.with.counts])
  colnames(ith.datmat) <- tt
  
  Counts.list[[i]] <- ith.datmat
  
  ith.isnamat <- 1-is.na(ith.datmat)
  NAS.list[[i]] <- ith.isnamat
  
  ith.shortind.list <- list()
  for(j in 1:7){
    
    ith.shortind.list[[j]] <- which(ith.isnamat[j,]==1,arr.ind=TRUE)
    
  }
  names(ith.shortind.list) <- deme.names
  short.inds.list[[i]] <- ith.shortind.list 
  
}
names(short.inds.list) <- c("short.inds4","short.inds5", "short.inds6")


guess.calc <- function(Yobs,Tvec){
  
  #  For calculations, time starts at zero.
  T.t <-Tvec-Tvec[1];
  #  Number of time series transitions, q.
  q <- length(Yobs)-1;
  #  q+1 gives the length of the time series.
  qp1 <- q+1;
  #  Time intervals
  S.t <- T.t[2:qp1]-T.t[1:q]; 
  #  Mean of the observations
  Ybar <- mean(Yobs);        
  # Variance of the observations
  Yvar <- sum((Yobs-Ybar)*(Yobs-Ybar))/q;
  # Initial mu estimate (at stationary distribution)
  mu1 <- Ybar;
  
  # Kludge an initial value for theta based on mean of Y(t+s) given Y(t).
  th1 <- -mean(log(abs((Yobs[2:qp1]-mu1)/(Yobs[1:q]-mu1)))/S.t);            
  # Moment estimate using stationary
  bsq1 <- 2*th1*Yvar/(1+2*th1);         
  # Observation error variance, assumed as first guess as betasq=tausq.
  tsq1 <- bsq1;                         
  
  # What to do if initial guesses is three 0's (or NAs)? Assume arbitrary values
  three0s <- sum(c(th1,bsq1,tsq1))
  if(three0s==0|is.na(three0s)){
    th1 <- 0.5;
    bsq1 <- 0.09; 
    tsq1 <- 0.23;
  }
  
  out1 <- c(th1,bsq1,tsq1);
  
  # What to do if initial guesses are too little? Assume arbitrary values
  if(sum(out1<1e-7)>=1){
    out1 <- c(0.5,0.09,0.23)
  }
  
  out <- c(mu1,out1);
  return(abs(out))
}


guess.calc2.0<- function(TimeAndNs){
  
  newmat <- TimeAndNs 
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

randmvn <- function(n, mu.vec, cov.mat){
  
  # Save the length of the mean vector of the multivariate normal distribution to sample
  p         <- length(mu.vec);
  # The Cholesky decomposition
  #(factorization of a real symmetric positive-definite sqr matrix)
  Tau       <- chol(cov.mat, pivot=TRUE);
  # generate normal deviates outside loop
  Zmat      <- matrix(rnorm(n=p*n,mean=0,sd=1),nrow=p,ncol=n);
  
  # empty matrix
  out       <- matrix(0,nrow=p,ncol=n);
  # iterate
  for(i in 1:n){
    Z       <- Zmat[,i];
    out[,i] <- t(Tau)%*%Z + mu.vec
  }
  
  return(out)
}

S4mat <- Counts.list$S4mat[1:6,1:28] # Cut the last time step: no data anywhere
S5mat <- Counts.list$S5mat[1:6,1:28]
S6mat <- Counts.list$S6mat[1:6,1:28]

NA4mat <- NAS.list$NA4mat[1:6,1:28]
NA5mat <- NAS.list$NA5mat[1:6,1:28]
NA6mat <- NAS.list$NA6mat[1:6,1:28]

ts.4guess  <- S4mat[3,]
tvec4guess  <- as.numeric(colnames(S4mat))
onets4guess <- cbind(tvec4guess, ts.4guess) # No NAS for KRV # log(ts.4guess) changed to ts.4guess
naive.guess <- guess.calc2.0(TimeAndNs = onets4guess)
naive.guess

# Doing the Univariate model fit
p <- nrow(S4mat)
myinits <- list(lmu = log(naive.guess[1]), # mu ≈ mean(log(ts.4guess)),
                ltheta = log(naive.guess[2]), # why theta in log?
                lsigma = log(sqrt(naive.guess[3])))


MGSS_DC6.0 <- function(){
  ##### ------------------------------------------------------------------------
  ##### Priors
  lmu    ~ dnorm(6, 1/10);
  ltheta ~ dnorm(0.5, 1/10);
  lsigma ~ dnorm(0, 1);
  
  mu <- exp(lmu);
  theta <- exp(ltheta);
  sigma <- exp(lsigma);
  sigmasq <- pow(sigma,2);
  
  
  #tausq ~ dchisq(df = 1);
  car.cap <- exp(mu);
  cc <- exp(-theta);
  a  <- mu*(1-cc);
  Vinf <- sigmasq/(1-exp(-2*theta));
  InfPrec <- 1/Vinf
  
  
  ##### ----------------------------------------- -------------------------------
  ##### Likelihood
  for( h in 1:H){
    #Process No
    X[1,h] ~ dnorm(mu, InfPrec)
    
    # First observation
    N4[1,h] ~ dpois(pdet*exp(X[1,h]))
    N5[1,h] ~ dpois(pdet*exp(X[1,h]))      
    N6[1,h] ~ dpois(pdet*exp(X[1,h]))
    
    for(j in 2:len){
      X[j,h] ~ dnorm(a + cc*X[(j-1),h], 1/sigmasq)

      # jth observations
      N4[j,h] ~ dpois(pdet*exp(X[j,h]))
      N5[j,h] ~ dpois(pdet*exp(X[j,h]))      
      N6[j,h] ~ dpois(pdet*exp(X[j,h]))
    }
    
  }
  
}

len <- ncol(S4mat)
deme.names <- row.names(S4mat)

cool.palette <- as.vector(met.brewer("Homer2",p))


all.ps <- read.csv(file="Code/p.estimates.wt.7.24.25.csv")
names(all.ps)
#View(all.ps)

# I am going to leave it as a matrix of ps in case we want to make them 
# time dependent later
phats <- matrix(nrow=6,ncol=28)
row.names(phats) <- deme.names
for(i in 1:6){
  subdat <- all.ps[all.ps$stratum==deme.names[i],]
  phats[i,] <- subdat$estimate
}


# KRV
pop.num <- 3

ith.pdet <- phats[pop.num,1] # change to [1:6,1] one at a time for each population


data4dclone <- list(H=1, 
                    N4=dcdim(data.matrix(S4mat[pop.num,])), 
                    N5=dcdim(data.matrix(S5mat[pop.num,])), 
                    N6=dcdim(data.matrix(S6mat[pop.num,])), 
                    len=len,
                    pdet = ith.pdet
) 


cl.seq <- c(1,4,16)
n.iter<-100000;n.adapt<-50000;n.update<-100;thin<-10;n.chains<-3;

out.parms <- c("a", "cc", "sigmasq", "mu")
mgssdclone6<- dc.fit(data4dclone, 
                     params=out.parms, 
                     model=MGSS_DC6.0, 
                     n.clones=cl.seq,
                     multiply="H", 
                     unchanged = c("len","pdet"),
                     n.chains = n.chains, 
                     n.adapt=n.adapt, 
                     n.update=n.update,
                     n.iter = n.iter, 
                     thin=thin,
                     inits=myinits
                     ) 

# First pdetect version is saved here:
#saveRDS(mgssdclone2, "Code/mgssdclone2.RDS")

# Second pdetect version with two constant ps over time:
#saveRDS(mgssdclone6, "Code/mgssdclone6-twopdetects.RDS")

krv.mcmc <- mgssdclone6
dcdiag(krv.mcmc)
summary(krv.mcmc)
confint(krv.mcmc)

# > dcdiag(krv.mcmc)
# n.clones lambda.max   ms.error r.squared  r.hat
# 1        1 4.83318884 483.573963 0.9077828 1.0002
# 2        4 0.15768869  43.319421 0.5788588     NA
# 3       16 0.03681289   7.900878 0.1983969     NA
# > summary(krv.mcmc)
# 
# Iterations = 50110:150100
# Thinning interval = 10 
# Number of chains = 3 
# Sample size per chain = 10000 
# Number of clones = 16
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean      SD   DC SD  Naive SE Time-series SE R hat
# a       1.4726 0.18661 0.74644 0.0010774      0.0010680     1
# cc      0.7023 0.03659 0.14638 0.0002113      0.0002092     1
# mu      4.9453 0.08082 0.32326 0.0004666      0.0004666     1
# sigmasq 0.2839 0.01956 0.07824 0.0001129      0.0001183     1
# 
# 2. Quantiles for each variable:
#   
#   2.5%    25%    50%    75%  97.5%
# a       1.1124 1.3467 1.4721 1.5973 1.8416
# cc      0.6301 0.6778 0.7025 0.7271 0.7733
# mu      4.7825 4.8930 4.9464 5.0001 5.1009
# sigmasq 0.2483 0.2704 0.2831 0.2964 0.3249
# 
# > confint(krv.mcmc)
# 2.5 %    97.5 %
#   a       0.009641183 2.9356365
# cc      0.415422246 0.9892179
# mu      4.311711694 5.5788855
# sigmasq 0.130584447 0.4372972

### OKEE

pop.num <- 4

ith.pdet <- phats[pop.num,1] # change to [1:6,1] one at a time for each population


data4dclone <- list(H=1, 
                    N4=dcdim(data.matrix(S4mat[pop.num,])), 
                    N5=dcdim(data.matrix(S5mat[pop.num,])), 
                    N6=dcdim(data.matrix(S6mat[pop.num,])), 
                    len=len,
                    pdet = ith.pdet
) 


cl.seq <- c(1,4,16)
n.iter<-100000;n.adapt<-50000;n.update<-100;thin<-10;n.chains<-3;

out.parms <- c("a", "cc", "sigmasq", "mu")
mgssdclone6<- dc.fit(data4dclone, 
                     params=out.parms, 
                     model=MGSS_DC6.0, 
                     n.clones=cl.seq,
                     multiply="H", 
                     unchanged = c("len","pdet"),
                     n.chains = n.chains, 
                     n.adapt=n.adapt, 
                     n.update=n.update,
                     n.iter = n.iter, 
                     thin=thin,
                     inits=myinits
) 


okee.mcmc <- mgssdclone6
dcdiag(okee.mcmc)
summary(okee.mcmc)
confint(okee.mcmc)

# > dcdiag(okee.mcmc)
# n.clones lambda.max   ms.error r.squared r.hat
# 1        1 6.78069993 453.798229 0.8893946    NA
# 2        4 0.31463207  29.459569 0.4995579    NA
# 3       16 0.07267919   7.372268 0.1928953    NA
# > summary(okee.mcmc)
# 
# Iterations = 50110:150100
# Thinning interval = 10 
# Number of chains = 3 
# Sample size per chain = 10000 
# Number of clones = 16
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean      SD  DC SD  Naive SE Time-series SE R hat
# a       1.4442 0.15546 0.6218 0.0008975      0.0008927     1
# cc      0.6057 0.03969 0.1588 0.0002292      0.0002367     1
# mu      3.6640 0.17363 0.6945 0.0010024      0.0009886     1
# sigmasq 2.2004 0.19028 0.7611 0.0010986      0.0011118     1
# 
# 2. Quantiles for each variable:
#   
#   2.5%    25%    50%    75%  97.5%
# a       1.1430 1.3390 1.4422 1.5474 1.7543
# cc      0.5271 0.5792 0.6061 0.6328 0.6825
# mu      3.3234 3.5477 3.6639 3.7797 4.0062
# sigmasq 1.8571 2.0681 2.1895 2.3226 2.5983
# 
# > confint(okee.mcmc)
# 2.5 %    97.5 %
#   a       0.2254601 2.6629653
# cc      0.2945165 0.9168684
# mu      2.3028251 5.0252659
# sigmasq 0.7086449 3.6921353


#### EVER

pop.num <- 2

ith.pdet <- phats[pop.num,1] # change to [1:6,1] one at a time for each population


data4dclone <- list(H=1, 
                    N4=dcdim(data.matrix(S4mat[pop.num,])), 
                    N5=dcdim(data.matrix(S5mat[pop.num,])), 
                    N6=dcdim(data.matrix(S6mat[pop.num,])), 
                    len=len,
                    pdet = ith.pdet
) 


cl.seq <- c(1,4,16)
n.iter<-100000;n.adapt<-50000;n.update<-100;thin<-10;n.chains<-3;

out.parms <- c("a", "cc", "sigmasq", "mu")
mgssdclone6<- dc.fit(data4dclone, 
                     params=out.parms, 
                     model=MGSS_DC6.0, 
                     n.clones=cl.seq,
                     multiply="H", 
                     unchanged = c("len","pdet"),
                     n.chains = n.chains, 
                     n.adapt=n.adapt, 
                     n.update=n.update,
                     n.iter = n.iter, 
                     thin=thin,
                     inits=myinits
) 

# First pdetect version is saved here:
#saveRDS(mgssdclone2, "Code/mgssdclone2.RDS")

# Second pdetect version with two constant ps over time:
#saveRDS(mgssdclone6, "Code/mgssdclone6-twopdetects.RDS")

ever.mcmc <- mgssdclone6
dcdiag(ever.mcmc)
summary(ever.mcmc)
confint(ever.mcmc)

# > dcdiag(ever.mcmc)
# n.clones lambda.max   ms.error r.squared    r.hat
# 1        1 16.0784394 201.821616 0.8061898       NA
# 2        4  1.1022877   9.714287 0.2404850 1.000238
# 3       16  0.1206531  16.635965 0.3689867       NA
# > summary(ever.mcmc)
# 
# Iterations = 50110:150100
# Thinning interval = 10 
# Number of chains = 3 
# Sample size per chain = 10000 
# Number of clones = 16
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean      SD  DC SD  Naive SE Time-series SE  R hat
# a       0.6697 0.18290 0.7316 0.0010560      0.0010580 1.0000
# cc      0.8660 0.03255 0.1302 0.0001879      0.0001888 0.9999
# mu      4.9669 0.32702 1.3081 0.0018881      0.0018722 1.0000
# sigmasq 0.9104 0.06681 0.2672 0.0003857      0.0003863 1.0001
# 
# 2. Quantiles for each variable:
#   
#   2.5%    25%    50%    75%  97.5%
# a       0.3399 0.5419 0.6597 0.7877 1.0533
# cc      0.7981 0.8449 0.8676 0.8888 0.9253
# mu      4.2541 4.7717 4.9917 5.1890 5.5388
# sigmasq 0.7885 0.8639 0.9071 0.9530 1.0511
# 
# > confint(ever.mcmc)
# 2.5 %   97.5 %
#   a       -0.7642529 2.103617
# cc       0.6108245 1.121237
# mu       2.4030725 7.530665
# sigmasq  0.3866487 1.434136


##### EAST

pop.num <- 1

ith.pdet <- phats[pop.num,1] # change to [1:6,1] one at a time for each population


data4dclone <- list(H=1, 
                    N4=dcdim(data.matrix(S4mat[pop.num,])), 
                    N5=dcdim(data.matrix(S5mat[pop.num,])), 
                    N6=dcdim(data.matrix(S6mat[pop.num,])), 
                    len=len,
                    pdet = ith.pdet
) 


cl.seq <- c(1,4,16)
n.iter<-100000;n.adapt<-50000;n.update<-100;thin<-10;n.chains<-3;

out.parms <- c("a", "cc", "sigmasq", "mu")
mgssdclone6<- dc.fit(data4dclone, 
                     params=out.parms, 
                     model=MGSS_DC6.0, 
                     n.clones=cl.seq,
                     multiply="H", 
                     unchanged = c("len","pdet"),
                     n.chains = n.chains, 
                     n.adapt=n.adapt, 
                     n.update=n.update,
                     n.iter = n.iter, 
                     thin=thin,
                     inits=myinits
) 


east.mcmc <- mgssdclone6
dcdiag(east.mcmc)
summary(east.mcmc)
confint(east.mcmc)


save.image("RegionByRegionFits.RData")


# > dcdiag(east.mcmc)
# n.clones lambda.max   ms.error r.squared    r.hat
# 1        1 1.65406848 365.831196 0.8761601 1.000406
# 2        4 0.16439338  16.116497 0.3512584       NA
# 3       16 0.03789413   6.184788 0.1586350       NA
# > summary(east.mcmc)
# 
# Iterations = 50110:150100
# Thinning interval = 10 
# Number of chains = 3 
# Sample size per chain = 10000 
# Number of clones = 16
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean      SD  DC SD  Naive SE Time-series SE R hat
# a       2.1670 0.18876 0.7550 0.0010898      0.0010843     1
# cc      0.4866 0.04367 0.1747 0.0002521      0.0002508     1
# mu      4.2209 0.08193 0.3277 0.0004730      0.0004750     1
# sigmasq 0.7973 0.05974 0.2390 0.0003449      0.0003449     1
# 
# 2. Quantiles for each variable:
#   
#   2.5%    25%    50%    75%  97.5%
# a       1.7987 2.0399 2.1654 2.2936 2.5381
# cc      0.4001 0.4573 0.4867 0.5162 0.5718
# mu      4.0605 4.1662 4.2209 4.2753 4.3827
# sigmasq 0.6890 0.7557 0.7946 0.8357 0.9229
# 
# > confint(east.mcmc)
# 2.5 %    97.5 %
#   a       0.6871306 3.6468591
# cc      0.1442512 0.8289559
# mu      3.5786307 4.8632497
# sigmasq 0.3288946 1.2656555
# > 