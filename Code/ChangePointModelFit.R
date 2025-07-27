# libraries:
library(MASS); 
# R functions for Data Cloning (maximum likelihood estimation using Bayesian MCMC)
library(dclone); 
# Create plots for MCMC output
library(mcmcplots)
# Data manipulation and elegant figures
library(tidyverse)
# Palettes inspired by works at the Metropolitan Museum of Art in New York
library(MetBrewer)


# Is there a detectable change in the population dynamics and if so, what changed?

# m1 => no change.  a1=a2, c1=c2, sigsq1=sigsq2 (3 parms + 1 obs error param = 4 params) NULL MODEL
# m2 => a1=a2, c1=c2,sigsq1, sigs2. Only env. noise changes. (4 parms + 1 =5 parms)
# m3 => a1=a2, c1, c2, sigsq. 4+1 params =5 params
# m4 => a1, a2, c1=c2, sigsq. 4+1 params = 5 params
# m5 => a1, a2, c1, c2, sigsq. 5+1 parms = 6 params
# m6 => a1=a2,c1,c2,sigsq1,sigsq2 = 5 + 1 params = 6 params
# m7 => a1,a2,c1=c2, sigsq1, sigsq2, = 5 + 1 params = 6 params
# m8 => a1,a2, c1,c2, sigsq1,sigsq2. 6 params + 1 obs error=7 params FULL MODEL

# aesthetic stuff:

cool.palette <- as.vector(met.brewer("Homer2",6))

source("BPGSSToolkit.R")

kalman.ns <- read.csv(file="Code/TimeSeriesMLEs-pdetect.csv")
kalman.ns$years <- 1997:2024
NsLCL <- kalman.ns$NsLCL.2
NsUCL <- kalman.ns$NsUCL.2
years <- kalman.ns$years
Nhats <- kalman.ns$NsMLE.2
lNhats <- log(Nhats)

# First we check where the most likely breakpoint is

guess8 <- c(log(1.2), log(1.2), atanh(0.7934),atanh(0.7934), log(0.09725), log(0.09725), log(0.2315))
tau1s <- 7:13
neglls <- rep(0, length(tau1s))
for(i in 1:length(tau1s)){
  
  ith.opt <- optim(par=guess8, fn=model.profmle, method="Nelder-Mead", 
                   data.vec=lNhats, M=tau1s[i], model.flag="m8")
  neglls[i] <- ith.opt$value
  
}

print(cbind(1:length(tau1s), tau1s,years[tau1s],neglls))


par(mfrow=c(1,1))
plot(years, Nhats, pch=16, type="b", col=cool.palette[1],
     main ="Estimated total local population sizes", 
     ylim= c(300,2200),
     xlab="Years", 
     ylab="Estimated pooled local population sizes", bty="l")
polygon(c(years,rev(years)), c(NsLCL,rev(NsUCL)), 
        col=scales::alpha(cool.palette[1],0.1), border=NA)
points(years, Nhats,pch=16, type="b", col=cool.palette[1])
abline(v=2007, lwd=2, col="blue")



year.break <- years[tau1s[which.min(neglls)]] # 2006
position.break <-tau1s[which.min(neglls)]



# Fitting the 8 models using the MLE for the breakpoint = 2007 (11th position in the vector:

tau1 <- position.break
guess8 <- c(log(1.2), log(1.2), atanh(0.7934),atanh(0.7934), log(0.09725), log(0.09725), log(0.2315))
trial.optim <- optim(par=guess8, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m8")
mles <- c(exp(trial.optim$par[1:2]), tanh(trial.optim$par[3:4]), exp(trial.optim$par[5:7]))
mles
#  4.558943432 1.190043826 0.336167225 0.827966890 0.051355877 0.044231718 0.000681639
mu1.hat <- mles[1]/(1-mles[3])
mu2.hat <- mles[2]/(1-mles[4])
t.half.hat <- log(0.5)/log(abs(mles[4]))
m8.bic <- 2*trial.optim$value + 7*log(length(kalman.ns$years)-1)
m8.bic # 18.07312

guess7 <- c(log(0.3929), log(0.3929), atanh(0.7934), log(0.09725), log(0.09725), log(0.2315))
trial7.optim <- optim(par=guess7, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m7")
mles7 <- c(exp(trial7.optim$par[1:2]), tanh(trial7.optim$par[3]), exp(trial7.optim$par[4:6]))
mles7
#  21.64012805 1.64529591 0.76064004 0.02867769 0.03000849 0.01296997
m7.bic <- 2*trial7.optim$value + 6*log(length(kalman.ns$years)-1)
m7.bic # 15.65349

guess6 <- c(log(0.3929), atanh(0.7934),atanh(0.7934), log(0.09725), log(0.09725), log(0.2315))
trial6.optim <- optim(par=guess6, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m6")
mles6 <- c(exp(trial6.optim$par[1]), tanh(trial6.optim$par[2:3]), exp(trial6.optim$par[4:6]))
mles6
# 0.62668298 0.90762240 0.91121956 0.02709834 0.02474586 0.01514553
m6.bic <- 2*trial6.optim$value + 6*log(length(kalman.ns$years)-1)
m6.bic # 16.45684

guess5 <- c(log(0.3929),log(0.3929), atanh(0.7934),atanh(0.7934),  log(0.09725), log(0.2315))
trial5.optim <- optim(par=guess5, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m5")
mles5 <- c(exp(trial5.optim$par[1:2]), tanh(trial5.optim$par[3:4]), exp(trial5.optim$par[5:6]))
mles5
#  1.98182958 0.85648546 0.70969128 0.87616071 0.02510353 0.01544866
m5.bic <- 2*trial5.optim$value + 6*log(length(kalman.ns$years)-1)
m5.bic # 15.2313

guess4 <- c(log(0.3929),log(0.3929), atanh(0.7934), log(0.09725), log(0.2315))
trial4.optim <- optim(par=guess4, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m4")
mles4 <- c(exp(trial4.optim$par[1:2]), tanh(trial4.optim$par[3]), exp(trial4.optim$par[4:5]))
mles4
#  0.66134211 0.68451609 0.90277469 0.01959724 0.01989169
m4.bic <- 2*trial4.optim$value + 5*log(length(kalman.ns$years)-1)
m4.bic # 12.86567

guess3 <- c(log(0.3929),atanh(0.7934), atanh(0.7934), log(0.09725), log(0.2315))
trial3.optim <- optim(par=guess3, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m3")
mles3 <- c(exp(trial3.optim$par[1]), tanh(trial4.optim$par[2:3]), exp(trial4.optim$par[4:5]))
mles3
#  1.03447555 -0.36187619  0.90277469  0.01959724  0.01989169
m3.bic <- 2*trial3.optim$value + 5*log(length(kalman.ns$years)-1)
m3.bic # 12.29269

guess2 <- c(log(0.3929),atanh(0.7934), log(0.09725),log(0.09725), log(0.2315))
trial2.optim <- optim(par=guess2, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m2")
mles2 <- c(exp(trial2.optim$par[1]), tanh(trial2.optim$par[2]), exp(trial2.optim$par[3:5]))
mles2
#  1.31784271 0.80795941 0.02293212 0.02815556 0.01628879
m2.bic <- 2*trial2.optim$value + 5*log(length(kalman.ns$years)-1)
m2.bic # 12.23604

guess1 <- c(log(0.3929),atanh(0.7934), log(0.09725), log(0.2315))
trial1.optim <- optim(par=guess1, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m1")
mles1 <- c(exp(trial1.optim$par[1]), tanh(trial1.optim$par[2]), exp(trial1.optim$par[3:4]))
mles1
#  1.34823483 0.80362599 0.02695655 0.01543765
m1.bic <- 2*trial1.optim$value + 4*log(length(kalman.ns$years)-1)
m1.bic # 8.978214

all.optims <- list(opt8=trial.optim, opt7=trial7.optim, opt6=trial6.optim, 
                   opt5=trial5.optim, opt4=trial4.optim, opt3=trial3.optim,
                   opt2=trial2.optim, opt1=trial1.optim)
all.bics <- c(m8.bic,m7.bic,m6.bic,m5.bic,m4.bic,m3.bic,m2.bic,m1.bic)

# Summary:  no change is best followed by change in only sigsq and change in only
# the maximum growth rate and not the density dependence nor the environmental noise variance/

mod.flags <- paste0("m",8:1)
summary.table <- matrix(0,nrow=8, ncol=7+2)
colnames(summary.table) <- c("a_1", "a_2", "c_1", "c_2", "sigma^{2}_1", "sigma^{2}_2","tau^{2}", "-logL","BIC")
for(i in 1:8){
  
  summary.table[i,1:7] <- optim.out.organizer(opt.object=all.optims[[i]], model.flag=mod.flags[i])
  summary.table[i,8] <- all.optims[[i]]$value
  summary.table[i,9] <- all.bics[i]
  
}

print(summary.table)

# > print(summary.table)
# a_1    a_2    c_1    c_2 sigma^{2}_1 sigma^{2}_2 tau^{2}     -logL       BIC
# [1,] 4.5589 1.1900 0.3362 0.8280      0.0514      0.0442  0.0007 -2.498869 18.073121
# [2,] 1.6401 1.6453 0.7606 0.7606      0.0287      0.0300  0.0130 -2.060763 15.653495
# [3,] 0.6267 0.6267 0.9076 0.9112      0.0271      0.0247  0.0151 -1.659090 16.456841
# [4,] 1.9818 0.8565 0.7097 0.8762      0.0251      0.0251  0.0154 -2.271859 15.231302
# [5,] 0.6613 0.6845 0.9028 0.9028      0.0196      0.0196  0.0199 -1.806759 12.865667
# [6,] 1.0345 1.0345 0.8482 0.8511      0.0256      0.0256  0.0162 -2.093249 12.292687
# [7,] 1.3178 1.3178 0.8080 0.8080      0.0229      0.0282  0.0163 -2.121574 12.236037
# [8,] 1.3482 1.3482 0.8036 0.8036      0.0270      0.0270  0.0154 -2.102567  8.978214



mod.descriptors <- c("m8= all different", "m7 = two a's, one c, two sigsq's", "m6 = one a,two cs, two sigsqs",
                     "m5 = two as, two cs, onsigsq", "m4 = two as,one c, one sigsq", "m3 = one a, two cs, one sigsq",
                     "m2 = one a, one c, two sigsqs", "m1 = one a, one c, one sigsq")
summarydf <- data.frame(model=mod.descriptors, summary.table)

save.image("Code/ChangePointOnKalmans.RData")

write.csv(summarydf, file="Code/SummmaryChangePointOnKalmans.csv")
