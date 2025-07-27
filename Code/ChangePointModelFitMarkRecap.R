cool.palette <- as.vector(met.brewer("Homer2",6))

source("BPGSSToolkit.R")


mark.dat <- read.csv("Code/2024_1997_popsize_ci.csv", header=TRUE)
years <- mark.dat$Year
Nhats <- mark.dat$Nsuper
lNhats <- log(Nhats)

NsLCL <- mark.dat$SuperLCL
NsUCL <- mark.dat$SuperUCL

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
     main ="Mark-Recapture Super Population sizes", 
     ylim= c(300,3500),
     xlab="Years", 
     ylab="Estimated Super population sizes", bty="l")
polygon(c(years,rev(years)), c(NsLCL,rev(NsUCL)), 
        col=scales::alpha(cool.palette[1],0.1), border=NA)
points(years, Nhats,pch=16, type="b", col=cool.palette[1])
abline(v=2007, lwd=2, col="blue")


year.break <- years[tau1s[which.min(neglls)]] # 2006
position.break <-tau1s[which.min(neglls)]


tau1 <- position.break
guess8 <- c(log(1.2), log(1.2), atanh(0.7934),atanh(0.7934), log(0.09725), log(0.09725), log(0.2315))
trial.optim <- optim(par=guess8, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m8")
mles <- c(exp(trial.optim$par[1:2]), tanh(trial.optim$par[3:4]), exp(trial.optim$par[5:7]))
mles
#  0.2887660878 1.7699237630 0.9577937521 0.7773560141 0.0596112130 0.0375399686 0.0007079917
mu1.hat <- mles[1]/(1-mles[3])
mu2.hat <- mles[2]/(1-mles[4])
t.half.hat <- log(0.5)/log(abs(mles[4]))
m8.bic <- 2*trial.optim$value + 7*log(length(years)-1)
m8.bic # 18.44991

guess7 <- c(log(0.3929), log(0.3929), atanh(0.7934), log(0.09725), log(0.09725), log(0.2315))
trial7.optim <- optim(par=guess7, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m7")
mles7 <- c(exp(trial7.optim$par[1:2]), tanh(trial7.optim$par[3]), exp(trial7.optim$par[4:6]))
mles7
#  1.401031e+00 1.501953e+00 8.100736e-01 7.128784e-02 3.058283e-02 2.421166e-05
m7.bic <- 2*trial7.optim$value + 6*log(length(years)-1)
m7.bic # 14.55828

guess6 <- c(log(0.3929), atanh(0.7934),atanh(0.7934), log(0.09725), log(0.09725), log(0.2315))
trial6.optim <- optim(par=guess6, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m6")
mles6 <- c(exp(trial6.optim$par[1]), tanh(trial6.optim$par[2:3]), exp(trial6.optim$par[4:6]))
mles6
# 1.587441e+00 7.748851e-01 8.024916e-01 1.039282e-01 5.125618e-11 2.864744e-02
m6.bic <- 2*trial6.optim$value + 6*log(length(years)-1)
m6.bic # 16.2696

guess5 <- c(log(0.3929),log(0.3929), atanh(0.7934),atanh(0.7934),  log(0.09725), log(0.2315))
trial5.optim <- optim(par=guess5, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m5")
mles5 <- c(exp(trial5.optim$par[1:2]), tanh(trial5.optim$par[3:4]), exp(trial5.optim$par[5:6]))
mles5
#  0.108236766 0.511129546 0.985556648 0.937204711 0.043692743 0.002763388
m5.bic <- 2*trial5.optim$value + 6*log(length(years)-1)
m5.bic # 19.63967

guess4 <- c(log(0.3929),log(0.3929), atanh(0.7934), log(0.09725), log(0.2315))
trial4.optim <- optim(par=guess4, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m4")
mles4 <- c(exp(trial4.optim$par[1:2]), tanh(trial4.optim$par[3]), exp(trial4.optim$par[4:5]))
mles4
#  1.080944e+00 1.173539e+00 8.532776e-01 5.012582e-02 3.798358e-05
m4.bic <- 2*trial4.optim$value + 5*log(length(years)-1)
m4.bic # 13.43785

guess3 <- c(log(0.3929),atanh(0.7934), atanh(0.7934), log(0.09725), log(0.2315))
trial3.optim <- optim(par=guess3, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m3")
mles3 <- c(exp(trial3.optim$par[1]), tanh(trial4.optim$par[2:3]), exp(trial4.optim$par[4:5]))
mles3
#  1.045821e-01 1.586715e-01 8.532776e-01 5.012582e-02 3.798358e-05
m3.bic <- 2*trial3.optim$value + 5*log(length(years)-1)
m3.bic # 17.39172

guess2 <- c(log(0.3929),atanh(0.7934), log(0.09725),log(0.09725), log(0.2315))
trial2.optim <- optim(par=guess2, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m2")
mles2 <- c(exp(trial2.optim$par[1]), tanh(trial2.optim$par[2]), exp(trial2.optim$par[3:5]))
mles2
#  1.314647e+00 8.292714e-01 7.532300e-02 3.278982e-02 8.152236e-20
m2.bic <- 2*trial2.optim$value + 5*log(length(years)-1)
m2.bic # 13.08028

guess1 <- c(log(0.3929),atanh(0.7934), log(0.09725), log(0.2315))
trial1.optim <- optim(par=guess1, fn=model.profmle, method="Nelder-Mead", data.vec=lNhats, M=tau1, model.flag="m1")
mles1 <- c(exp(trial1.optim$par[1]), tanh(trial1.optim$par[2]), exp(trial1.optim$par[3:4]))
mles1
#  9.814724e-01 8.709476e-01 5.258244e-02 1.284786e-10
m1.bic <- 2*trial1.optim$value + 4*log(length(years)-1)
m1.bic # 11.59899

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

# print(summary.table)
# a_1    a_2    c_1    c_2 sigma^{2}_1 sigma^{2}_2 tau^{2}       -logL      BIC
# [1,] 0.2888 1.7699 0.9578 0.7774      0.0596      0.0375  0.0007 -2.31047182 18.44991
# [2,] 1.4010 1.5020 0.8101 0.8101      0.0713      0.0306  0.0000 -2.60837106 14.55828
# [3,] 1.5874 1.5874 0.7749 0.8025      0.1039      0.0000  0.0286 -1.75271304 16.26960
# [4,] 0.1082 0.5111 0.9856 0.9372      0.0437      0.0437  0.0028 -0.06767679 19.63967
# [5,] 1.0809 1.1735 0.8533 0.8533      0.0501      0.0501  0.0000 -1.52066819 13.43785
# [6,] 0.1046 0.1046 0.9857 0.9961      0.0494      0.0494  0.0023  0.45626620 17.39172
# [7,] 1.3146 1.3146 0.8293 0.8293      0.0753      0.0328  0.0000 -1.69945110 13.08028
# [8,] 0.9815 0.9815 0.8709 0.8709      0.0526      0.0526  0.0000 -0.79217715 11.59899

mod.descriptors <- c("m8= all different", "m7 = two a's, one c, two sigsq's", "m6 = one a,two cs, two sigsqs",
                     "m5 = two as, two cs, onsigsq", "m4 = two as,one c, one sigsq", "m3 = one a, two cs, one sigsq",
                     "m2 = one a, one c, two sigsqs", "m1 = one a, one c, one sigsq")
summarydf <- data.frame(model=mod.descriptors, summary.table)

write.csv(summarydf, file="Code/SummmaryChangePointOnMark-RecapNs.csv")

save.image("Code/ChangePointOnMark-RecapNs.RData")



