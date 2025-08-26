GSM.sim <- function(a,cc,sigsq,no, len,nreps){
  
  b <- cc-1
  pop.mat <- matrix(0,nrow=nreps,ncol=len)
  pop.mat[,1] <- rep(no,nreps)
  
  for(i in 2:len){
    
    for(j in 1:nreps){
      
      if(pop.mat[j,(i-1)]<=0){pop.mat[j,(i-1)] <- 0}else{
        pop.mat[j,i] <- pop.mat[j,(i-1)]*exp(a+b*log(pop.mat[j,(i-1)]) + rnorm(n=1,mean=0,sd=sqrt(sigsq)) ) 
      }
      
    }
    
  }
  
  return(t(pop.mat))
  
}



a <- 0.75
cc <- 0.88
K <- exp(a/(1-cc))
print(K)
sigsq <- 0.06#0.09725
statio.var <- sigsq/(1-cc^2)

no <- 15
len   <- 100
nreps <- 10 

trial <- GSM.sim(a=a, cc=cc, sigsq=sigsq, no=no, len=len,nreps=nreps)
trial.det <- GSM.sim(a=a, cc=cc, sigsq=0, no=15, len=len,nreps=1)
E.trial <- GSM.sim(a=a, cc=cc, sigsq=sigsq, no=no, len=len,nreps=50000)
E.trial.mean <- apply(E.trial,1,mean)

up.sup <- 1700
stop.plot <- 75


# Version 1: including multiple realizations of the process

par(mfrow=c(2,2), oma=c(5,5,4,1), mar=c(2,2,2,1))
# Upper left
plot(1:stop.plot,trial.det[1:stop.plot], type="l", lty=1,  col="blue", lwd=2,
     xlab="", ylab="", cex.lab=0.9, main="", cex.main=0.9,
     bty="l", xlim=c(0,len), cex.lab=1.25, axes=FALSE, ylim=c(0,up.sup))
axis(side=1, at=c(0,seq(20,len,by=20)))
axis(side=2, at=c(0,seq(200,up.sup,by=200)))
text(stop.plot+3,K-8, expression(K[b]), cex=1.15)
points(stop.plot,K,pch=16, col="blue", cex=0.8)
text(5,1400,"A.",cex=1.5)


# Upper right
matplot(1:len,trial, type="l",col=rep("darkgrey",5), lty=rep(1,5), 
xlab="", ylab="", cex.lab=0.9, main="", cex.main=0.9,
bty="l", xlim=c(0,len), cex.lab=1.25, axes=FALSE, ylim=c(0,up.sup))
axis(side=1, at=c(0,seq(20,len,by=20)))
axis(side=2, at=c(0,seq(200,up.sup,by=200)))
text(5,1400,"B.",cex=1.5)


# Lower left
matplot(1:stop.plot,trial[1:stop.plot,], type="l",col=rep("darkgrey",5), lty=rep(1,5), 
        xlab="", ylab="", cex.lab=0.9, main="", cex.main=0.9,
        bty="l", xlim=c(0,len), cex.lab=1.25, axes=FALSE, ylim=c(0,up.sup))
axis(side=1, at=c(0,seq(20,len,by=20)))
axis(side=2, at=c(0,seq(200,up.sup,by=200)))
# Support and density of the side-wise pdf:
supp1 <- seq(0,up.sup,by=5)
dens1 <- 15000*dlnorm(supp1,meanlog=a/(1-cc), sdlog=sqrt(statio.var))
points(c(stop.plot,stop.plot), c(0,up.sup), type="l", lty=1)
points(stop.plot+dens1, supp1, type="l", lwd=1, col="red")
text(5,1400,"C.",cex=1.5)


# Lower right
matplot(1:stop.plot,trial[1:stop.plot,], type="l",col=rep("darkgrey",5), lty=rep(1,5), 
        xlab="", ylab="", cex.lab=0.9, main="", cex.main=0.9,
        bty="l", xlim=c(0,len), cex.lab=1.25, axes=FALSE, ylim=c(0,up.sup))
axis(side=1, at=c(0,seq(20,len,by=20)))
axis(side=2, at=c(0,seq(200,up.sup,by=200)))
# Support and density of the side-wise pdf:
supp1 <- seq(0,up.sup,by=5)
dens1 <- 15000*dlnorm(supp1,meanlog=a/(1-cc), sdlog=sqrt(statio.var))
points(c(stop.plot,stop.plot), c(0,up.sup), type="l", lty=1)
points(stop.plot+dens1, supp1, type="l", lwd=1, col="red")
# Adding the deterministic trajectory
points(1:stop.plot,trial.det[1:stop.plot], type="l", lty=1,  col="blue", lwd=2)
text(stop.plot+3,K-25, expression(K[b]), cex=1.15)
points(stop.plot,K,pch=16, col="blue", cex=0.8)
# Adding the trajectory of the Expected Value
points(1:stop.plot,E.trial.mean[1:stop.plot], type="l", lty=1,  col="red", lwd=2)
points(stop.plot,exp((a/(1-cc)) + statio.var/2 ), pch=16, cex=0.8,col="red")
text(stop.plot+5,25+exp((a/(1-cc)) + statio.var/2 ),expression(paste("E[",N[t],"]", sep="")),cex=1.15 )
text(5,1400,"D.",cex=1.5)


mtext("Time", side=1, outer=TRUE, cex=1.15, line=1)
mtext("Population abundance", side=2, outer=TRUE, cex=1.15, line=1)


# Version 2: including only one realization of the process

no <- 15
len   <- 200
nreps <- 10 

trial <- GSM.sim(a=a, cc=cc, sigsq=sigsq, no=no, len=len,nreps=nreps)
trial.det <- GSM.sim(a=a, cc=cc, sigsq=0, no=no, len=len,nreps=1)
E.trial <- GSM.sim(a=a, cc=cc, sigsq=sigsq, no=no, len=len,nreps=50000)
E.trial.mean <- apply(E.trial,1,mean)

up.sup <- 1700
stop.plot <- floor((3/4)*len)



sim4plot <- 6

par(mfrow=c(2,2), oma=c(5,5,4,3), mar=c(2,2,2,2))
# Upper left
plot(1:stop.plot,trial.det[1:stop.plot], type="l", lty=1,  col="blue", lwd=2,
     xlab="", ylab="", cex.lab=0.9, main="", cex.main=0.9,
     bty="l", xlim=c(0,len), cex.lab=1.25, axes=FALSE, ylim=c(0,up.sup))
axis(side=1, at=c(0,seq(20,len,by=20)))
axis(side=2, at=c(0,seq(200,up.sup,by=200)))
text(stop.plot+3,K-8, expression(K[b]), cex=1.15)
points(stop.plot,K,pch=16, col="blue", cex=0.8)
text(5,1400,"A.",cex=1.5)


# Upper right
matplot(1:stop.plot,trial[1:stop.plot,sim4plot], type="l",col="darkgrey", lty=1, 
        xlab="", ylab="", cex.lab=0.9, main="", cex.main=0.9,
        bty="l", xlim=c(0,len), cex.lab=1.25, axes=FALSE, ylim=c(0,up.sup))
axis(side=1, at=c(0,seq(20,len,by=20)))
axis(side=2, at=c(0,seq(200,up.sup,by=200)))
# Support and density of the side-wise pdf:
supp1 <- seq(0,up.sup,by=5)
dens1 <- 25000*dlnorm(supp1,meanlog=a/(1-cc), sdlog=sqrt(statio.var))
points(c(stop.plot,stop.plot), c(0,up.sup), type="l", lty=1)
points(stop.plot+dens1, supp1, type="l", lwd=1, col="red")
text(5,1400,"B.",cex=1.5)


# Lower left
matplot(1:stop.plot,trial[1:stop.plot,sim4plot], type="l",col="darkgrey", lty=1, 
        xlab="", ylab="", cex.lab=0.9, main="", cex.main=0.9,
        bty="l", xlim=c(0,len), cex.lab=1.25, axes=FALSE, ylim=c(0,up.sup))
axis(side=1, at=c(0,seq(20,len,by=20)))
axis(side=2, at=c(0,seq(200,up.sup,by=200)))
# Support and density of the side-wise pdf:
supp1 <- seq(0,up.sup,by=5)
dens1 <- 25000*dlnorm(supp1,meanlog=a/(1-cc), sdlog=sqrt(statio.var))
points(c(stop.plot,stop.plot), c(0,up.sup), type="l", lty=1)
points(stop.plot+dens1, supp1, type="l", lwd=1, col="red")
# Adding the deterministic trajectory
points(1:stop.plot,trial.det[1:stop.plot], type="l", lty=1,  col="blue", lwd=2)
text(stop.plot+8,K-25, expression(K[b]), cex=1.15)
points(stop.plot,K,pch=16, col="blue", cex=0.8)
# Adding the trajectory of the Expected Value
points(1:stop.plot,E.trial.mean[1:stop.plot], type="l", lty=1,  col="red", lwd=2)
points(stop.plot,exp((a/(1-cc)) + statio.var/2 ), pch=16, cex=0.8,col="red")
text(stop.plot+11,25+exp((a/(1-cc)) + statio.var/2 ),expression(paste("E[",N[t],"]", sep="")),cex=1.15 )
text(5,1400,"C.",cex=1.5)



# Lower right
matplot(1:stop.plot,trial[1:stop.plot,sim4plot], type="l",col="darkgrey", lty=1, 
        xlab="", ylab="", cex.lab=0.9, main="", cex.main=0.9,
        bty="l", xlim=c(0,len), cex.lab=1.25, axes=FALSE, ylim=c(0,up.sup))
axis(side=1, at=c(0,seq(20,len,by=20)))
axis(side=2, at=c(0,seq(200,up.sup,by=200)))
# Support and density of the side-wise pdf:
supp1 <- seq(0,up.sup,by=5)
dens1 <- 25000*dlnorm(supp1,meanlog=a/(1-cc), sdlog=sqrt(statio.var))
points(c(stop.plot,stop.plot), c(0,up.sup), type="l", lty=1)
points(stop.plot+dens1, supp1, type="l", lwd=1, col="red")
text(5,1400,"D.",cex=1.5)

# Adding the third axis and polygon with p(extinction)
nc <- 300
axis(side=4, at=c(0,seq(200,up.sup,by=200)))
points(c(stop.plot,len),c(nc,nc), type="l", lty=2, lwd=1,col="red")
# Find index of intersection of pdf and cuttoff point nc:
int.ind <- which.min(((supp1)-nc)^2)
#dens1.intersect <- round((stop.plot+dens1)[int.ind],digits=2)
xs4polygon <- stop.plot+dens1[1:int.ind]   
ys.low4polygon <- supp1[1:int.ind]
ys.hi4polygon <- rep(nc,int.ind)
polygon(x=c(xs4polygon,rev(xs4polygon)), y=c(ys.low4polygon,ys.hi4polygon),
        col=scales::alpha("red",0.1), border=NA)
text(stop.plot-7,nc, expression(n[c]),cex=1.15)


mtext("Time", side=1, outer=TRUE, cex=1.15, line=1)
mtext("Population abundance", side=2, outer=TRUE, cex=1.15, line=1)



save.image("StatioPlot2.0.RData")

# Adjusting from comments in the report ####
par(mfrow=c(1,1), oma=c(4,4,3,2), mar=c(2,2,2,2))

# Left
matplot(1:stop.plot,trial[1:stop.plot, sim4plot], type="l",
        col="darkgrey", lty=1, 
        xlab="", ylab="", cex.lab=0.9, main="", cex.main=0.9,
        bty="l", xlim=c(0,len), cex.lab=1.25, axes=FALSE, ylim=c(0,up.sup-250))

# adding three axis
axis(side=1, at=c(0,seq(20,len,by=20)))
axis(side=2, at=c(0,seq(250,up.sup - 250,by=250)))
#axis(side=4, at=c(0,seq(250,up.sup + 250,by=250)))

# Support and density of the side-wise pdf:
supp1 <- seq(0,up.sup,by=5)
dens1 <- 25000*dlnorm(supp1,meanlog=a/(1-cc), sdlog=sqrt(statio.var))
points(c(stop.plot,stop.plot), c(0,up.sup), type="l", lty=1)
points(stop.plot+dens1, supp1, type="l", lwd=1, col="black")
# Adding the deterministic trajectory
#points(1:stop.plot,trial.det[1:stop.plot], type="l", lty=1,  col="black", lwd=2)
#text(stop.plot+10,K-50, expression(K[b]), cex=1.15, col = "black")
#points(stop.plot,K,pch=16, col="black", cex=0.8)
# Adding the trajectory of the Expected Value
points(1:stop.plot,E.trial.mean[1:stop.plot], type="l", lty=1,  col="blue", lwd=2)
points(stop.plot,exp((a/(1-cc)) + statio.var/2 ), pch=16, cex=0.8,col="blue")
text(stop.plot+12,50+exp((a/(1-cc)) + statio.var/2),
     expression(paste("E[",N[t],"]", sep="")), cex=1.15, col = "blue")

# Adding polygon with p(extinction)
nc <- 300
points(c(1,len),c(nc,nc), type="l", lty=2, lwd=1,col="red")
# Find index of intersection of pdf and cuttoff point nc:
int.ind <- which.min(((supp1)-nc)^2)
#dens1.intersect <- round((stop.plot+dens1)[int.ind],digits=2)
xs4polygon <- stop.plot+dens1[1:int.ind]   
ys.low4polygon <- supp1[1:int.ind]
ys.hi4polygon <- rep(nc,int.ind)
polygon(x=c(xs4polygon,rev(xs4polygon)), y=c(ys.low4polygon,ys.hi4polygon),
        col=scales::alpha("red",0.2), border = NA)
text(stop.plot-5, nc-100, expression(n[c]),cex=1.15, col = "red")

mtext("Time", side=1, outer=TRUE, cex=1.15, line=1)
mtext("Population abundance", side=2, outer=TRUE, cex=1.15, line=1)
