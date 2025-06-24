# Running the GSSM in Nimble 


mu <- 5
theta <- 0.23
sigsq <- 0.09

cc <- exp(-theta)
a <- mu*(1-cc)
Vinf <- sigsq/(1-cc^2)

len <- 30
Onesim <- Pois.gssm(pars=c(mu=mu, theta=theta, sigsq=sigsq), len=len)
plot(1:len, Onesim[,2], pch=16, type="b")


ncl <- 10
data <- list(N=nimMatrix(Onesim[,2], nrow= nrow(Onesim), ncol=ncl))
constants <- list(qp1 = nrow(Onesim), ncl=ncl)
monitors <- c("mu", "theta", "sigmasq", "a", "cc", "car.cap", "Vinf")

tmp.model <- nimbleMCMC( code=GSSM, data = data, 
                         constants=constants, #inits = inits,
                         monitors = monitors, thin = 10,
                         niter = 50000, nburnin = 10000, nchains = 3,
                         summary = TRUE, samplesAsCodaMCMC = TRUE)
#tmp.model$summary$all.chains

cbind(tmp.model$summary$all.chains[,1], c(
Vinf,a,car.cap=exp(mu),cc,mu,sigsq, theta))


