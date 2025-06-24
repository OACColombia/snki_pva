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