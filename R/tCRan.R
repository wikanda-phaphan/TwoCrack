#' @title The Two-parameter Crack Distribution
#'
#' @description The two-parameter crack distribution is a positive skewness model, which is extensively used to model
#' failure times of fatiguing materials. This distribution proposed by Saengthong and Bodhisuwan (2014).
#'
#' @param X vector of data
#' @param n a number of observations
#' @param lambda a value of the parameter lambda
#' @param theta a value of the parameter theta
#'
#' @return ?
#'
#' @examples
#' lambda <- 2
#' theta <- 2
#' X <- c(7,-2)
#' dIG(X,lambda,theta)
#' dLBIG(X,lambda,theta)
#' dBS(X,theta,lambda)
#' dTwoCrack(X,theta,lambda)
#'
###############Pack TwoCrack distribution#####################
##########################dTG################################
dIG=function(X,lambda,theta){
  nn<-length(X)
  p<-rep(0,nn)
  for (i in 1:nn) {

    if (X[c(i)]>0) {
      p[c(i)]<-(lambda/(theta*sqrt(2*pi)))*((theta/X[c(i)])^(3/2))*exp(-0.5*(sqrt(X[c(i)]/theta)-(lambda*sqrt(theta/X[c(i)]))))
    } else {
      p[c(i)]<-0
    }
  }
  return(p)
}
#dIG(X,lambda,theta)

##################################dLBIG####################
dLBIG=function(X,lambda,theta)
{
  nn<-length(X)
  p<-rep(0,nn)
  for (i in 1:nn) {

    if (X[c(i)]>0) {
      p[c(i)]<-(1/(theta*sqrt(2*pi)))*((theta/X[c(i)])^(1/2))*exp(-0.5*(sqrt(X[c(i)]/theta)-(lambda*sqrt(theta/X[c(i)]))))
    } else {
      p[c(i)]<-0
    }
  }
  return(p)
}
#dLBIG(X,lambda,theta)
#########################dBS########################
dBS=function(X,theta,lambda){
  nn<-length(X)
  p<-rep(0,nn)
  for (i in 1:nn) {

    if (X[c(i)]>0) {
      p[c(i)]<-0.5*dIG(X[c(i)],theta,lambda)+0.5*dLBIG(X[c(i)],theta,lambda)
    } else {
      p[c(i)]<-0
    }
  }
  return(p)
}
#dBS(X,theta,lambda)
############################the probability density function =dTwoCrack############################################
dTwoCrack=function(X,theta,lambda){
  nn<-length(X)
  p<-rep(0,nn)
  for (i in 1:nn) {

    if (X[c(i)]>0) {
      p[c(i)]<-(theta/(theta+1))*dIG(X[c(i)],theta,lambda)+(1/(theta+1))*dLBIG(X[c(i)],theta,lambda)
    } else {
      p[c(i)]<-0
    }
  }
  return(p)
}

#dTwoCrack(X,theta,lambda)
#########################the random generation=rTwoCrack############################################
fM=function(X,lambda,theta){
  dTwoCrack(X,theta,lambda)/dBS(X,theta,lambda)
}
#fM(X,lambda,theta)

### BS-random numbers generation procedure=rBS ### 1 time
rBS=function(n,lambda,theta)
{
  alpha <- array(0,dim=c(n,1))
  BS <- array(0,dim=c(n,1))

  for (i in 1:n){
    alpha[i,1]<- rnorm(1,0,1)
    BS[i,1] <- theta*(((alpha[i,1]/2)+sqrt(((alpha[i,1]^2)/4)+lambda))^2)
  }
  return(BS)
}
#rBS(n,lambda,theta)

### TC-random numbers generation procedure=rTwoCrack ### 1 time
rTwoCrack=function(n,lambda,theta)
{
  M <-2*max((theta/(theta+1)),(1/(theta+1)))
  X=NULL
  while(length(X)<n){
    y=rBS(n*M,lambda,theta)
    u=runif(n*M,0,M)
    X=c(X,y[u<fM(y,lambda,theta)])
  }
  X=X[1:n]
  MGG<-matrix(X, ncol=1)
  return(MGG)
}

#X<-rTwoCrack(n,lambda,theta)
#####################the distribution function =pTwoCrack#########################
Xx=function(X,lambda,theta)
{
  sqrt(X/theta)+(lambda*sqrt(theta/X))
}
#Xx<-Xx(X,lambda,theta)
Xy=function(X,lambda,theta)
{
  sqrt(X/theta)-(lambda*sqrt(theta/X))
}
#Xy<-Xy(X,lambda,theta)
pTwoCrack=function(X,lambda,theta)
{
  nn<-length(X)
  p<-rep(0,nn)
  for (i in 1:nn) {

    if (X[c(i)]>0) {
      Xx<-Xx(X[c(i)],lambda,theta)
      Xy<-Xy(X[c(i)],lambda,theta)
      p[c(i)]<-dnorm(Xy,mean=0,sd=1)-((1-theta)/(theta+1))*exp(2*lambda*(1-dnorm(Xx,mean=0,sd=1)))
    } else {
      p[c(i)]<-0
    }
  }
  return(p)

}
#pTwoCrack(X,lambda,theta)


### IG-random numbers generation procedure ### 1 time = rIG

rIG=function(n,lambda,theta)
{
  alpha1 <- array(0,dim=c(n,1))
  alpha <- array(0,dim=c(n,1))
  IG <- array(0,dim=c(n,1))
  u <- array(0,dim=c(n,1))
  jib <- array(0,dim=c(n,1))
  for (i in 1:n){
    alpha1[i,1]<- runif(1,0,1)
    alpha[i,1]<- rnorm(1,0,1)
    u[i,1] <- (lambda*theta)+(theta/2)*((alpha[i,1]^2)-sqrt((alpha[i,1]^4)+(4*lambda*(alpha[i,1]^2))))
    jib[i,1] <- (lambda*theta)/((lambda*theta)+u[i,1])
    if (alpha1[i,1] < jib[i,1]) {
      IG[i,1]=u[i,1]
    } else {
      IG[i,1]=(lambda^2)*(theta^2)/u[i,1]
    }
  }
  #return(list(IG,u,alpha))
  return(IG)
}
#X<-rIG(n,lambda,theta)

### LBIG-random numbers generation procedure ### 1 time = rLBIG
rLBIG=function(n,lambda,theta)
{
  alpha1 <- array(0,dim=c(n,1))
  alpha <- array(0,dim=c(n,1))
  alpha0 <- array(0,dim=c(n,1))
  LBIG <- array(0,dim=c(n,1))
  IG <- array(0,dim=c(n,1))
  u <- array(0,dim=c(n,1))
  jib <- array(0,dim=c(n,1))
  for (i in 1:n){
    alpha1[i,1]<- runif(1,0,1)
    alpha[i,1]<- rnorm(1,0,1)
    alpha0[i,1]<- rnorm(1,0,1)
    u[i,1] <- (lambda*theta)+(theta/2)*((alpha[i,1]^2)-sqrt((alpha[i,1]^4)+(4*lambda*(alpha[i,1]^2))))
    jib[i,1] <- (lambda*theta)/((lambda*theta)+u[i,1])
    if (alpha1[i,1] < jib[i,1]) {
      IG[i,1]=u[i,1]
    } else {
      IG[i,1]=(lambda^2)*(theta^2)/u[i,1]
    }
    LBIG[i,1]=IG[i,1]+(theta*(alpha0[i,1]^2))
  }
  #return(list(LBIG, IG,u,alpha))
  return(LBIG)
}
#X<-rLBIG(n,lambda,theta)
