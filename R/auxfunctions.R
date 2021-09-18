#library(purrr)
#library(mvtnorm)
#library(MASS)
#library(dplyr)
#library(ggplot2)
#library(nlme)
#library(numDeriv)
#library(moments)
#source("InformationMatrix-v3.R")
is.wholenumber <- function(x, tol1 = .Machine$double.eps^0.5)  abs(x - round(x)) < tol1
################################################################
#Root of a symmetric matrix
################################################################
errorVar<- function(times,object=NULL,sigma2=NULL,depStruct=NULL,phi=NULL) {
  if((!is.null(object))&(!inherits(object,c("SMSN","SMN")))) stop("object must inherit from class SMSN or SMN")
  if (is.null(object)&is.null(depStruct)) stop("object or depStruct must be provided")
  if (is.null(object)&is.null(sigma2)) stop("object or sigma2 must be provided")
  if (is.null(depStruct)) depStruct<-object$depStruct
  if (depStruct=="CI") depStruct = "UNC"
  if (depStruct!="UNC" & is.null(object)&is.null(phi)) stop("object or phi must be provided")
  if (!(depStruct %in% c("UNC","ARp","CS","DEC","CAR1"))) stop("accepted depStruct: UNC, ARp, CS, DEC or CAR1")
  if (is.null(sigma2)) sigma2<-object$estimates$sigma2
  if (is.null(phi)&depStruct!="UNC") phi<-object$estimates$phi
  if (depStruct=="ARp" & (any(!is.wholenumber(times))|any(times<=0))) stop("times must contain positive integer numbers when using ARp dependency")
  if (depStruct=="ARp" & any(tphitopi(phi)< -1|tphitopi(phi)>1)) stop("AR(p) non stationary, choose other phi")
  #
  if (depStruct=="UNC") var.out<- sigma2*diag(length(times))
  if (depStruct=="ARp") var.out<- sigma2*CovARp(phi,times)
  if (depStruct=="CS") var.out<- sigma2*CovCS(phi,length(times))
  if (depStruct=="DEC") var.out<- sigma2*CovDEC(phi[1],phi[2],times)
  if (depStruct=="CAR1") var.out<- sigma2*CovDEC(phi,1,times)
  var.out
}

matrix.sqrt <- function(A)
{
  if (length(A)==1) return(sqrt(A))
  else{
    sva <- svd(A)
    if (min(sva$d)>=0) {
      Asqrt <- sva$u%*%diag(sqrt(sva$d))%*%t(sva$v) # svd e decomposi??o espectral
      if (all(abs(Asqrt%*%Asqrt-A)<1e-4)) return(Asqrt)
      else stop("Matrix square root is not defined/not real")
    }
    else stop("Matrix square root is not defined/not real")
  }
}
################################################################
#trace of a matrix of dim >=1
################################################################
traceM <- function(Mat){
  if(length(Mat)==1) tr<- as.numeric(Mat)
  else tr<-sum(diag(Mat))
  tr
}
################################################################
#inverter ordem de hierarquia de uma lista com nomes
################################################################
revert_list <- function(ls) { # @Josh O'Brien
  # get sub-elements in same order
  x <- lapply(ls, `[`, names(ls[[1]]))
  # stack and reslice
  apply(do.call(rbind, x), 2, as.list)
}

#transformation function: pi to phi
estphit <- function(pit) {
  p <- length(pit)
  Phi <- matrix(0,ncol=p,nrow=p)
  if (p>1) {
    diag(Phi) <- pit
    for (j in 2:p) {
      for (k in 1:(j-1)) {
        Phi[j,k] <- Phi[j-1,k] - pit[j]*Phi[j-1,j-k]
      }
    }
    return(Phi[p,])
  }
  else return(pit)
}
#transformation function: phi to pi
tphitopi <- function(phit) {
  p <- length(phit)
  Phi <- matrix(0,ncol=p,nrow=p)
  Phi[p,] <- phit
  if (p>1) {
    for (k in p:2) {
      for (i in 1:(k-1)) {
        Phi[k-1,i] <- (Phi[k,i] + Phi[k,k]*Phi[k,k-i])/(1-Phi[k,k]^2)
      }
    }
    return(diag(Phi))
  }
  else return(phit)
}

# matriz Mn = 1/sigma2 * autocov(arp) in function of vector of times
CovARp<-function(phi,ti) {
  p <- length(phi)
  n <- max(ti)
  if (n==1) Rn <- matrix(1)
  else Rn<- toeplitz(ARMAacf(ar=phi, ma=0, lag.max = n-1))
  rhos <- ARMAacf(ar=phi, ma=0, lag.max = p)[-1]
  Rn<- Rn/(1-sum(rhos*phi))
  return(as.matrix(Rn[ti,ti]))
}

# corr matriz of comp symmetry
CovCS<-function(phi,n) {
  if (n==1) Rn <- matrix(1)
  else Rn<- toeplitz(c(1,rep(phi,n-1)))
  return(Rn)
}

# corr matriz of CAR1
# CovDEC(phi,theta=1)

# corr matriz of DEC
CovDEC<-function(phi1,phi2,ti) {
  ni <- length(ti)
  Rn <- diag(ni)
  if (ni==1) Rn <- matrix(1)
  else {
    for (i in 1:(ni-1)) for (j in (i+1):ni) Rn[i,j] <- phi1^(abs(ti[i]-ti[j])^phi2)
    Rn[lower.tri(Rn)] <-  t(Rn)[lower.tri(Rn)]
  }
  return(Rn)
}

Dmatrix <- function(dd) {
  q2 <- length(dd)
  q1 <- -.5+sqrt(1+8*q2)/2
  if (q1%%1 != 0) stop("wrong dimension of dd")
  D1 <- matrix(nrow = q1,ncol = q1)
  D1[upper.tri(D1,diag = T)] <- as.numeric(dd)
  D1[lower.tri(D1)] <- t(D1)[lower.tri(D1)]
  return(D1)
}

# gerando smsn para um individuo usando rep hierarquica
gerar_ind_smsn = function(ni,Sig,Di,beta,lambda,distr="sn",nu=NULL) {
  if (distr=="sn") {ui=1; c.=-sqrt(2/pi)}
  if (distr=="st") {ui=rgamma(1,nu/2,nu/2)
                    c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {ui=rbeta(1,nu,1)
                    c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {ui=ifelse(runif(1)<nu[1],nu[2],1)
                    c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%(lambda)))
  Delta = matrix.sqrt(Di)%*%delta
  Gammab = Di - Delta%*%t(Delta)
  Xi = cbind(1,runif(ni,0,2))
  Zi = matrix(1,nrow=ni)
  Beta = matrix(beta,ncol=1)
  q1 = nrow(Di)
  ti = c.+abs(rnorm(1,0,ui^-.5))
  bi = t(rmvnorm(1,Delta*ti,sigma=ui^(-1)*Gammab))
  Yi = t(rmvnorm(1,Xi%*%Beta+Zi%*%bi,sigma=ui^(-1)*Sig))
  return(data.frame(y=Yi,x=Xi[,2],tempo=1:ni,ui=ui))
}
################################################################
#Log-likelihood - AR(p)
################################################################
ljnormalAR <-function(j,y,x,z,time,beta1,Gammab,Deltab,sigmae,phiAR){
  c. = -sqrt(2/pi)
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovARp(phiAR,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))
}
#
ljtAR <-function(j,nu,y,x,z,time,beta1,Gammab,Deltab,sigmae,phiAR){
  c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovARp(phiAR,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med))
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  #log(2*dmvt(y1,delta = med, sigma = Psi, df = nu,log=F)*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))#veroST1(Psi,Ajj,dj,nu,pp=njj))
  log(2*dtj*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))
}
#
ljsAR <-function(j,nu,y,x,z,time,beta1,Gammab,Deltab,sigmae,phiAR){
  c.=-sqrt(2/pi)*nu/(nu-.5)
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovARp(phiAR,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med))
  #f <- function(u) u^(nu - 1)*dmvnorm(y1,med,Psi/u)*pnorm(u^(1/2)*Ajj)
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
  resp <- integrate(Vectorize(f2),0,1)$value
  log(2*nu*resp)
}
#
ljcnAR <-function(j,nu,y,x,z,time,beta1,Gammab,Deltab,sigmae,phiAR){
  c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovARp(phiAR,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(sqrt(nu[2])*Ajj,0,1)+
           (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
}

logveroAR = function(y,x,z,time,ind,beta1,sigmae,phiAR,D1,lambda,distr,nu){ #ind = indicadora de individuo

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda));
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  N <-length(ind)

  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalAR,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtAR,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsAR,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnAR,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiAR=phiAR))
  lv
}
logveroARpi = function(y,x,z,time,ind,beta1,sigmae,piAR,D1,lambda,distr,nu){ #ind = indicadora de individuo

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  phiAR <- estphit(piAR)
  N <-length(ind)

  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalAR,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtAR,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsAR,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnAR,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiAR=phiAR))
  lv
}

##############################################################################
# EM - AR(p)
##############################################################################
calcbi_emjAR <- function(jseq,y,x,z,time,beta1,Gammab, Deltab,sigmae,piAR,zeta,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  #
  y1=y[jseq]
  t1=time[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  nj = length(y1)
  Sigma = sigmae*CovARp(phi = estphit(piAR),t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ajj<-as.numeric(mutj/sqrt(Mtj2))
  D1<- Gammab+Deltab%*%t(Deltab)
  #
  mediab<-D1%*%t(z1)%*%solve(Psi)%*%(y1-med)+c.*Deltab
  Lambda<-solve(solve(D1)+t(z1)%*%solve(Sigma)%*%z1)
  #
  if  (distr=="sn"){
    bi<-mediab+Lambda%*%zeta/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))*as.numeric(dnorm(Ajj,0,1))/as.numeric(pnorm(Ajj,0,1))
  } else if (distr=="st"){
    dtj = gamma((nu+nj)/2)/gamma(nu/2)/pi^(nj/2)/sqrt(det(Psi))*nu^(-nj/2)*(dj/nu+1)^(-(nj+nu)/2)
    denST = 2*dtj*pt(sqrt(nu+nj)*Ajj/sqrt(dj+nu),nu+nj)
    esper2<- as.numeric(dmvt(y1,delta=med,sigma=Psi,df=nu,log=F)*gamma((nu+nj-1)/2)*(nu+dj)^((nu+nj)/2)/
                          (denST*pi^.5*gamma((nu+nj)/2)*(nu+dj+Ajj^2)^((nu+nj-1)/2)))
    bi<-mediab+Lambda%*%zeta/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))*esper2
  } else if (distr=="ss"){
    f2 <- function(u) u^(nu - 1)*((2*pi)^(-nj/2))*(u^(nj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
    denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
    esper2<-2^(nu)*nu*gamma(nu-.5+nj/2)*pgamma(1,nu-.5+nj/2,(dj+Ajj^2)/2)/
      (denSS*(dj+Ajj^2)^(nu-.5+nj/2)*pi^(nj/2+.5)*det(Psi)^.5)
    bi<-mediab+Lambda%*%zeta*esper2/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))
  } else if (distr=="scn"){
    fy<-as.numeric(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                        (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
    esper2<-as.numeric(2*(nu[1]*nu[2]^(-1/2)*dmvnorm(y1,med,Psi/nu[2])*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                            (1-nu[1])*dmvnorm(y1,med,Psi)*dnorm(Ajj,0,1))/fy)
    bi<-mediab+Lambda%*%zeta*esper2/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))
  }
  bi
}
# calcui_emjAR = function(jseq, y, x, z,time, beta1, Gammab, Deltab, sigmae,piAR,distr,nu) {
#   if  (distr=="sn"){
#     return(1)
#   }
#   if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
#   if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
#   if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
#   #
#   y1=y[jseq]
#   t1=time[jseq]
#   p= ncol(x);q1=ncol(z)
#   x1=matrix(x[jseq,  ],ncol=p)
#   z1=matrix(z[jseq,  ],ncol=q1)
#   med<-x1%*%beta1 + c.*z1%*%Deltab
#   nj = length(y1)
#   Sigma = sigmae*CovARp(phi = estphit(piAR),t1)
#   Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
#   dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
#   Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
#   mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
#   Ajj<-as.numeric(mutj/sqrt(Mtj2))
#   D1<- Gammab+Deltab%*%t(Deltab)
#   #
#   if (distr=="st"){
#     dtj = gamma((nu+nj)/2)/gamma(nu/2)/pi^(nj/2)/sqrt(det(Psi))*nu^(-nj/2)*(dj/nu+1)^(-(nj+nu)/2)
#     denST = 2*dtj*pt(sqrt(nu+nj)*Ajj/sqrt(dj+nu),nu+nj)
#     uj<-as.numeric(2^2*(nu+dj)^(-(nu+nj+2)/2)*gamma((nj+nu+2)/2)*nu^(nu/2)*
#                      pt(sqrt((nj+nu+2)/(dj+nu))*Ajj,nj+nu+2)/(gamma(nu/2)*det(Psi)^.5*denST*pi^(nj/2)))
#   } else if (distr=="ss"){
#     f2esp <- function(s) pnorm(s^.5*Ajj)*dgamma(s,nu+1+nj/2,dj/2)/pgamma(1,nu+1+nj/2,dj/2)
#     EspVal <- integrate(f2esp,0,1)$value#mean(pnorm(S^(1/2)*Ajj))#
#     f2 <- function(u) u^(nu - 1)*((2*pi)^(-nj/2))*(u^(nj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
#     denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
#     uj<-2^(2+nu)*nu*gamma(nu+1+nj/2)*pgamma(1,nu+1+nj/2,dj/2)*EspVal/
#       (denSS*dj^(nu+1+nj/2)*pi^(nj/2)*det(Psi)^.5)
#   } else if (distr=="scn"){
#     fy<-as.numeric(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
#                         (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
#     uj<-2*(nu[1]*nu[2]*dmvnorm(y1,med,Psi/nu[2])*pnorm(nu[2]^(1/2)*Ajj,0,1)+
#              (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))/fy
#   }
#   return(uj)
# }
calc_ui = function(jseq, y, x, z,time=NULL, beta1, Gammab, Deltab, sigmae,phi=NULL,distr,
                   depStruct,nu) {
  if  (distr=="sn"){
    return(1)
  }
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  #
  y1=y[jseq]
  nj = length(y1)
  if (is.null(time)) {
    t1=seq_len(nj)
  } else t1=time[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  Sigma = errorVar(times=t1,sigma2=sigmae,depStruct=depStruct,phi=phi)
  #sigmae*CovARp(phi = estphit(piAR),t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ajj<-as.numeric(mutj/sqrt(Mtj2))
  D1<- Gammab+Deltab%*%t(Deltab)
  #
  if (distr=="st"){
    dtj = gamma((nu+nj)/2)/gamma(nu/2)/pi^(nj/2)/sqrt(det(Psi))*nu^(-nj/2)*(dj/nu+1)^(-(nj+nu)/2)
    denST = 2*dtj*pt(sqrt(nu+nj)*Ajj/sqrt(dj+nu),nu+nj)
    uj<-as.numeric(2^2*(nu+dj)^(-(nu+nj+2)/2)*gamma((nj+nu+2)/2)*nu^(nu/2)*
                     pt(sqrt((nj+nu+2)/(dj+nu))*Ajj,nj+nu+2)/(gamma(nu/2)*det(Psi)^.5*denST*pi^(nj/2)))
  } else if (distr=="ss"){
    f2esp <- function(s) pnorm(s^.5*Ajj)*dgamma(s,nu+1+nj/2,dj/2)/pgamma(1,nu+1+nj/2,dj/2)
    EspVal <- integrate(f2esp,0,1)$value#mean(pnorm(S^(1/2)*Ajj))#
    f2 <- function(u) u^(nu - 1)*((2*pi)^(-nj/2))*(u^(nj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
    denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
    uj<-2^(2+nu)*nu*gamma(nu+1+nj/2)*pgamma(1,nu+1+nj/2,dj/2)*EspVal/
      (denSS*dj^(nu+1+nj/2)*pi^(nj/2)*det(Psi)^.5)
  } else if (distr=="scn"){
    fy<-as.numeric(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                        (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
    uj<-2*(nu[1]*nu[2]*dmvnorm(y1,med,Psi/nu[2])*pnorm(nu[2]^(1/2)*Ajj,0,1)+
             (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))/fy
  }
  return(uj)
}
calcs_ui = function(jseq, y, x, z,time=NULL, beta1, D1, sigmae,phi=NULL,distr,
                   depStruct,nu) {
  if  (distr=="sn"){
    uj<-1
    return(uj)
  }
  y1=y[jseq]
  nj = length(y1)
  if (is.null(time)) {
    t1=seq_len(nj)
  } else t1=time[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1
  Sigma = errorVar(times=t1,sigma2=sigmae,depStruct=depStruct,phi=phi)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  #
  if (distr=="st"){
    uj<-(nj+nu)/(dj+nu)
  } else if (distr=="ss"){
    uj<-pgamma(1,nj/2+nu+1,dj/2)/pgamma(1,nj/2+nu,dj/2)*(nj+2*nu)/dj
  } else if (distr=="scn"){
    fy<-as.numeric((nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))
    uj<-as.numeric((nu[1]*nu[2]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))/fy
  }
  return(uj)
}

emjAR = function(jseq, y, x, z,time, beta1, Gammab, Deltab, sigmae,piAR,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  #
  y1=y[jseq]
  t1=time[jseq]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  nj = length(y1)
  Sigma = sigmae*CovARp(phi = estphit(piAR),t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ajj<-as.numeric(mutj/sqrt(Mtj2))
  D1<- Gammab+Deltab%*%t(Deltab)
  #
  if  (distr=="sn"){
    uj<-1
    esper<-as.numeric(dnorm(Ajj,0,1)/pnorm(Ajj,0,1))
  }

  if (distr=="st"){
    dtj = gamma((nu+nj)/2)/gamma(nu/2)/pi^(nj/2)/sqrt(det(Psi))*nu^(-nj/2)*(dj/nu+1)^(-(nj+nu)/2)
    denST = 2*dtj*pt(sqrt(nu+nj)*Ajj/sqrt(dj+nu),nu+nj)
    uj<-as.numeric(2^2*(nu+dj)^(-(nu+nj+2)/2)*gamma((nj+nu+2)/2)*nu^(nu/2)*
                     pt(sqrt((nj+nu+2)/(dj+nu))*Ajj,nj+nu+2)/(gamma(nu/2)*det(Psi)^.5*denST*pi^(nj/2)))
    esper <-as.numeric(2*nu^(nu/2)*gamma((nj+nu+1)/2)/(denST*pi^((nj+1)/2)*gamma(nu/2)*det(Psi)^.5*(nu+dj+Ajj^2)^((nu+nj+1)/2)))
  }

  if (distr=="ss"){
    f2esp <- function(s) pnorm(s^.5*Ajj)*dgamma(s,nu+1+nj/2,dj/2)/pgamma(1,nu+1+nj/2,dj/2)
    EspVal <- integrate(f2esp,0,1)$value#mean(pnorm(S^(1/2)*Ajj))#
    f2 <- function(u) u^(nu - 1)*((2*pi)^(-nj/2))*(u^(nj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
    denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
    uj<-2^(2+nu)*nu*gamma(nu+1+nj/2)*pgamma(1,nu+1+nj/2,dj/2)*EspVal/
      (denSS*dj^(nu+1+nj/2)*pi^(nj/2)*det(Psi)^.5)
    esper <- 2^(1+nu)*nu*gamma(nu+.5+nj/2)*pgamma(1,nu+.5+nj/2,(dj+Ajj^2)/2)/
      (denSS*(dj+Ajj^2)^(nu+.5+nj/2)*pi^(nj/2+.5)*det(Psi)^.5)
  }

  if (distr=="scn"){
    fy<-as.numeric(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                        (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
    uj<-2*(nu[1]*nu[2]*dmvnorm(y1,med,Psi/nu[2])*pnorm(nu[2]^(1/2)*Ajj,0,1)+
             (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))/fy
    esper<-as.numeric(2*(nu[1]*nu[2]^(1/2)*dmvnorm(y1,med,Psi/nu[2])*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                           (1-nu[1])*dmvnorm(y1,med,Psi)*dnorm(Ajj,0,1))/fy)
  }
  sSigma = solve(Sigma)
  sRi = sSigma*sigmae
  Tbj<-solve(solve(Gammab)+t(z1)%*%sSigma%*%z1)
  r<-Tbj%*%t(z1)%*%sSigma%*%(y1-x1%*%beta1)
  s1<-(diag(q1)-Tbj%*%t(z1)%*%sSigma%*%z1)%*%Deltab
  utj<-as.numeric(uj*(mutj+c.)+sqrt(Mtj2)*esper)
  ut2j<-as.numeric(uj*(mutj+c.)^2+Mtj2+sqrt(Mtj2)*(mutj+2*c.)*esper)
  ub<-uj*r+s1*utj
  utbj<- r*utj+s1*ut2j
  ub2j<-Tbj+uj*r%*%t(r)+s1%*%t(r)*utj+r%*%t(s1)*utj+s1%*%t(s1)*ut2j
  #
  sum1<-uj*t(x1)%*%sRi%*%x1 #denom beta
  sum2<-(t(x1)%*%sRi%*%(uj*y1-z1%*%ub)) #num beta
  sum3<-uj*t(y1-x1%*%beta1)%*%sRi%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%sRi%*%z1%*%ub-
    t(ub)%*%t(z1)%*%sRi%*%(y1-x1%*%beta1)+traceM(sRi%*%z1%*%ub2j%*%t(z1)) #soma do sig2
  sum4<-ub2j-utbj%*%t(Deltab)-Deltab%*%t(utbj)+
    ut2j*Deltab%*%t(Deltab) #soma do Gamma
  sum5<-utbj #num do delta
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,sum5=sum5,ut2j=ut2j,
                 uj=uj,ubj=ub,ub2j=ub2j)
  #if (calcbi) obj.out$bi=bi
  return(obj.out)
}

#função para maximizar em pi
lcAR <- function(piAR,beta1,sigmae,y,x,z,time,ind,u,ub,ub2) {
  m<-n_distinct(ind)
  N<-length(ind)
  indlevels <- levels(ind)
  soma <-0
  for (i in seq_len(m)) {
    jseq <- which(ind==indlevels[i])
    y1=y[jseq]
    t1=time[jseq]
    p= ncol(x)
    q1=ncol(z)
    x1=matrix(x[jseq,  ],ncol=p)
    z1=matrix(z[jseq,  ],ncol=q1)
    med = x1%*%beta1
    nj = length(y1)
    Sigma = CovARp(phi = estphit(piAR),t1)*sigmae
    sSigma = solve(Sigma)
    indi = which(names(u)==indlevels[i])
    uj = u[[indi]]
    ubj = ub[[indi]]
    ub2j = ub2[[indi]]
    soma = soma + as.numeric(-.5*log(det(Sigma))-.5*uj*t(y1-med)%*%sSigma%*%(y1-med)+
                               t(y1-med)%*%sSigma%*%z1%*%ubj-.5*traceM(sSigma%*%z1%*%ub2j%*%t(z1)))
  }
  soma
}

EM.SkewAR<- function(formFixed,formRandom,data,groupVar,pAR,timeVar,
                     distr,beta1,sigmae,phiAR,D1,lambda,nu,lb,lu,
                     precisao,informa,max.iter,showiter,showerroriter){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]

  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x)
  q1<-ncol(z)
  #
  if ((!is.null(phiAR)) & pAR!=length(phiAR)) stop("initial value from phi must be in agreement with pAR")
  if ((pAR%%1)!=0|pAR==0) stop("pAR must be integer greater than 1")

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  #zeta<-matrix.sqrt(solve(D1))%*%lambda

  if (is.null(phiAR)) {
    lmeAR = try(lme(formFixed,random=~1|ind,data=data,correlation=corARMA(p=pAR,q=0)),silent=T)
    if (class(lmeAR)=="try-error") piAR =as.numeric(pacf(y-x%*%beta1,lag.max=pAR,plot=F)$acf)
    else {
      phiAR = capture.output(lmeAR$modelStruct$corStruct)[3]
      phiAR = as.numeric(strsplit(phiAR, " ")[[1]])
      phiAR = phiAR[!is.na(phiAR)]
      piAR = tphitopi(phiAR)
    }
  } else piAR = tphitopi(phiAR)
  if (any(piAR< -1 | piAR>1)) stop("invalid initial value from phi")

  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,piAR,nu)

  criterio<-10
  count<-0
  llji = logveroARpi(y, x, z, time,ind, beta1, sigmae,piAR, D1, lambda, distr, nu)
  if (is.nan(llji)|is.infinite(abs(llji))) stop("NaN/infinity initial likelihood")

  while((criterio > precisao)&(count<max.iter)){
    #print(nu)

    count <- count + 1
    res_emj = revert_list(tapply(1:N,ind,emjAR,y=y, x=x, z=z,time=time, beta1=beta1, Gammab=Gammab,
                                 Deltab=Deltab, sigmae=sigmae,piAR=piAR, distr=distr,nu=nu))
    sum1 = Reduce("+",res_emj$sum1)
    sum2 = Reduce("+",res_emj$sum2)
    sum3 = sum(unlist(res_emj$sum3))
    sum4 = Reduce("+",res_emj$sum4)
    sum5 = Reduce("+",res_emj$sum5)
    ut2j = unlist(res_emj$ut2j,use.names = F)

    #if (calcbi) bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
    beta1<-solve(sum1)%*%sum2
    sigmae<-as.numeric(sum3)/N
    Gammab<-sum4/m
    Deltab<-sum5/sum(ut2j)
    #
    D1<-Gammab+Deltab%*%t(Deltab)
    lambda<-matrix.sqrt(solve(D1))%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%solve(D1)%*%Deltab))
    #zeta<-matrix.sqrt(solve(D1))%*%lambda
    #
    piAR<- optim(piAR,lcAR,gr = NULL,method = "L-BFGS-B", lower =rep(-.9999,pAR),
                 upper = rep(.9999,pAR),control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                 y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
    #
    logvero1<-function(nu){logveroARpi(y, x, z,time, ind, beta1, sigmae,piAR, D1, lambda, distr, nu)}

    if (distr=="sn"){ nu<-NULL} else
    {nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par}
    #
    param <- teta
    teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,piAR,nu)
    criterio2 <- as.numeric(sqrt((teta-param)%*%(teta-param)))
    llj1<-llji
    llji <- logveroARpi(y, x, z, time,ind, beta1, sigmae,piAR, D1, lambda, distr, nu)
    criterio <- abs((llji-llj1)/llj1)
    if (showiter&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio,"\r") #  criterium ",criterio," or ",criterio2,"\r")
  }
  if (count==max.iter) message("\n maximum number of iterations reachead")
  cat("\n")
  zeta<-matrix.sqrt(solve(D1))%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjAR,y=y, x=x, z=z, time=time, beta1=beta1, Gammab=Gammab,
                             Deltab=Deltab, sigmae=sigmae,piAR=piAR,zeta=zeta, distr=distr,nu=nu,simplify = FALSE)),ncol=q1,byrow = T)
  # bi = t(list.cbind(tapply(1:N,ind,calcbi_emjAR,y=y, x=x, z=z, time=time, beta1=beta1, Gammab=Gammab,
  #                         Deltab=Deltab, sigmae=sigmae,piAR=piAR,zeta=zeta, distr=distr,nu=nu,simplify = FALSE)))
  dd<-matrix.sqrt(D1)[upper.tri(D1, diag = T)]
  phiAR=estphit(piAR)
  theta = c(beta1,sigmae,phiAR,dd,lambda,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",paste0("phiAR",1:length(piAR)),paste0("Dsqrt",1:length(dd)),paste0("lambda",1:q1))
  else names(theta)<- c(colnames(x),"sigma2",paste0("phiAR",1:length(piAR)),paste0("Dsqrt",1:length(dd)),paste0("lambda",1:q1),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                            phi=phiAR,dsqrt=dd,D=D1,lambda=as.numeric(lambda)),
                  uhat=unlist(res_emj$uj)) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixAR(y,x,z,time,ind,beta1,sigmae,phiAR,D1,lambda,distr = distr,nu = nu),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      q2<-q1*(q1+1)/2
      desvios[(p+pAR+q2+2):(p+pAR+1+q2+q1)] <- rep(NA,q1)
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(llji)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=criterio
  obj.out
}

predictf.skewAR<- function(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr,pAR,theta){
  dataPred[,all.vars(formFixed)[1]] <- 0
  dataFit$ind <-dataFit[,groupVar]
  dataPred$ind <-dataPred[,groupVar]
  dataPred$ind <- droplevels(dataPred$ind)
  #
  #theta = beta1,sigmae,phiAR,D1,lambda,nu
  p <- ncol(model.matrix(formFixed,data=dataPred))
  q1 <- ncol(model.matrix(formRandom,data=dataPred))
  q2 <- q1*(q1+1)/2
  beta1 <- matrix(theta[1:p],ncol=1)
  sigmae <- as.numeric(theta[p+1])
  phiAR <- as.numeric(theta[(p+2):(p+pAR+1)])
  dd <- theta[(p+pAR+2):(p+pAR+1+q2)]
  lambda <- matrix(theta[(p+pAR+q2+2):(p+pAR+1+q2+q1)],ncol=1)
  if (distr=="sn") {nu<- NULL
                    c. = -sqrt(2/pi)}
  if (distr=="st") {nu<- theta[p+pAR+q2+q1+2]
                    c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {nu<- theta[p+pAR+q2+q1+2]
                    c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {nu<- theta[(p+pAR+q2+q1+2):(p+pAR+q2+q1+3)]
                    c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  if ((p+pAR+1+q2+q1+length(nu))!=length(theta)) stop("theta misspecified")
  D1sqrt <- Dmatrix(dd)
  D1 <- D1sqrt%*%D1sqrt
  #
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-D1sqrt%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  zeta<-matrix.sqrt(solve(D1))%*%lambda
  ypred <- numeric(length = nrow(dataPred))
  timepred <- numeric(length = nrow(dataPred))
  xpred<-matrix(nrow= nrow(dataPred),ncol=p)
  #
  for (indj in levels(dataPred$ind)) {
    #indj = levels(dataPred$ind)[1]
    dataFitj <- subset(dataFit,dataFit$ind==indj,select = c("ind",all.vars(formFixed),all.vars(formRandom),timeVar))
    dataPredj <- subset(dataPred,dataPred$ind==indj,select = c("ind",all.vars(formFixed),all.vars(formRandom),timeVar))
    if (!is.null(timeVar)) {
      dataFitj$time <- dataFitj[,timeVar]
      dataPredj$time <- dataPredj[,timeVar]
    }
    njFit = nrow(dataFitj)
    njPred = nrow(dataPredj)
    seqFit = 1:njFit
    seqPred = njFit+1:njPred
    #
    if (is.null(timeVar)) {
      dataFitj$time<- seqFit
      dataPredj$time<- seqPred
    }
    dataPlus <- rbind(dataFitj,dataPredj)
    #
    xPlus1 <- model.matrix(formFixed,data=dataPlus)
    zPlus1<-model.matrix(formRandom,data=dataPlus)
    z1 <- matrix(zPlus1[seqFit,],ncol=ncol(zPlus1))
    x1 <- matrix(xPlus1[seqFit,],ncol=ncol(xPlus1))
    z1Pred <- matrix(zPlus1[seqPred,],ncol=ncol(zPlus1))
    x1Pred <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
    #
    medFit <- x1%*%beta1 + c.*z1%*%Deltab
    medPred <- x1Pred%*%beta1 + c.*z1Pred%*%Deltab
    #
    y1=dataFitj[,all.vars(formFixed)[1]]
    SigmaPlus = sigmae*CovARp(phi = phiAR,c(dataFitj$time,dataPredj$time))
    PsiPlus<-(zPlus1)%*%(D1)%*%t(zPlus1)+SigmaPlus
    dj<-as.numeric(t(y1-medFit)%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit))
    LambdaPlus <- solve(solve(D1)+ t(zPlus1)%*%solve(SigmaPlus)%*%zPlus1)
    sPsiPlus <- solve(PsiPlus)
    lambdaBarPlus <- matrix.sqrt(sPsiPlus)%*%zPlus1%*%D1%*%zeta/as.numeric(sqrt(1+t(zeta)%*%LambdaPlus%*%zeta))
    vj <- matrix.sqrt(sPsiPlus)%*%lambdaBarPlus
    Psi22.1 <- PsiPlus[seqPred,seqPred]- PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]
    vjtil <- (vj[seqFit] + solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]%*%vj[seqPred])/
      as.numeric(sqrt(1+t(vj[seqPred])%*%Psi22.1%*%vj[seqPred]))
    Ajj<-as.numeric(t(vjtil)%*%(y1-medFit)) #as.numeric(mutj/sqrt(Mtj2))
    mu2.1 <- medPred + PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit)
    if (distr=="sn") tau1 <- dnorm(Ajj)/pnorm(Ajj)
    if (distr=="st") {
      dtj = gamma((nu+njFit)/2)/gamma(nu/2)/pi^(njFit/2)/sqrt(det(as.matrix(PsiPlus[seqFit,seqFit])))*nu^(-njFit/2)*
        (dj/nu+1)^(-(njFit+nu)/2)
      denST = 2*dtj*pt(sqrt(nu+njFit)*Ajj/sqrt(dj+nu),nu+njFit)
      tau1 <- as.numeric(dmvt(y1,delta=medFit,sigma=as.matrix(PsiPlus[seqFit,seqFit]),df=nu,log=F)*gamma((nu+njFit-1)/2)*(nu+dj)^((nu+njFit)/2)/
                           (denST*pi^.5*gamma((nu+njFit)/2)*(nu+dj+Ajj^2)^((nu+njFit-1)/2)))
    }
    if (distr=="ss") {
      f2 <- function(u) u^(nu - 1)*((2*pi)^(-njFit/2))*(u^(njFit/2))*((det(as.matrix(PsiPlus[seqFit,seqFit])))^(-1/2))*
        exp(-0.5*u*dj)*pnorm(u^(1/2)*Ajj)
      denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
      tau1<-2^(nu)*nu*gamma(nu-.5+njFit/2)*pgamma(1,nu-.5+njFit/2,(dj+Ajj^2)/2)/
        (denSS*(dj+Ajj^2)^(nu-.5+njFit/2)*pi^(njFit/2+.5)*det(as.matrix(PsiPlus[seqFit,seqFit]))^.5)
    }
    if (distr=="scn") {
      fy<-as.numeric(2*(nu[1]*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                          (1-nu[1])*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]))*pnorm(Ajj,0,1)))
      tau1<-as.numeric(2*(nu[1]*nu[2]^(-1/2)*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]/nu[2]))*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                            (1-nu[1])*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]))*dnorm(Ajj,0,1))/fy)
    }
    ypredj <- mu2.1 + tau1/as.numeric(sqrt(1+t(vj[seqPred])%*%Psi22.1%*%vj[seqPred]))*Psi22.1%*%vj[seqPred]
    ypred[dataPred$ind==indj] <- ypredj
    xpred[dataPred$ind==indj,] <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
    timepred[dataPred$ind==indj] <- dataPredj$time
  }
  xpred = as.data.frame(xpred)
  colnames(xpred) = colnames(xPlus1)
  if (all(xpred[,1]==1)) xpred=xpred[-1]
  if (is.null(timeVar)) dfout = data.frame(groupVar=dataPred$ind,xpred,ypred)
  else dfout = data.frame(groupVar=dataPred$ind,time=timepred,xpred,ypred)
  dfout
}

################################################################
#Log-likelihood - independent
################################################################
ljnormal <-function(j,y,x,z,beta1,Gammab,Deltab,sigmae){
  c. = -sqrt(2/pi)
  y1=y[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj) #z1 D1 z1^t + sig2*I_nj ??
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))
}
#
ljt <-function(j,nu,y,x,z,beta1,Gammab,Deltab,sigmae){
  c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  y1=y[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med))
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  #log(2*dmvt(y1,delta = med, sigma = Psi, df = nu,log=F)*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))#veroST1(Psi,Ajj,dj,nu,pp=njj))
  log(2*dtj*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))
}
#
ljs <-function(j,nu,y,x,z,beta1,Gammab,Deltab,sigmae){
  c.=-sqrt(2/pi)*nu/(nu-.5)
  y1=y[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med))
  #f <- function(u) u^(nu - 1)*dmvnorm(y1,med,Psi/u)*pnorm(u^(1/2)*Ajj)
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
  resp <- integrate(Vectorize(f2),0,1)$value
  log(2*nu*resp)
}
#
ljcn <-function(j,nu,y,x,z,beta1,Gammab,Deltab,sigmae){
  c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(sqrt(nu[2])*Ajj,0,1)+
           (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
}
logvero = function(y,x,z,ind,beta1,sigmae,D1,lambda,distr,nu){ #ind = indicadora de individuo
  m<-n_distinct(ind)
  N<-length(ind)
  p<-dim(x)[2]
  q1<-dim(z)[2]

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)

  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormal,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljt,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljs,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcn,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  lv
}

##############################################################################
# EM - independent
##############################################################################
calcbi_emj <- function(jseq,y,x,z,beta1,Gammab, Deltab,sigmae,zeta,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  #
  y1=y[jseq]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  nj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(nj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(nj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(nj)+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ajj<-as.numeric(mutj/sqrt(Mtj2))
  D1<- Gammab+Deltab%*%t(Deltab)
  #
  mediab<-D1%*%t(z1)%*%solve(Psi)%*%(y1-med)+c.*Deltab
  Lambda<-solve(solve(D1)+t(z1)%*%z1/sigmae)
  #
  if  (distr=="sn"){
    bi<-mediab+Lambda%*%zeta/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))*as.numeric(dnorm(Ajj,0,1))/as.numeric(pnorm(Ajj,0,1))
  }

  if (distr=="st"){
    dtj = gamma((nu+nj)/2)/gamma(nu/2)/pi^(nj/2)/sqrt(det(Psi))*nu^(-nj/2)*(dj/nu+1)^(-(nj+nu)/2)
    denST = 2*dtj*pt(sqrt(nu+nj)*Ajj/sqrt(dj+nu),nu+nj)
    esper2<- as.numeric(dmvt(y1,delta=med,sigma=Psi,df=nu,log=F)*gamma((nu+nj-1)/2)*(nu+dj)^((nu+nj)/2)/
                            (denST*pi^.5*gamma((nu+nj)/2)*(nu+dj+Ajj^2)^((nu+nj-1)/2)))
    bi<-mediab+Lambda%*%zeta/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))*esper2
  }

  if (distr=="ss"){
    f2 <- function(u) u^(nu - 1)*((2*pi)^(-nj/2))*(u^(nj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
    denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
    esper2<-2^(nu)*nu*gamma(nu-.5+nj/2)*pgamma(1,nu-.5+nj/2,(dj+Ajj^2)/2)/
        (denSS*(dj+Ajj^2)^(nu-.5+nj/2)*pi^(nj/2+.5)*det(Psi)^.5)
    bi<-mediab+Lambda%*%zeta*esper2/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))
  }

  if (distr=="scn"){
    fy<-as.numeric(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                        (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
    esper2<-as.numeric(2*(nu[1]*nu[2]^(-1/2)*dmvnorm(y1,med,Psi/nu[2])*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                              (1-nu[1])*dmvnorm(y1,med,Psi)*dnorm(Ajj,0,1))/fy)
    bi<-mediab+Lambda%*%zeta*esper2/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))
  }
  bi
}

emj = function(jseq, y, x, z, beta1, Gammab, Deltab, sigmae,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  #
  y1=y[jseq]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  nj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(nj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(nj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(nj)+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ajj<-as.numeric(mutj/sqrt(Mtj2))
  D1<- Gammab+Deltab%*%t(Deltab)
  #
  if  (distr=="sn"){
    uj<-1
    esper<-as.numeric(dnorm(Ajj,0,1)/pnorm(Ajj,0,1))
  }

  if (distr=="st"){
    dtj = gamma((nu+nj)/2)/gamma(nu/2)/pi^(nj/2)/sqrt(det(Psi))*nu^(-nj/2)*(dj/nu+1)^(-(nj+nu)/2)
    denST = 2*dtj*pt(sqrt(nu+nj)*Ajj/sqrt(dj+nu),nu+nj)
    uj<-as.numeric(2^2*(nu+dj)^(-(nu+nj+2)/2)*gamma((nj+nu+2)/2)*nu^(nu/2)*
                     pt(sqrt((nj+nu+2)/(dj+nu))*Ajj,nj+nu+2)/(gamma(nu/2)*det(Psi)^.5*denST*pi^(nj/2)))
    esper <-as.numeric(2*nu^(nu/2)*gamma((nj+nu+1)/2)/(denST*pi^((nj+1)/2)*gamma(nu/2)*det(Psi)^.5*(nu+dj+Ajj^2)^((nu+nj+1)/2)))
  }

  if (distr=="ss"){
    f2esp <- function(s) pnorm(s^.5*Ajj)*dgamma(s,nu+1+nj/2,dj/2)/pgamma(1,nu+1+nj/2,dj/2)
    EspVal <- integrate(f2esp,0,1)$value#mean(pnorm(S^(1/2)*Ajj))#
    f2 <- function(u) u^(nu - 1)*((2*pi)^(-nj/2))*(u^(nj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
    denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
    uj<-2^(2+nu)*nu*gamma(nu+1+nj/2)*pgamma(1,nu+1+nj/2,dj/2)*EspVal/
      (denSS*dj^(nu+1+nj/2)*pi^(nj/2)*det(Psi)^.5)
    esper <- 2^(1+nu)*nu*gamma(nu+.5+nj/2)*pgamma(1,nu+.5+nj/2,(dj+Ajj^2)/2)/
      (denSS*(dj+Ajj^2)^(nu+.5+nj/2)*pi^(nj/2+.5)*det(Psi)^.5)
  }

  if (distr=="scn"){
    fy<-as.numeric(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                        (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
    uj<-2*(nu[1]*nu[2]*dmvnorm(y1,med,Psi/nu[2])*pnorm(nu[2]^(1/2)*Ajj,0,1)+
             (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))/fy
    esper<-as.numeric(2*(nu[1]*nu[2]^(1/2)*dmvnorm(y1,med,Psi/nu[2])*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                           (1-nu[1])*dmvnorm(y1,med,Psi)*dnorm(Ajj,0,1))/fy)
  }

  Tbj<-solve(solve(Gammab)+t(z1)%*%z1/sigmae)
  r<-Tbj%*%t(z1)%*%(y1-x1%*%beta1)/sigmae
  s1<-(diag(q1)-Tbj%*%t(z1)%*%z1/sigmae)%*%Deltab
  utj<-as.numeric(uj*(mutj+c.)+sqrt(Mtj2)*esper)
  ut2j<-as.numeric(uj*(mutj+c.)^2+Mtj2+sqrt(Mtj2)*(mutj+2*c.)*esper)
  ub<-uj*r+s1*utj
  utbj<- r*utj+s1*ut2j
  ub2j<-Tbj+uj*r%*%t(r)+s1%*%t(r)*utj+r%*%t(s1)*utj+s1%*%t(s1)*ut2j
  #
  sum1<-uj*t(x1)%*%x1 #denom beta
  sum2<-(t(x1)%*%(uj*y1-z1%*%ub)) #num beta
  sum3<-uj*t(y1-x1%*%beta1)%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%z1%*%ub-
    t(ub)%*%t(z1)%*%(y1-x1%*%beta1)+traceM(ub2j%*%t(z1)%*%z1) #soma do sig2
  sum4<-ub2j-utbj%*%t(Deltab)-Deltab%*%t(utbj)+
    ut2j*Deltab%*%t(Deltab) #soma do Gamma
  sum5<-utbj #num do delta
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,sum5=sum5,ut2j=ut2j,uj=uj)
  return(obj.out)
}

EM.Skew<- function(formFixed,formRandom,data,groupVar,distr,beta1,sigmae,D1,lambda,nu,lb,lu,
                   precisao,informa,max.iter,showiter,showerroriter){
  ti = Sys.time()
  x <- model.matrix(formFixed,data=data)
  varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  #
  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x)
  q1<-ncol(z)
  #
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  #zeta<-matrix.sqrt(solve(D1))%*%lambda

  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,nu)

  criterio<-10
  count<-0
  llji = logvero(y, x, z, ind, beta1, sigmae, D1, lambda, distr, nu)
  if (is.nan(llji)|is.infinite(abs(llji))) stop("NaN/infinity initial likelihood")

  while((criterio > precisao)&(count<max.iter)){

    count <- count + 1
    res_emj = revert_list(tapply(1:N,ind,emj,y=y, x=x, z=z, beta1=beta1, Gammab=Gammab,
                                 Deltab=Deltab, sigmae=sigmae, distr=distr,nu=nu))
    sum1 = Reduce("+",res_emj$sum1)
    sum2 = Reduce("+",res_emj$sum2)
    sum3 = sum(unlist(res_emj$sum3))
    sum4 = Reduce("+",res_emj$sum4)
    sum5 = Reduce("+",res_emj$sum5)
    ut2j = unlist(res_emj$ut2j,use.names = F)
    uj = unlist(res_emj$uj,use.names = F)
    #if (calcbi) bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))

    beta1<-solve(sum1)%*%sum2
    sigmae<-as.numeric(sum3)/N
    Gammab<-sum4/m
    Deltab<-sum5/sum(ut2j)
    #
    D1<-Gammab+Deltab%*%t(Deltab)
    lambda<-matrix.sqrt(solve(D1))%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%solve(D1)%*%Deltab))
    #
    logvero1<-function(nu){logvero(y, x, z, ind, beta1, sigmae, D1, lambda, distr, nu)}

    if (distr=="sn"){ nu<-NULL} else
    {
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par
    }
    param <- teta
    teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,nu)
    criterio2 <- as.numeric(sqrt((teta-param)%*%(teta-param)))
    llj1 <- llji
    llji <- logvero(y, x, z, ind, beta1, sigmae, D1, lambda, distr, nu)
    criterio <- abs((llji-llj1)/llj1)
    if (showiter&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio,"\r") #  criterium ",criterio," or ",criterio2,"\r")
  }
  if (count==max.iter) message("\n maximum number of iterations reachead")
  cat("\n")
  zeta<-matrix.sqrt(solve(D1))%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emj,y=y, x=x, z=z, beta1=beta1, Gammab=Gammab,
                             Deltab=Deltab, sigmae=sigmae,zeta=zeta, distr=distr,nu=nu,simplify = FALSE)),ncol=q1,byrow = T)
  # bi = t(list.cbind(tapply(1:N,ind,calcbi_emj,y=y, x=x, z=z, beta1=beta1, Gammab=Gammab,
  #                         Deltab=Deltab, sigmae=sigmae,zeta=zeta, distr=distr,nu=nu,simplify = FALSE)))
  dd<-matrix.sqrt(D1)[upper.tri(D1, diag = T)]
  theta = c(beta1,sigmae,dd,lambda,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",paste0("Dsqrt",1:length(dd)),paste0("lambda",1:q1))
  else names(theta)<- c(colnames(x),"sigma2",paste0("Dsqrt",1:length(dd)),paste0("lambda",1:q1),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                      dsqrt=dd,D=D1,lambda=as.numeric(lambda)),
                  uhat=unlist(res_emj$uj))
  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(Infmatrix(y,x,z,ind,beta1,sigmae,D1,lambda,distr = distr,nu = nu),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      q2<-q1*(q1+1)/2
      desvios[(p+q2+2):(p+1+q2+q1)] <- rep(NA,q1)
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(llji)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=criterio
  #class(obj.out) <- "EM.Skew"
  obj.out
}

predictf.skew<- function(formFixed,formRandom,dataFit,dataPred,groupVar,distr,theta){
  dataPred[,all.vars(formFixed)[1]] <- 0
  dataFit$ind <-dataFit[,groupVar]
  dataPred$ind <-dataPred[,groupVar]
  #
  p <- ncol(model.matrix(formFixed,data=dataPred))
  q1 <- ncol(model.matrix(formRandom,data=dataPred))
  q2 <- q1*(q1+1)/2
  beta1 <- matrix(theta[1:p],ncol=1)
  sigmae <- as.numeric(theta[p+1])
  dd <- theta[(p+2):(p+1+q2)]
  lambda <- matrix(theta[(p+q2+2):(p+1+q2+q1)],ncol=1)
  if (distr=="sn") {nu<- NULL
                    c. = -sqrt(2/pi)}
  if (distr=="st") {nu<- theta[p+q2+q1+2]
                    c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {nu<- theta[p+q2+q1+2]
                    c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {nu<- theta[(p+q2+q1+2):(p+q2+q1+3)]
                    c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  if ((p+1+q2+q1+length(nu))!=length(theta)) stop("theta misspecified")
  D1sqrt <- Dmatrix(dd)
  D1 <- D1sqrt%*%D1sqrt
  #
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-D1sqrt%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  zeta<-matrix.sqrt(solve(D1))%*%lambda
  ypred <- numeric(length = nrow(dataPred))
  dataPred$ind <- droplevels(dataPred$ind)
  xpred<-matrix(nrow= nrow(dataPred),ncol=p)
  #
  for (indj in levels(dataPred$ind)) {
    #indj = levels(dataPred$ind)[1]
    dataFitj <- subset(dataFit,dataFit$ind==indj,select = c("ind",all.vars(formFixed),all.vars(formRandom)))
    dataPredj <- subset(dataPred,dataPred$ind==indj,select = c("ind",all.vars(formFixed),all.vars(formRandom)))
    dataPlus <- rbind(dataFitj,dataPredj)
    njFit = nrow(dataFitj)
    njPred = nrow(dataPredj)
    seqFit = 1:njFit
    seqPred = njFit+1:njPred
    #
    xPlus1 <- model.matrix(formFixed,data=dataPlus)
    zPlus1<-model.matrix(formRandom,data=dataPlus)
    z1 <- matrix(zPlus1[seqFit,],ncol=ncol(zPlus1))
    x1 <- matrix(xPlus1[seqFit,],ncol=ncol(xPlus1))
    z1Pred <- matrix(zPlus1[seqPred,],ncol=ncol(zPlus1))
    x1Pred <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
    #
    medFit <- x1%*%beta1 + c.*z1%*%Deltab
    medPred <- x1Pred%*%beta1 + c.*z1Pred%*%Deltab
    #
    y1=dataFitj[,all.vars(formFixed)[1]]
    SigmaPlus = sigmae*diag(njFit+njPred)
    PsiPlus<-(zPlus1)%*%(D1)%*%t(zPlus1)+SigmaPlus
    dj<-as.numeric(t(y1-medFit)%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit))
    LambdaPlus <- solve(solve(D1)+ t(zPlus1)%*%solve(SigmaPlus)%*%zPlus1)
    sPsiPlus <- solve(PsiPlus)
    lambdaBarPlus <- matrix.sqrt(sPsiPlus)%*%zPlus1%*%D1%*%zeta/as.numeric(sqrt(1+t(zeta)%*%LambdaPlus%*%zeta))
    vj <- matrix.sqrt(sPsiPlus)%*%lambdaBarPlus
    Psi22.1 <- PsiPlus[seqPred,seqPred]- PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]
    vjtil <- (vj[seqFit] + solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]%*%vj[seqPred])/as.numeric(sqrt(1+t(vj[seqPred])%*%Psi22.1%*%vj[seqPred]))
    Ajj<-as.numeric(t(vjtil)%*%(y1-medFit)) #as.numeric(mutj/sqrt(Mtj2))
    mu2.1 <- medPred + PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit)
    if (distr=="sn") tau1 <- dnorm(Ajj)/pnorm(Ajj)
    if (distr=="st") {
      dtj = gamma((nu+njFit)/2)/gamma(nu/2)/pi^(njFit/2)/sqrt(det(as.matrix(PsiPlus[seqFit,seqFit])))*nu^(-njFit/2)*
        (dj/nu+1)^(-(njFit+nu)/2)
      denST = 2*dtj*pt(sqrt(nu+njFit)*Ajj/sqrt(dj+nu),nu+njFit)
      tau1 <- as.numeric(dmvt(y1,delta=medFit,sigma=as.matrix(PsiPlus[seqFit,seqFit]),df=nu,log=F)*gamma((nu+njFit-1)/2)*(nu+dj)^((nu+njFit)/2)/
                           (denST*pi^.5*gamma((nu+njFit)/2)*(nu+dj+Ajj^2)^((nu+njFit-1)/2)))
    }
    if (distr=="ss") {
      f2 <- function(u) u^(nu - 1)*((2*pi)^(-njFit/2))*(u^(njFit/2))*((det(as.matrix(PsiPlus[seqFit,seqFit])))^(-1/2))*
        exp(-0.5*u*dj)*pnorm(u^(1/2)*Ajj)
      denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
      tau1<-2^(nu)*nu*gamma(nu-.5+njFit/2)*pgamma(1,nu-.5+njFit/2,(dj+Ajj^2)/2)/
        (denSS*(dj+Ajj^2)^(nu-.5+njFit/2)*pi^(njFit/2+.5)*det(as.matrix(PsiPlus[seqFit,seqFit]))^.5)
    }
    if (distr=="scn") {
      fy<-as.numeric(2*(nu[1]*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                          (1-nu[1])*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]))*pnorm(Ajj,0,1)))
      tau1<-as.numeric(2*(nu[1]*nu[2]^(-1/2)*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]/nu[2]))*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                            (1-nu[1])*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]))*dnorm(Ajj,0,1))/fy)
    }
    ypredj <- mu2.1 + tau1/as.numeric(sqrt(1+t(vj[seqPred])%*%Psi22.1%*%vj[seqPred]))*Psi22.1%*%vj[seqPred]
    ypred[dataPred$ind==indj] <- ypredj
    xpred[dataPred$ind==indj,] <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
  }
  xpred = as.data.frame(xpred)
  colnames(xpred) = colnames(xPlus1)
  if (all(xpred[,1]==1)) xpred=xpred[-1]
  data.frame(groupVar=dataPred$ind,xpred,ypred)
}

################################################################
#Log-likelihood - CS
################################################################
ljnormalCS <-function(j,y,x,z,beta1,Gammab,Deltab,sigmae,phiCS){
  c. = -sqrt(2/pi)
  y1=y[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovCS(phiCS,njj)#CovARp(phiAR,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))
}
#
ljtCS <-function(j,nu,y,x,z,beta1,Gammab,Deltab,sigmae,phiCS){
  c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  y1=y[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovCS(phiCS,njj)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med))
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  #log(2*dmvt(y1,delta = med, sigma = Psi, df = nu,log=F)*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))#veroST1(Psi,Ajj,dj,nu,pp=njj))
  log(2*dtj*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))
}
#
ljsCS <-function(j,nu,y,x,z,beta1,Gammab,Deltab,sigmae,phiCS){
  c.=-sqrt(2/pi)*nu/(nu-.5)
  y1=y[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovCS(phiCS,njj)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med))
  #f <- function(u) u^(nu - 1)*dmvnorm(y1,med,Psi/u)*pnorm(u^(1/2)*Ajj)
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
  resp <- integrate(Vectorize(f2),0,1)$value
  log(2*nu*resp)
}
#
ljcnCS <-function(j,nu,y,x,z,beta1,Gammab,Deltab,sigmae,phiCS){
  c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovCS(phiCS,njj)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(sqrt(nu[2])*Ajj,0,1)+
           (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
}

logveroCS = function(y,x,z,ind,beta1,sigmae,phiCS,D1,lambda,distr,nu){ #ind = indicadora de individuo

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  N <-length(ind)

  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalCS,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiCS=phiCS))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtCS,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiCS=phiCS))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsCS,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiCS=phiCS))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnCS,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiCS=phiCS))
  lv
}

##############################################################################
# EM - CS
##############################################################################
calcbi_emjCS <- function(jseq,y,x,z,beta1,Gammab, Deltab,sigmae,phiCS,zeta,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  #
  y1=y[jseq]
  #t1=time[jseq]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  nj = length(y1)
  Sigma = sigmae*CovCS(phiCS,nj)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ajj<-as.numeric(mutj/sqrt(Mtj2))
  D1<- Gammab+Deltab%*%t(Deltab)
  #
  mediab<-D1%*%t(z1)%*%solve(Psi)%*%(y1-med)+c.*Deltab
  Lambda<-solve(solve(D1)+t(z1)%*%solve(Sigma)%*%z1)
  #
  if  (distr=="sn"){
    bi<-mediab+Lambda%*%zeta/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))*as.numeric(dnorm(Ajj,0,1))/as.numeric(pnorm(Ajj,0,1))
  }

  if (distr=="st"){
    dtj = gamma((nu+nj)/2)/gamma(nu/2)/pi^(nj/2)/sqrt(det(Psi))*nu^(-nj/2)*(dj/nu+1)^(-(nj+nu)/2)
    denST = 2*dtj*pt(sqrt(nu+nj)*Ajj/sqrt(dj+nu),nu+nj)
    esper2<- as.numeric(dmvt(y1,delta=med,sigma=Psi,df=nu,log=F)*gamma((nu+nj-1)/2)*(nu+dj)^((nu+nj)/2)/
                          (denST*pi^.5*gamma((nu+nj)/2)*(nu+dj+Ajj^2)^((nu+nj-1)/2)))
    bi<-mediab+Lambda%*%zeta/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))*esper2
  }

  if (distr=="ss"){
    f2 <- function(u) u^(nu - 1)*((2*pi)^(-nj/2))*(u^(nj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
    denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
    esper2<-2^(nu)*nu*gamma(nu-.5+nj/2)*pgamma(1,nu-.5+nj/2,(dj+Ajj^2)/2)/
      (denSS*(dj+Ajj^2)^(nu-.5+nj/2)*pi^(nj/2+.5)*det(Psi)^.5)
    bi<-mediab+Lambda%*%zeta*esper2/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))
  }

  if (distr=="scn"){
    fy<-as.numeric(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                        (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
    esper2<-as.numeric(2*(nu[1]*nu[2]^(-1/2)*dmvnorm(y1,med,Psi/nu[2])*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                            (1-nu[1])*dmvnorm(y1,med,Psi)*dnorm(Ajj,0,1))/fy)
    bi<-mediab+Lambda%*%zeta*esper2/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))
  }
  bi
}
emjCS = function(jseq, y, x, z, beta1, Gammab, Deltab, sigmae,phiCS,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  #
  y1=y[jseq]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  nj = length(y1)
  Sigma = sigmae*CovCS(phiCS,nj)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ajj<-as.numeric(mutj/sqrt(Mtj2))
  D1<- Gammab+Deltab%*%t(Deltab)
  #
  if  (distr=="sn"){
    uj<-1
    esper<-as.numeric(dnorm(Ajj,0,1)/pnorm(Ajj,0,1))
  }

  if (distr=="st"){
    dtj = gamma((nu+nj)/2)/gamma(nu/2)/pi^(nj/2)/sqrt(det(Psi))*nu^(-nj/2)*(dj/nu+1)^(-(nj+nu)/2)
    denST = 2*dtj*pt(sqrt(nu+nj)*Ajj/sqrt(dj+nu),nu+nj)
    uj<-as.numeric(2^2*(nu+dj)^(-(nu+nj+2)/2)*gamma((nj+nu+2)/2)*nu^(nu/2)*
                     pt(sqrt((nj+nu+2)/(dj+nu))*Ajj,nj+nu+2)/(gamma(nu/2)*det(Psi)^.5*denST*pi^(nj/2)))
    esper <-as.numeric(2*nu^(nu/2)*gamma((nj+nu+1)/2)/(denST*pi^((nj+1)/2)*gamma(nu/2)*det(Psi)^.5*(nu+dj+Ajj^2)^((nu+nj+1)/2)))
  }

  if (distr=="ss"){
    f2esp <- function(s) pnorm(s^.5*Ajj)*dgamma(s,nu+1+nj/2,dj/2)/pgamma(1,nu+1+nj/2,dj/2)
    EspVal <- integrate(f2esp,0,1)$value#mean(pnorm(S^(1/2)*Ajj))#
    f2 <- function(u) u^(nu - 1)*((2*pi)^(-nj/2))*(u^(nj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
    denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
    uj<-2^(2+nu)*nu*gamma(nu+1+nj/2)*pgamma(1,nu+1+nj/2,dj/2)*EspVal/
      (denSS*dj^(nu+1+nj/2)*pi^(nj/2)*det(Psi)^.5)
    esper <- 2^(1+nu)*nu*gamma(nu+.5+nj/2)*pgamma(1,nu+.5+nj/2,(dj+Ajj^2)/2)/
      (denSS*(dj+Ajj^2)^(nu+.5+nj/2)*pi^(nj/2+.5)*det(Psi)^.5)
  }

  if (distr=="scn"){
    fy<-as.numeric(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                        (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
    uj<-2*(nu[1]*nu[2]*dmvnorm(y1,med,Psi/nu[2])*pnorm(nu[2]^(1/2)*Ajj,0,1)+
             (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))/fy
    esper<-as.numeric(2*(nu[1]*nu[2]^(1/2)*dmvnorm(y1,med,Psi/nu[2])*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                           (1-nu[1])*dmvnorm(y1,med,Psi)*dnorm(Ajj,0,1))/fy)
  }

  sSigma = solve(Sigma)
  sRi = sSigma*sigmae
  Tbj<-solve(solve(Gammab)+t(z1)%*%sSigma%*%z1)
  r<-Tbj%*%t(z1)%*%sSigma%*%(y1-x1%*%beta1)
  s1<-(diag(q1)-Tbj%*%t(z1)%*%sSigma%*%z1)%*%Deltab
  utj<-as.numeric(uj*(mutj+c.)+sqrt(Mtj2)*esper)
  ut2j<-as.numeric(uj*(mutj+c.)^2+Mtj2+sqrt(Mtj2)*(mutj+2*c.)*esper)
  ub<-uj*r+s1*utj
  utbj<- r*utj+s1*ut2j
  ub2j<-Tbj+uj*r%*%t(r)+s1%*%t(r)*utj+r%*%t(s1)*utj+s1%*%t(s1)*ut2j
  #
  sum1<-uj*t(x1)%*%sRi%*%x1 #denom beta
  sum2<-(t(x1)%*%sRi%*%(uj*y1-z1%*%ub)) #num beta
  sum3<-uj*t(y1-x1%*%beta1)%*%sRi%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%sRi%*%z1%*%ub-
    t(ub)%*%t(z1)%*%sRi%*%(y1-x1%*%beta1)+traceM(sRi%*%z1%*%ub2j%*%t(z1)) #soma do sig2
  sum4<-ub2j-utbj%*%t(Deltab)-Deltab%*%t(utbj)+
    ut2j*Deltab%*%t(Deltab) #soma do Gamma
  sum5<-utbj #num do delta
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,sum5=sum5,ut2j=ut2j,
                 uj=uj,ubj=ub,ub2j=ub2j)
  #if (calcbi) obj.out$bi=bi
  return(obj.out)
}

#função para maximizar em pi
lcCS <- function(phiCS,beta1,sigmae,y,x,z,ind,u,ub,ub2) {
  m<-n_distinct(ind)
  N<-length(ind)
  indlevels <- levels(ind)
  soma <-0
  for (i in seq_len(m)) {
    jseq <- which(ind==indlevels[i])
    y1=y[jseq]
    p= ncol(x)
    q1=ncol(z)
    x1=matrix(x[jseq,  ],ncol=p)
    z1=matrix(z[jseq,  ],ncol=q1)
    med = x1%*%beta1
    nj = length(y1)
    Sigma = CovCS(phiCS,nj)*sigmae
    sSigma = solve(Sigma)
    indi = which(names(u)==indlevels[i])
    uj = u[[indi]]
    ubj = ub[[indi]]
    ub2j = ub2[[indi]]
    soma = soma + as.numeric(-.5*log(det(Sigma))-.5*uj*t(y1-med)%*%sSigma%*%(y1-med)+
                               t(y1-med)%*%sSigma%*%z1%*%ubj-.5*traceM(sSigma%*%z1%*%ub2j%*%t(z1)))
  }
  soma
}

EM.SkewCS<- function(formFixed,formRandom,data,groupVar,
                     distr,beta1,sigmae,phiCS,D1,lambda,nu,lb,lu,
                     precisao,informa,max.iter,showiter,showerroriter){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  #
  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x)
  q1<-ncol(z)
  #
  if (!is.null(phiCS) && length(phiCS)!=1) stop ("initial value from phi must have length 1 or be NULL")
  if (!is.null(phiCS)) if (phiCS<=0 | phiCS>=1) stop("0<initialValue$phi<1 needed")

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  #zeta<-matrix.sqrt(solve(D1))%*%lambda

  if (is.null(phiCS)) {
    phiCS = abs(as.numeric(pacf(y-x%*%beta1,lag.max=1,plot=F)$acf))
  }

  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,phiCS,nu)

  criterio<-10
  count<-0
  llji = logveroCS(y, x, z,ind, beta1, sigmae,phiCS, D1, lambda, distr, nu)
  if (is.nan(llji)|is.infinite(abs(llji))) stop("NaN/infinity initial likelihood")

  while((criterio > precisao)&(count<max.iter)){
    #print(nu)

    count <- count + 1
    res_emj = revert_list(tapply(1:N,ind,emjCS,y=y, x=x, z=z, beta1=beta1, Gammab=Gammab,
                                 Deltab=Deltab, sigmae=sigmae,phiCS=phiCS, distr=distr,nu=nu))
    sum1 = Reduce("+",res_emj$sum1)
    sum2 = Reduce("+",res_emj$sum2)
    sum3 = sum(unlist(res_emj$sum3))
    sum4 = Reduce("+",res_emj$sum4)
    sum5 = Reduce("+",res_emj$sum5)
    ut2j = unlist(res_emj$ut2j,use.names = F)

    #if (calcbi) bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))

    beta1<-solve(sum1)%*%sum2
    sigmae<-as.numeric(sum3)/N
    Gammab<-sum4/m
    Deltab<-sum5/sum(ut2j)
    #
    D1<-Gammab+Deltab%*%t(Deltab)
    lambda<-matrix.sqrt(solve(D1))%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%solve(D1)%*%Deltab))
    #zeta<-matrix.sqrt(solve(D1))%*%lambda
    #
    phiCS <- optim(phiCS,lcCS,gr = NULL,method = "L-BFGS-B", lower =0,
                   upper = .9999,control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
    #
    logvero1<-function(nu){logveroCS(y, x, z, ind, beta1, sigmae,phiCS, D1, lambda, distr, nu)}

    if (distr=="sn"){ nu<-NULL} else
    {nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par}
    #
    param <- teta
    teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,phiCS,nu)
    criterio2 <- as.numeric(sqrt((teta-param)%*%(teta-param)))
    llj1<-llji
    llji <- logveroCS(y, x, z, ind, beta1, sigmae,phiCS, D1, lambda, distr, nu)
    criterio <- abs((llji-llj1)/llj1)
    if (showiter&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio,"\r") #  criterium ",criterio," or ",criterio2,"\r")
  }
  if (count==max.iter) message("\n maximum number of iterations reachead")
  cat("\n")
  zeta<-matrix.sqrt(solve(D1))%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjCS,y=y, x=x, z=z, beta1=beta1, Gammab=Gammab,
                             Deltab=Deltab, sigmae=sigmae,phiCS=phiCS,zeta=zeta, distr=distr,nu=nu,simplify = FALSE)),ncol=q1,byrow = T)
  # bi = t(list.cbind(tapply(1:N,ind,calcbi_emjCS,y=y, x=x, z=z, beta1=beta1, Gammab=Gammab,
  #                         Deltab=Deltab, sigmae=sigmae,phiCS=phiCS,zeta=zeta, distr=distr,nu=nu,simplify = FALSE)))
  dd<-matrix.sqrt(D1)[upper.tri(D1, diag = T)]
  theta = c(beta1,sigmae,phiCS,dd,lambda,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phiCS",paste0("Dsqrt",1:length(dd)),paste0("lambda",1:q1))
  else names(theta)<- c(colnames(x),"sigma2","phiCS",paste0("Dsqrt",1:length(dd)),paste0("lambda",1:q1),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                phi=phiCS,dsqrt=dd,D=D1,lambda=as.numeric(lambda)),
                  uhat=unlist(res_emj$uj))
  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixCS(y,x,z,ind,beta1,sigmae,phiCS,D1,lambda,distr = distr,nu = nu),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      q2<-q1*(q1+1)/2
      desvios[(p+q2+3):(p+2+q2+q1)] <- rep(NA,q1)
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(llji)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=criterio
  obj.out
}

predictf.skewCS<- function(formFixed,formRandom,dataFit,dataPred,groupVar,distr,theta){
  dataPred[,all.vars(formFixed)[1]] <- 0
  dataFit$ind <-dataFit[,groupVar]
  dataPred$ind <-dataPred[,groupVar]
  dataPred$ind <- droplevels(dataPred$ind)
  #
  #theta = beta1,sigmae,phiAR,D1,lambda,nu
  p <- ncol(model.matrix(formFixed,data=dataPred))
  q1 <- ncol(model.matrix(formRandom,data=dataPred))
  q2 <- q1*(q1+1)/2
  beta1 <- matrix(theta[1:p],ncol=1)
  sigmae <- as.numeric(theta[p+1])
  phiCS <- as.numeric(theta[(p+2)])
  dd <- theta[(p+3):(p+2+q2)]
  lambda <- matrix(theta[(p+q2+3):(p+2+q2+q1)],ncol=1)
  if (distr=="sn") {nu<- NULL
                    c. = -sqrt(2/pi)}
  if (distr=="st") {nu<- theta[p+q2+q1+3]
                    c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {nu<- theta[p+q2+q1+3]
                    c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {nu<- theta[(p+q2+q1+3):(p+q2+q1+4)]
                    c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  if ((p+2+q2+q1+length(nu))!=length(theta)) stop("theta misspecified")
  D1sqrt <- Dmatrix(dd)
  D1 <- D1sqrt%*%D1sqrt
  #
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-D1sqrt%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  zeta<-matrix.sqrt(solve(D1))%*%lambda
  ypred <- numeric(length = nrow(dataPred))
  xpred<-matrix(nrow= nrow(dataPred),ncol=p)
  #
  for (indj in levels(dataPred$ind)) {
    #indj = levels(dataPred$ind)[1]
    dataFitj <- subset(dataFit,dataFit$ind==indj,select = c("ind",all.vars(formFixed),all.vars(formRandom)))
    dataPredj <- subset(dataPred,dataPred$ind==indj,select = c("ind",all.vars(formFixed),all.vars(formRandom)))
    njFit = nrow(dataFitj)
    njPred = nrow(dataPredj)
    seqFit = 1:njFit
    seqPred = njFit+1:njPred
    #
    dataPlus <- rbind(dataFitj,dataPredj)
    #
    xPlus1 <- model.matrix(formFixed,data=dataPlus)
    zPlus1<-model.matrix(formRandom,data=dataPlus)
    z1 <- matrix(zPlus1[seqFit,],ncol=ncol(zPlus1))
    x1 <- matrix(xPlus1[seqFit,],ncol=ncol(xPlus1))
    z1Pred <- matrix(zPlus1[seqPred,],ncol=ncol(zPlus1))
    x1Pred <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
    #
    medFit <- x1%*%beta1 + c.*z1%*%Deltab
    medPred <- x1Pred%*%beta1 + c.*z1Pred%*%Deltab
    #
    y1=dataFitj[,all.vars(formFixed)[1]]
    SigmaPlus = sigmae*CovCS(phiCS,njFit+njPred)
    PsiPlus<-(zPlus1)%*%(D1)%*%t(zPlus1)+SigmaPlus
    dj<-as.numeric(t(y1-medFit)%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit))
    LambdaPlus <- solve(solve(D1)+ t(zPlus1)%*%solve(SigmaPlus)%*%zPlus1)
    sPsiPlus <- solve(PsiPlus)
    lambdaBarPlus <- matrix.sqrt(sPsiPlus)%*%zPlus1%*%D1%*%zeta/as.numeric(sqrt(1+t(zeta)%*%LambdaPlus%*%zeta))
    vj <- matrix.sqrt(sPsiPlus)%*%lambdaBarPlus
    Psi22.1 <- PsiPlus[seqPred,seqPred]- PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]
    vjtil <- (vj[seqFit] + solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]%*%vj[seqPred])/
      as.numeric(sqrt(1+t(vj[seqPred])%*%Psi22.1%*%vj[seqPred]))
    Ajj<-as.numeric(t(vjtil)%*%(y1-medFit)) #as.numeric(mutj/sqrt(Mtj2))
    mu2.1 <- medPred + PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit)
    if (distr=="sn") tau1 <- dnorm(Ajj)/pnorm(Ajj)
    if (distr=="st") {
      dtj = gamma((nu+njFit)/2)/gamma(nu/2)/pi^(njFit/2)/sqrt(det(as.matrix(PsiPlus[seqFit,seqFit])))*nu^(-njFit/2)*
        (dj/nu+1)^(-(njFit+nu)/2)
      denST = 2*dtj*pt(sqrt(nu+njFit)*Ajj/sqrt(dj+nu),nu+njFit)
      tau1 <- as.numeric(dmvt(y1,delta=medFit,sigma=as.matrix(PsiPlus[seqFit,seqFit]),df=nu,log=F)*gamma((nu+njFit-1)/2)*(nu+dj)^((nu+njFit)/2)/
                           (denST*pi^.5*gamma((nu+njFit)/2)*(nu+dj+Ajj^2)^((nu+njFit-1)/2)))
    }
    if (distr=="ss") {
      f2 <- function(u) u^(nu - 1)*((2*pi)^(-njFit/2))*(u^(njFit/2))*((det(as.matrix(PsiPlus[seqFit,seqFit])))^(-1/2))*
        exp(-0.5*u*dj)*pnorm(u^(1/2)*Ajj)
      denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
      tau1<-2^(nu)*nu*gamma(nu-.5+njFit/2)*pgamma(1,nu-.5+njFit/2,(dj+Ajj^2)/2)/
        (denSS*(dj+Ajj^2)^(nu-.5+njFit/2)*pi^(njFit/2+.5)*det(as.matrix(PsiPlus[seqFit,seqFit]))^.5)
    }
    if (distr=="scn") {
      fy<-as.numeric(2*(nu[1]*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                          (1-nu[1])*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]))*pnorm(Ajj,0,1)))
      tau1<-as.numeric(2*(nu[1]*nu[2]^(-1/2)*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]/nu[2]))*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                            (1-nu[1])*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]))*dnorm(Ajj,0,1))/fy)
    }
    ypredj <- mu2.1 + tau1/as.numeric(sqrt(1+t(vj[seqPred])%*%Psi22.1%*%vj[seqPred]))*Psi22.1%*%vj[seqPred]
    ypred[dataPred$ind==indj] <- ypredj
    xpred[dataPred$ind==indj,] <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
  }
  xpred = as.data.frame(xpred)
  colnames(xpred) = colnames(xPlus1)
  if (all(xpred[,1]==1)) xpred=xpred[-1]
  data.frame(groupVar=dataPred$ind,xpred,ypred)
}


################################################################
#Log-likelihood - DEC
################################################################
ljnormalDEC <-function(j,y,x,z,time,beta1,Gammab,Deltab,sigmae,phiDEC,thetaDEC){
  c. = -sqrt(2/pi)
  y1=y[j]
  t1= time[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,thetaDEC,t1)#CovARp(phiAR,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))
}
#
ljtDEC <-function(j,nu,y,x,z,time,beta1,Gammab,Deltab,sigmae,phiDEC,thetaDEC){
  c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  y1=y[j]
  t1= time[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,thetaDEC,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med))
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  #log(2*dmvt(y1,delta = med, sigma = Psi, df = nu,log=F)*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))#veroST1(Psi,Ajj,dj,nu,pp=njj))
  log(2*dtj*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))
}
#
ljsDEC <-function(j,nu,y,x,z,time,beta1,Gammab,Deltab,sigmae,phiDEC,thetaDEC){
  c.=-sqrt(2/pi)*nu/(nu-.5)
  y1=y[j]
  t1= time[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,thetaDEC,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med))
  #f <- function(u) u^(nu - 1)*dmvnorm(y1,med,Psi/u)*pnorm(u^(1/2)*Ajj)
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
  resp <- integrate(Vectorize(f2),0,1)$value
  log(2*nu*resp)
}
#
ljcnDEC <-function(j,nu,y,x,z,time,beta1,Gammab,Deltab,sigmae,phiDEC,thetaDEC){
  c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[j]
  t1= time[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,thetaDEC,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(sqrt(nu[2])*Ajj,0,1)+
           (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
}

logveroDEC = function(y,x,z,time,ind,beta1,sigmae,phiDEC,thetaDEC,D1,lambda,distr,nu){ #ind = indicadora de individuo

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  N <-length(ind)

  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalDEC,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtDEC,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsDEC,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnDEC,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  lv
}

##############################################################################
# EM - DEC
##############################################################################
calcbi_emjDEC <- function(jseq,y,x,z,time,beta1,Gammab, Deltab,sigmae,phiDEC,thetaDEC,zeta,
                          distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  #
  y1=y[jseq]
  t1=time[jseq]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  nj = length(y1)
  Sigma = sigmae*CovDEC(phiDEC,thetaDEC,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ajj<-as.numeric(mutj/sqrt(Mtj2))
  D1<- Gammab+Deltab%*%t(Deltab)
  #
  mediab<-D1%*%t(z1)%*%solve(Psi)%*%(y1-med)+c.*Deltab
  Lambda<-solve(solve(D1)+t(z1)%*%solve(Sigma)%*%z1)
  #
  if  (distr=="sn"){
    bi<-mediab+Lambda%*%zeta/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))*as.numeric(dnorm(Ajj,0,1))/as.numeric(pnorm(Ajj,0,1))
  }

  if (distr=="st"){
    dtj = gamma((nu+nj)/2)/gamma(nu/2)/pi^(nj/2)/sqrt(det(Psi))*nu^(-nj/2)*(dj/nu+1)^(-(nj+nu)/2)
    denST = 2*dtj*pt(sqrt(nu+nj)*Ajj/sqrt(dj+nu),nu+nj)
    esper2<- as.numeric(dmvt(y1,delta=med,sigma=Psi,df=nu,log=F)*gamma((nu+nj-1)/2)*(nu+dj)^((nu+nj)/2)/
                          (denST*pi^.5*gamma((nu+nj)/2)*(nu+dj+Ajj^2)^((nu+nj-1)/2)))
    bi<-mediab+Lambda%*%zeta/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))*esper2
  }

  if (distr=="ss"){
    f2 <- function(u) u^(nu - 1)*((2*pi)^(-nj/2))*(u^(nj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
    denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
    esper2<-2^(nu)*nu*gamma(nu-.5+nj/2)*pgamma(1,nu-.5+nj/2,(dj+Ajj^2)/2)/
      (denSS*(dj+Ajj^2)^(nu-.5+nj/2)*pi^(nj/2+.5)*det(Psi)^.5)
    bi<-mediab+Lambda%*%zeta*esper2/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))
  }

  if (distr=="scn"){
    fy<-as.numeric(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                        (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
    esper2<-as.numeric(2*(nu[1]*nu[2]^(-1/2)*dmvnorm(y1,med,Psi/nu[2])*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                            (1-nu[1])*dmvnorm(y1,med,Psi)*dnorm(Ajj,0,1))/fy)
    bi<-mediab+Lambda%*%zeta*esper2/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))
  }
  bi
}
emjDEC = function(jseq, y, x, z,time, beta1, Gammab, Deltab, sigmae,phiDEC,thetaDEC,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  #
  y1=y[jseq]
  t1=time[jseq]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  nj = length(y1)
  Sigma = sigmae*CovDEC(phiDEC,thetaDEC,t1)#CovARp(phi = estphit(piAR),t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ajj<-as.numeric(mutj/sqrt(Mtj2))
  D1<- Gammab+Deltab%*%t(Deltab)
  #
  if  (distr=="sn"){
    uj<-1
    esper<-as.numeric(dnorm(Ajj,0,1)/pnorm(Ajj,0,1))
  }

  if (distr=="st"){
    dtj = gamma((nu+nj)/2)/gamma(nu/2)/pi^(nj/2)/sqrt(det(Psi))*nu^(-nj/2)*(dj/nu+1)^(-(nj+nu)/2)
    denST = 2*dtj*pt(sqrt(nu+nj)*Ajj/sqrt(dj+nu),nu+nj)
    uj<-as.numeric(2^2*(nu+dj)^(-(nu+nj+2)/2)*gamma((nj+nu+2)/2)*nu^(nu/2)*
                     pt(sqrt((nj+nu+2)/(dj+nu))*Ajj,nj+nu+2)/(gamma(nu/2)*det(Psi)^.5*denST*pi^(nj/2)))
    esper <-as.numeric(2*nu^(nu/2)*gamma((nj+nu+1)/2)/(denST*pi^((nj+1)/2)*gamma(nu/2)*det(Psi)^.5*(nu+dj+Ajj^2)^((nu+nj+1)/2)))
  }

  if (distr=="ss"){
    f2esp <- function(s) pnorm(s^.5*Ajj)*dgamma(s,nu+1+nj/2,dj/2)/pgamma(1,nu+1+nj/2,dj/2)
    EspVal <- integrate(f2esp,0,1)$value#mean(pnorm(S^(1/2)*Ajj))#
    f2 <- function(u) u^(nu - 1)*((2*pi)^(-nj/2))*(u^(nj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
    denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
    uj<-2^(2+nu)*nu*gamma(nu+1+nj/2)*pgamma(1,nu+1+nj/2,dj/2)*EspVal/
      (denSS*dj^(nu+1+nj/2)*pi^(nj/2)*det(Psi)^.5)
    esper <- 2^(1+nu)*nu*gamma(nu+.5+nj/2)*pgamma(1,nu+.5+nj/2,(dj+Ajj^2)/2)/
      (denSS*(dj+Ajj^2)^(nu+.5+nj/2)*pi^(nj/2+.5)*det(Psi)^.5)
  }

  if (distr=="scn"){
    fy<-as.numeric(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                        (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
    uj<-2*(nu[1]*nu[2]*dmvnorm(y1,med,Psi/nu[2])*pnorm(nu[2]^(1/2)*Ajj,0,1)+
             (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))/fy
    esper<-as.numeric(2*(nu[1]*nu[2]^(1/2)*dmvnorm(y1,med,Psi/nu[2])*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                           (1-nu[1])*dmvnorm(y1,med,Psi)*dnorm(Ajj,0,1))/fy)
  }

  sSigma = solve(Sigma)
  sRi = sSigma*sigmae
  Tbj<-solve(solve(Gammab)+t(z1)%*%sSigma%*%z1)
  r<-Tbj%*%t(z1)%*%sSigma%*%(y1-x1%*%beta1)
  s1<-(diag(q1)-Tbj%*%t(z1)%*%sSigma%*%z1)%*%Deltab
  utj<-as.numeric(uj*(mutj+c.)+sqrt(Mtj2)*esper)
  ut2j<-as.numeric(uj*(mutj+c.)^2+Mtj2+sqrt(Mtj2)*(mutj+2*c.)*esper)
  ub<-uj*r+s1*utj
  utbj<- r*utj+s1*ut2j
  ub2j<-Tbj+uj*r%*%t(r)+s1%*%t(r)*utj+r%*%t(s1)*utj+s1%*%t(s1)*ut2j
  #
  sum1<-uj*t(x1)%*%sRi%*%x1 #denom beta
  sum2<-(t(x1)%*%sRi%*%(uj*y1-z1%*%ub)) #num beta
  sum3<-uj*t(y1-x1%*%beta1)%*%sRi%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%sRi%*%z1%*%ub-
    t(ub)%*%t(z1)%*%sRi%*%(y1-x1%*%beta1)+traceM(sRi%*%z1%*%ub2j%*%t(z1)) #soma do sig2
  sum4<-ub2j-utbj%*%t(Deltab)-Deltab%*%t(utbj)+
    ut2j*Deltab%*%t(Deltab) #soma do Gamma
  sum5<-utbj #num do delta
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,sum5=sum5,ut2j=ut2j,
                 uj=uj,ubj=ub,ub2j=ub2j)
  #if (calcbi) obj.out$bi=bi
  return(obj.out)
}

#função para maximizar
lcDEC <- function(parDEC,beta1,sigmae,y,x,z,time,ind,u,ub,ub2) {
  #print(parDEC)
  phiDEC <- parDEC[1]
  thetaDEC<-parDEC[2]#2.128601
  m<-n_distinct(ind)
  N<-length(ind)
  indlevels <- levels(ind)
  soma <-0
  for (i in seq_len(m)) {
    jseq <- which(ind==indlevels[i])
    y1=y[jseq]
    t1=time[jseq]
    p= ncol(x)
    q1=ncol(z)
    x1=matrix(x[jseq,  ],ncol=p)
    z1=matrix(z[jseq,  ],ncol=q1)
    med = x1%*%beta1
    nj = length(y1)
    covmat = CovDEC(phiDEC,thetaDEC,t1)
    Sigma = covmat*sigmae
    sSigma = solve(Sigma)
    indi = which(names(u)==indlevels[i])
    uj = u[[indi]]
    ubj = ub[[indi]]
    ub2j = ub2[[indi]]
    detSigma = det(Sigma)
    if (sum(eigen(solve(Sigma))$values<=0)>0) {soma = -thetaDEC^100;break}
    else {
      soma = soma + as.numeric(-.5*log(detSigma)-.5*uj*t(y1-med)%*%sSigma%*%(y1-med)+
                                 t(y1-med)%*%sSigma%*%z1%*%ubj-.5*traceM(sSigma%*%z1%*%ub2j%*%t(z1)))
    }
  }
  soma
}

EM.SkewDEC<- function(formFixed,formRandom,data,groupVar,timeVar,
                      beta1,sigmae,D1,lambda,distr,nu,parDEC,lb,lu,luDEC,
                      precisao,informa,max.iter,showiter,showerroriter){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]
  #
  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x)
  q1<-ncol(z)
  #
  if (!is.null(parDEC)) {
    if (length(parDEC)!=2) stop ("initial value from phi should have length 2 or NULL")
    if (parDEC[1]<=0|parDEC[1]>=1) stop("invalid initial value from phi1")
    if (parDEC[2]<=0) stop("invalid initial value from phi2")
    if (parDEC[2]>= luDEC) stop("initial value from phi2 must be smaller than luDEC")
  }

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  #zeta<-matrix.sqrt(solve(D1))%*%lambda

  if (is.null(parDEC)) {
    #cat("calculating initial values for DEC... \n")
    thetat<- seq(0.1,2,by=.1)
    phit <- seq(0.1,.9,by=.05)
    vect <-merge(phit,thetat,all=T)
    logveroDECv<-function(phitheta){logveroDEC(y, x, z,time,ind, beta1=beta1, sigmae=sigmae,
                                               phiDEC=phitheta[1],thetaDEC=phitheta[2],
                                               D1=D1,lambda=lambda, distr=distr, nu=nu)}
    logverovec <- apply(vect,1,logveroDECv)
    parDEC <- as.numeric(vect[which.max(logverovec),])
  }
  phiDEC=parDEC[1]
  thetaDEC=parDEC[2]

  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,phiDEC,thetaDEC,nu)

  criterio<-10
  count<-0
  llji = logveroDEC(y, x, z, time,ind, beta1, sigmae,phiDEC,thetaDEC, D1, lambda, distr, nu)
  if (is.nan(llji)|is.infinite(abs(llji))) stop("NaN/infinity initial likelihood")

  while((criterio > precisao)&(count<max.iter)){

    count <- count + 1
    res_emj = revert_list(tapply(1:N,ind,emjDEC,y=y, x=x, z=z,time=time, beta1=beta1, Gammab=Gammab,
                                 Deltab=Deltab, sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC, distr=distr,nu=nu))
    sum1 = Reduce("+",res_emj$sum1)
    sum2 = Reduce("+",res_emj$sum2)
    sum3 = sum(unlist(res_emj$sum3))
    sum4 = Reduce("+",res_emj$sum4)
    sum5 = Reduce("+",res_emj$sum5)
    ut2j = unlist(res_emj$ut2j,use.names = F)

    #if (calcbi) bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
    beta1<-solve(sum1)%*%sum2
    sigmae<-as.numeric(sum3)/N
    Gammab<-sum4/m
    Deltab<-sum5/sum(ut2j)
    #
    D1<-Gammab+Deltab%*%t(Deltab)
    lambda<-matrix.sqrt(solve(D1))%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%solve(D1)%*%Deltab))
    #zeta<-matrix.sqrt(solve(D1))%*%lambda
    #
    parDEC<- optim(c(phiDEC,thetaDEC),lcDEC,gr = NULL,method = "L-BFGS-B", lower =rep(0.0001,2),
                   upper = c(.9999,luDEC),control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
    phiDEC<-parDEC[1]
    thetaDEC<-parDEC[2]
    #
    logvero1<-function(nu){logveroDEC(y, x, z,time, ind, beta1, sigmae,phiDEC,thetaDEC, D1, lambda, distr, nu)}

    if (distr=="sn"){ nu<-NULL} else
    {nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par}
    #
    param <- teta
    teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,phiDEC,thetaDEC,nu)
    criterio2 <- as.numeric(sqrt((teta-param)%*%(teta-param)))
    llj1<-llji
    llji <- logveroDEC(y, x, z, time,ind, beta1, sigmae,phiDEC,thetaDEC, D1, lambda, distr, nu)
    criterio <- abs((llji-llj1)/llj1)
    if (showiter&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio,"\r") #  criterium ",criterio," or ",criterio2,"\r")
  }
  if (count==max.iter) message("\n maximum number of iterations reachead")
  cat("\n")
  zeta<-matrix.sqrt(solve(D1))%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjDEC,y=y, x=x, z=z, time=time, beta1=beta1, Gammab=Gammab,
                             Deltab=Deltab, sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC,zeta=zeta, distr=distr,nu=nu,simplify = FALSE)),ncol=q1,byrow = T)
  # bi = t(list.cbind(tapply(1:N,ind,calcbi_emjDEC,y=y, x=x, z=z, time=time, beta1=beta1, Gammab=Gammab,
  #                         Deltab=Deltab, sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC,zeta=zeta, distr=distr,nu=nu,simplify = FALSE)))
  dd<-matrix.sqrt(D1)[upper.tri(D1, diag = T)]
  theta = c(beta1,sigmae,phiDEC,thetaDEC,dd,lambda,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phi1DEC","phi2DEC",paste0("Dsqrt",1:length(dd)),paste0("lambda",1:q1))
  else names(theta)<- c(colnames(x),"sigma2","phi1DEC","phi2DEC",paste0("Dsqrt",1:length(dd)),paste0("lambda",1:q1),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                  phi=c(phiDEC,thetaDEC),dsqrt=dd,D=D1,lambda=as.numeric(lambda)),
                  uhat=unlist(res_emj$uj))
  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixDEC(y,x,z,time,ind,beta1,sigmae,phiDEC,thetaDEC,D1,lambda,distr = distr,nu = nu),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      q2<-q1*(q1+1)/2
      desvios[(p+q2+4):(p+3+q2+q1)] <- rep(NA,q1)
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }

  obj.out$loglik <-as.numeric(llji)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=criterio
  obj.out
}

predictf.skewDEC<- function(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr,theta){
  dataPred[,all.vars(formFixed)[1]] <- 0
  dataFit$ind <-dataFit[,groupVar]
  dataPred$ind <-dataPred[,groupVar]
  dataPred$ind <- droplevels(dataPred$ind)
  #
  #theta = beta1,sigmae,phiAR,D1,lambda,nu
  p <- ncol(model.matrix(formFixed,data=dataPred))
  q1 <- ncol(model.matrix(formRandom,data=dataPred))
  q2 <- q1*(q1+1)/2
  beta1 <- matrix(theta[1:p],ncol=1)
  sigmae <- as.numeric(theta[p+1])
  phiDEC <- as.numeric(theta[(p+2)])
  thetaDEC <- as.numeric(theta[(p+3)])
  dd <- theta[(p+4):(p+3+q2)]
  lambda <- matrix(theta[(p+q2+4):(p+3+q2+q1)],ncol=1)
  if (distr=="sn") {nu<- NULL
                    c. = -sqrt(2/pi)}
  if (distr=="st") {nu<- theta[p+q2+q1+4]
                    c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {nu<- theta[p+q2+q1+4]
                    c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {nu<- theta[(p+q2+q1+4):(p+q2+q1+5)]
                    c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  if ((p+3+q2+q1+length(nu))!=length(theta)) stop("theta misspecified")
  D1sqrt <- Dmatrix(dd)
  D1 <- D1sqrt%*%D1sqrt
  #
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-D1sqrt%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  zeta<-matrix.sqrt(solve(D1))%*%lambda
  ypred <- numeric(length = nrow(dataPred))
  timepred <- numeric(length = nrow(dataPred))
  xpred<-matrix(nrow= nrow(dataPred),ncol=p)
  #
  for (indj in levels(dataPred$ind)) {
    #indj = levels(dataPred$ind)[1]
    dataFitj <- subset(dataFit,dataFit$ind==indj,select = c("ind",all.vars(formFixed),all.vars(formRandom),timeVar))
    dataPredj <- subset(dataPred,dataPred$ind==indj,select = c("ind",all.vars(formFixed),all.vars(formRandom),timeVar))
    if (!is.null(timeVar)) {
      dataFitj$time <- dataFitj[,timeVar]
      dataPredj$time <- dataPredj[,timeVar]
    }
    njFit = nrow(dataFitj)
    njPred = nrow(dataPredj)
    seqFit = 1:njFit
    seqPred = njFit+1:njPred
    #
    if (is.null(timeVar)) {
      dataFitj$time<- seqFit
      dataPredj$time<- seqPred
    }
    dataPlus <- rbind(dataFitj,dataPredj)
    #
    xPlus1 <- model.matrix(formFixed,data=dataPlus)
    zPlus1<-model.matrix(formRandom,data=dataPlus)
    z1 <- matrix(zPlus1[seqFit,],ncol=ncol(zPlus1))
    x1 <- matrix(xPlus1[seqFit,],ncol=ncol(xPlus1))
    z1Pred <- matrix(zPlus1[seqPred,],ncol=ncol(zPlus1))
    x1Pred <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
    #
    medFit <- x1%*%beta1 + c.*z1%*%Deltab
    medPred <- x1Pred%*%beta1 + c.*z1Pred%*%Deltab
    #
    y1=dataFitj[,all.vars(formFixed)[1]]
    SigmaPlus = sigmae*CovDEC(phiDEC,thetaDEC,c(dataFitj$time,dataPredj$time))
    PsiPlus<-(zPlus1)%*%(D1)%*%t(zPlus1)+SigmaPlus
    dj<-as.numeric(t(y1-medFit)%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit))
    LambdaPlus <- solve(solve(D1)+ t(zPlus1)%*%solve(SigmaPlus)%*%zPlus1)
    sPsiPlus <- solve(PsiPlus)
    lambdaBarPlus <- matrix.sqrt(sPsiPlus)%*%zPlus1%*%D1%*%zeta/as.numeric(sqrt(1+t(zeta)%*%LambdaPlus%*%zeta))
    vj <- matrix.sqrt(sPsiPlus)%*%lambdaBarPlus
    Psi22.1 <- PsiPlus[seqPred,seqPred]- PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]
    vjtil <- (vj[seqFit] + solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]%*%vj[seqPred])/
      as.numeric(sqrt(1+t(vj[seqPred])%*%Psi22.1%*%vj[seqPred]))
    Ajj<-as.numeric(t(vjtil)%*%(y1-medFit)) #as.numeric(mutj/sqrt(Mtj2))
    mu2.1 <- medPred + PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit)
    if (distr=="sn") tau1 <- dnorm(Ajj)/pnorm(Ajj)
    if (distr=="st") {
      dtj = gamma((nu+njFit)/2)/gamma(nu/2)/pi^(njFit/2)/sqrt(det(as.matrix(PsiPlus[seqFit,seqFit])))*nu^(-njFit/2)*
        (dj/nu+1)^(-(njFit+nu)/2)
      denST = 2*dtj*pt(sqrt(nu+njFit)*Ajj/sqrt(dj+nu),nu+njFit)
      tau1 <- as.numeric(dmvt(y1,delta=medFit,sigma=as.matrix(PsiPlus[seqFit,seqFit]),df=nu,log=F)*gamma((nu+njFit-1)/2)*(nu+dj)^((nu+njFit)/2)/
                           (denST*pi^.5*gamma((nu+njFit)/2)*(nu+dj+Ajj^2)^((nu+njFit-1)/2)))
    }
    if (distr=="ss") {
      f2 <- function(u) u^(nu - 1)*((2*pi)^(-njFit/2))*(u^(njFit/2))*((det(as.matrix(PsiPlus[seqFit,seqFit])))^(-1/2))*
        exp(-0.5*u*dj)*pnorm(u^(1/2)*Ajj)
      denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
      tau1<-2^(nu)*nu*gamma(nu-.5+njFit/2)*pgamma(1,nu-.5+njFit/2,(dj+Ajj^2)/2)/
        (denSS*(dj+Ajj^2)^(nu-.5+njFit/2)*pi^(njFit/2+.5)*det(as.matrix(PsiPlus[seqFit,seqFit]))^.5)
    }
    if (distr=="scn") {
      fy<-as.numeric(2*(nu[1]*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                          (1-nu[1])*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]))*pnorm(Ajj,0,1)))
      tau1<-as.numeric(2*(nu[1]*nu[2]^(-1/2)*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]/nu[2]))*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                            (1-nu[1])*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]))*dnorm(Ajj,0,1))/fy)
    }
    ypredj <- mu2.1 + tau1/as.numeric(sqrt(1+t(vj[seqPred])%*%Psi22.1%*%vj[seqPred]))*Psi22.1%*%vj[seqPred]
    ypred[dataPred$ind==indj] <- ypredj
    xpred[dataPred$ind==indj,] <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
    timepred[dataPred$ind==indj] <- dataPredj$time
  }
  xpred = as.data.frame(xpred)
  colnames(xpred) = colnames(xPlus1)
  if (all(xpred[,1]==1)) xpred=xpred[-1]
  if (is.null(timeVar)) dfout = data.frame(groupVar=dataPred$ind,xpred,ypred)
  else dfout = data.frame(groupVar=dataPred$ind,time=timepred,xpred,ypred)
  dfout
}


################################################################
#Log-likelihood - CAR(1)
################################################################
ljnormalCAR1 <-function(j,y,x,z,time,beta1,Gammab,Deltab,sigmae,phiDEC){
  c. = -sqrt(2/pi)
  y1=y[j]
  t1= time[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,1,t1)#CovARp(phiAR,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))
}
#
ljtCAR1 <-function(j,nu,y,x,z,time,beta1,Gammab,Deltab,sigmae,phiDEC){
  c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  y1=y[j]
  t1= time[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,1,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med))
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  #log(2*dmvt(y1,delta = med, sigma = Psi, df = nu,log=F)*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))#veroST1(Psi,Ajj,dj,nu,pp=njj))
  log(2*dtj*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))
}
#
ljsCAR1 <-function(j,nu,y,x,z,time,beta1,Gammab,Deltab,sigmae,phiDEC){
  c.=-sqrt(2/pi)*nu/(nu-.5)
  y1=y[j]
  t1= time[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,1,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med))
  #f <- function(u) u^(nu - 1)*dmvnorm(y1,med,Psi/u)*pnorm(u^(1/2)*Ajj)
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
  resp <- integrate(Vectorize(f2),0,1)$value
  log(2*nu*resp)
}
#
ljcnCAR1 <-function(j,nu,y,x,z,time,beta1,Gammab,Deltab,sigmae,phiDEC){
  c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[j]
  t1= time[j]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,1,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(sqrt(nu[2])*Ajj,0,1)+
           (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
}

logveroCAR1 = function(y,x,z,time,ind,beta1,sigmae,phiDEC,D1,lambda,distr,nu){ #ind = indicadora de individuo

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  N <-length(ind)

  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalCAR1,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtCAR1,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsCAR1,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnCAR1,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC))
  lv
}

##############################################################################
# EM - CAR(1)
##############################################################################
emjCAR1 = function(jseq, y, x, z,time, beta1, Gammab, Deltab, sigmae,phiDEC,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  #
  y1=y[jseq]
  t1=time[jseq]
  p= ncol(x)
  q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1 + c.*z1%*%Deltab
  nj = length(y1)
  Sigma = sigmae*CovDEC(phiDEC,1,t1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ajj<-as.numeric(mutj/sqrt(Mtj2))
  D1<- Gammab+Deltab%*%t(Deltab)
  #
  if  (distr=="sn"){
    uj<-1
    esper<-as.numeric(dnorm(Ajj,0,1)/pnorm(Ajj,0,1))
  }

  if (distr=="st"){
    dtj = gamma((nu+nj)/2)/gamma(nu/2)/pi^(nj/2)/sqrt(det(Psi))*nu^(-nj/2)*(dj/nu+1)^(-(nj+nu)/2)
    denST = 2*dtj*pt(sqrt(nu+nj)*Ajj/sqrt(dj+nu),nu+nj)
    uj<-as.numeric(2^2*(nu+dj)^(-(nu+nj+2)/2)*gamma((nj+nu+2)/2)*nu^(nu/2)*
                     pt(sqrt((nj+nu+2)/(dj+nu))*Ajj,nj+nu+2)/(gamma(nu/2)*det(Psi)^.5*denST*pi^(nj/2)))
    esper <-as.numeric(2*nu^(nu/2)*gamma((nj+nu+1)/2)/(denST*pi^((nj+1)/2)*gamma(nu/2)*det(Psi)^.5*(nu+dj+Ajj^2)^((nu+nj+1)/2)))
  }

  if (distr=="ss"){
    f2esp <- function(s) pnorm(s^.5*Ajj)*dgamma(s,nu+1+nj/2,dj/2)/pgamma(1,nu+1+nj/2,dj/2)
    EspVal <- integrate(f2esp,0,1)$value#mean(pnorm(S^(1/2)*Ajj))#
    f2 <- function(u) u^(nu - 1)*((2*pi)^(-nj/2))*(u^(nj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
    denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
    uj<-2^(2+nu)*nu*gamma(nu+1+nj/2)*pgamma(1,nu+1+nj/2,dj/2)*EspVal/
      (denSS*dj^(nu+1+nj/2)*pi^(nj/2)*det(Psi)^.5)
    esper <- 2^(1+nu)*nu*gamma(nu+.5+nj/2)*pgamma(1,nu+.5+nj/2,(dj+Ajj^2)/2)/
      (denSS*(dj+Ajj^2)^(nu+.5+nj/2)*pi^(nj/2+.5)*det(Psi)^.5)
  }

  if (distr=="scn"){
    fy<-as.numeric(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                        (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
    uj<-2*(nu[1]*nu[2]*dmvnorm(y1,med,Psi/nu[2])*pnorm(nu[2]^(1/2)*Ajj,0,1)+
             (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))/fy
    esper<-as.numeric(2*(nu[1]*nu[2]^(1/2)*dmvnorm(y1,med,Psi/nu[2])*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                           (1-nu[1])*dmvnorm(y1,med,Psi)*dnorm(Ajj,0,1))/fy)
  }

  sSigma = solve(Sigma)
  sRi = sSigma*sigmae
  Tbj<-solve(solve(Gammab)+t(z1)%*%sSigma%*%z1)
  r<-Tbj%*%t(z1)%*%sSigma%*%(y1-x1%*%beta1)
  s1<-(diag(q1)-Tbj%*%t(z1)%*%sSigma%*%z1)%*%Deltab
  utj<-as.numeric(uj*(mutj+c.)+sqrt(Mtj2)*esper)
  ut2j<-as.numeric(uj*(mutj+c.)^2+Mtj2+sqrt(Mtj2)*(mutj+2*c.)*esper)
  ub<-uj*r+s1*utj
  utbj<- r*utj+s1*ut2j
  ub2j<-Tbj+uj*r%*%t(r)+s1%*%t(r)*utj+r%*%t(s1)*utj+s1%*%t(s1)*ut2j
  #
  sum1<-uj*t(x1)%*%sRi%*%x1 #denom beta
  sum2<-(t(x1)%*%sRi%*%(uj*y1-z1%*%ub)) #num beta
  sum3<-uj*t(y1-x1%*%beta1)%*%sRi%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%sRi%*%z1%*%ub-
    t(ub)%*%t(z1)%*%sRi%*%(y1-x1%*%beta1)+traceM(sRi%*%z1%*%ub2j%*%t(z1)) #soma do sig2
  sum4<-ub2j-utbj%*%t(Deltab)-Deltab%*%t(utbj)+
    ut2j*Deltab%*%t(Deltab) #soma do Gamma
  sum5<-utbj #num do delta
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,sum5=sum5,ut2j=ut2j,
                 uj=uj,ubj=ub,ub2j=ub2j)
  #if (calcbi) obj.out$bi=bi
  return(obj.out)
}

#função para maximizar
lcCAR1 <- function(phiCAR,beta1,sigmae,y,x,z,time,ind,u,ub,ub2) {
  m<-n_distinct(ind)
  N<-length(ind)
  indlevels <- levels(ind)
  soma <-0
  for (i in seq_len(m)) {
    jseq <- which(ind==indlevels[i])
    y1=y[jseq]
    t1=time[jseq]
    p= ncol(x)
    q1=ncol(z)
    x1=matrix(x[jseq,  ],ncol=p)
    z1=matrix(z[jseq,  ],ncol=q1)
    med = x1%*%beta1
    nj = length(y1)
    Sigma = CovDEC(phiCAR,1,t1)*sigmae
    sSigma = solve(Sigma)
    indi = which(names(u)==indlevels[i])
    uj = u[[indi]]
    ubj = ub[[indi]]
    ub2j = ub2[[indi]]
    soma = soma + as.numeric(-.5*log(det(Sigma))-.5*uj*t(y1-med)%*%sSigma%*%(y1-med)+
                               t(y1-med)%*%sSigma%*%z1%*%ubj-.5*traceM(sSigma%*%z1%*%ub2j%*%t(z1)))
  }
  soma
}

EM.SkewCAR1<- function(formFixed,formRandom,data,groupVar,timeVar,
                       distr,beta1,sigmae,phiCAR1,D1,lambda,nu,lb,lu,
                       precisao,informa,max.iter,showiter,showerroriter){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]
  #
  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x)
  q1<-ncol(z)
  #
  if (!is.null(phiCAR1) & length(phiCAR1)!=1) stop("initial value from phi must have length 1 or be NULL")
  if (!is.null(phiCAR1)) if (phiCAR1>=1 | phiCAR1<=0) stop ("0<initialValue$phi<1 needed")

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  #zeta<-matrix.sqrt(solve(D1))%*%lambda

  if (is.null(phiCAR1)) {
    lmeCAR = try(lme(formFixed,random=~1|ind,data=data,correlation=corCAR1(form = ~time)),silent=T)
    if (class(lmeCAR)=="try-error") phiDEC =abs(as.numeric(pacf(y-x%*%beta1,lag.max=1,plot=F)$acf))
    else {
      phiDEC = capture.output(lmeCAR$modelStruct$corStruct)[3]
      phiDEC = as.numeric(strsplit(phiDEC, " ")[[1]])
    }
  } else phiDEC <- phiCAR1

  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,phiDEC,nu)

  criterio<-10
  count<-0
  llji = logveroCAR1(y, x, z, time,ind, beta1, sigmae,phiDEC, D1, lambda, distr, nu)
  if (is.nan(llji)|is.infinite(abs(llji))) stop("NaN/infinity initial likelihood")

  while((criterio > precisao)&(count<max.iter)){
    #print(nu)

    count <- count + 1
    res_emj = revert_list(tapply(1:N,ind,emjCAR1,y=y, x=x, z=z,time=time, beta1=beta1, Gammab=Gammab,
                                 Deltab=Deltab, sigmae=sigmae,phiDEC=phiDEC,distr=distr,nu=nu))
    sum1 = Reduce("+",res_emj$sum1)
    sum2 = Reduce("+",res_emj$sum2)
    sum3 = sum(unlist(res_emj$sum3))
    sum4 = Reduce("+",res_emj$sum4)
    sum5 = Reduce("+",res_emj$sum5)
    ut2j = unlist(res_emj$ut2j,use.names = F)

    #if (calcbi) bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
    beta1<-solve(sum1)%*%sum2
    sigmae<-as.numeric(sum3)/N
    Gammab<-sum4/m
    Deltab<-sum5/sum(ut2j)
    #
    D1<-Gammab+Deltab%*%t(Deltab)
    lambda<-matrix.sqrt(solve(D1))%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%solve(D1)%*%Deltab))
    #zeta<-matrix.sqrt(solve(D1))%*%lambda
    #
    phiDEC<- optim(phiDEC,lcCAR1,gr = NULL,method = "L-BFGS-B", lower =0.0001,
                   upper = .9999,control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
    #
    logvero1<-function(nu){logveroCAR1(y, x, z,time, ind, beta1, sigmae,phiDEC, D1, lambda, distr, nu)}

    if (distr=="sn"){ nu<-NULL} else
    {nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par}
    #
    param <- teta
    teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,phiDEC,nu)
    criterio2 <- as.numeric(sqrt((teta-param)%*%(teta-param)))
    llj1<-llji
    llji <- logveroCAR1(y, x, z, time,ind, beta1, sigmae,phiDEC, D1, lambda, distr, nu)
    criterio <- abs((llji-llj1)/llj1)
    if (showiter&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio,"\r") #  criterium ",criterio," or ",criterio2,"\r")
  }
  if (count==max.iter) message("\n maximum number of iterations reachead")
  cat("\n")
  zeta<-matrix.sqrt(solve(D1))%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjDEC,y=y, x=x, z=z, time=time, beta1=beta1, Gammab=Gammab,
                             Deltab=Deltab, sigmae=sigmae,phiDEC=phiDEC,thetaDEC=1,zeta=zeta, distr=distr,nu=nu,simplify = FALSE)),ncol=q1,byrow = T)
  # bi = t(list.cbind(tapply(1:N,ind,calcbi_emjDEC,y=y, x=x, z=z, time=time, beta1=beta1, Gammab=Gammab,
  #                         Deltab=Deltab, sigmae=sigmae,phiDEC=phiDEC,thetaDEC=1,zeta=zeta, distr=distr,nu=nu,simplify = FALSE)))
  dd<-matrix.sqrt(D1)[upper.tri(D1, diag = T)]
  theta = c(beta1,sigmae,phiDEC,dd,lambda,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phiCAR1",paste0("Dsqrt",1:length(dd)),paste0("lambda",1:q1))
  else names(theta)<- c(colnames(x),"sigma2","phiCAR1",paste0("Dsqrt",1:length(dd)),paste0("lambda",1:q1),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                            phi=phiDEC,dsqrt=dd,D=D1,lambda=as.numeric(lambda)),
                  uhat=unlist(res_emj$uj))
  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixCAR1(y,x,z,time,ind,beta1,sigmae,phiDEC,D1,lambda,distr = distr,nu = nu),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      q2<-q1*(q1+1)/2
      desvios[(p+q2+3):(p+2+q2+q1)] <- rep(NA,q1)
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(llji)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=criterio
  obj.out
}

predictf.skewCAR1<- function(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr,theta){
  dataPred[,all.vars(formFixed)[1]] <- 0
  dataFit$ind <-dataFit[,groupVar]
  dataPred$ind <-dataPred[,groupVar]
  dataPred$ind <- droplevels(dataPred$ind)
  #
  #theta = beta1,sigmae,phiAR,D1,lambda,nu
  p <- ncol(model.matrix(formFixed,data=dataPred))
  q1 <- ncol(model.matrix(formRandom,data=dataPred))
  q2 <- q1*(q1+1)/2
  beta1 <- matrix(theta[1:p],ncol=1)
  sigmae <- as.numeric(theta[p+1])
  phiDEC <- as.numeric(theta[(p+2)])
  dd <- theta[(p+3):(p+2+q2)]
  lambda <- matrix(theta[(p+q2+3):(p+2+q2+q1)],ncol=1)
  if (distr=="sn") {nu<- NULL
                    c. = -sqrt(2/pi)}
  if (distr=="st") {nu<- theta[p+q2+q1+3]
                    c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {nu<- theta[p+q2+q1+3]
                    c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {nu<- theta[(p+q2+q1+3):(p+q2+q1+4)]
                    c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  if ((p+2+q2+q1+length(nu))!=length(theta)) stop("theta misspecified")
  D1sqrt <- Dmatrix(dd)
  D1 <- D1sqrt%*%D1sqrt
  #
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-D1sqrt%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  zeta<-matrix.sqrt(solve(D1))%*%lambda
  ypred <- numeric(length = nrow(dataPred))
  timepred <- numeric(length = nrow(dataPred))
  xpred<-matrix(nrow= nrow(dataPred),ncol=p)
  #
  for (indj in levels(dataPred$ind)) {
    #indj = levels(dataPred$ind)[1]
    dataFitj <- subset(dataFit,dataFit$ind==indj,select = c("ind",all.vars(formFixed),all.vars(formRandom),timeVar))
    dataPredj <- subset(dataPred,dataPred$ind==indj,select = c("ind",all.vars(formFixed),all.vars(formRandom),timeVar))
    if (!is.null(timeVar)) {
      dataFitj$time <- dataFitj[,timeVar]
      dataPredj$time <- dataPredj[,timeVar]
    }
    njFit = nrow(dataFitj)
    njPred = nrow(dataPredj)
    seqFit = 1:njFit
    seqPred = njFit+1:njPred
    #
    if (is.null(timeVar)) {
      dataFitj$time<- seqFit
      dataPredj$time<- seqPred
    }
    dataPlus <- rbind(dataFitj,dataPredj)
    #
    xPlus1 <- model.matrix(formFixed,data=dataPlus)
    zPlus1<-model.matrix(formRandom,data=dataPlus)
    z1 <- matrix(zPlus1[seqFit,],ncol=ncol(zPlus1))
    x1 <- matrix(xPlus1[seqFit,],ncol=ncol(xPlus1))
    z1Pred <- matrix(zPlus1[seqPred,],ncol=ncol(zPlus1))
    x1Pred <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
    #
    medFit <- x1%*%beta1 + c.*z1%*%Deltab
    medPred <- x1Pred%*%beta1 + c.*z1Pred%*%Deltab
    #
    y1=dataFitj[,all.vars(formFixed)[1]]
    SigmaPlus = sigmae*CovDEC(phiDEC,1,c(dataFitj$time,dataPredj$time))
    PsiPlus<-(zPlus1)%*%(D1)%*%t(zPlus1)+SigmaPlus
    dj<-as.numeric(t(y1-medFit)%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit))
    LambdaPlus <- solve(solve(D1)+ t(zPlus1)%*%solve(SigmaPlus)%*%zPlus1)
    sPsiPlus <- solve(PsiPlus)
    lambdaBarPlus <- matrix.sqrt(sPsiPlus)%*%zPlus1%*%D1%*%zeta/as.numeric(sqrt(1+t(zeta)%*%LambdaPlus%*%zeta))
    vj <- matrix.sqrt(sPsiPlus)%*%lambdaBarPlus
    Psi22.1 <- PsiPlus[seqPred,seqPred]- PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]
    vjtil <- (vj[seqFit] + solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]%*%vj[seqPred])/
      as.numeric(sqrt(1+t(vj[seqPred])%*%Psi22.1%*%vj[seqPred]))
    Ajj<-as.numeric(t(vjtil)%*%(y1-medFit)) #as.numeric(mutj/sqrt(Mtj2))
    mu2.1 <- medPred + PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit)
    if (distr=="sn") tau1 <- dnorm(Ajj)/pnorm(Ajj)
    if (distr=="st") {
      dtj = gamma((nu+njFit)/2)/gamma(nu/2)/pi^(njFit/2)/sqrt(det(as.matrix(PsiPlus[seqFit,seqFit])))*nu^(-njFit/2)*
        (dj/nu+1)^(-(njFit+nu)/2)
      denST = 2*dtj*pt(sqrt(nu+njFit)*Ajj/sqrt(dj+nu),nu+njFit)
      tau1 <- as.numeric(dmvt(y1,delta=medFit,sigma=as.matrix(PsiPlus[seqFit,seqFit]),df=nu,log=F)*gamma((nu+njFit-1)/2)*(nu+dj)^((nu+njFit)/2)/
                           (denST*pi^.5*gamma((nu+njFit)/2)*(nu+dj+Ajj^2)^((nu+njFit-1)/2)))
    }
    if (distr=="ss") {
      f2 <- function(u) u^(nu - 1)*((2*pi)^(-njFit/2))*(u^(njFit/2))*((det(as.matrix(PsiPlus[seqFit,seqFit])))^(-1/2))*
        exp(-0.5*u*dj)*pnorm(u^(1/2)*Ajj)
      denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
      tau1<-2^(nu)*nu*gamma(nu-.5+njFit/2)*pgamma(1,nu-.5+njFit/2,(dj+Ajj^2)/2)/
        (denSS*(dj+Ajj^2)^(nu-.5+njFit/2)*pi^(njFit/2+.5)*det(as.matrix(PsiPlus[seqFit,seqFit]))^.5)
    }
    if (distr=="scn") {
      fy<-as.numeric(2*(nu[1]*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                          (1-nu[1])*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]))*pnorm(Ajj,0,1)))
      tau1<-as.numeric(2*(nu[1]*nu[2]^(-1/2)*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]/nu[2]))*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                            (1-nu[1])*dmvnorm(y1,medFit,as.matrix(PsiPlus[seqFit,seqFit]))*dnorm(Ajj,0,1))/fy)
    }
    ypredj <- mu2.1 + tau1/as.numeric(sqrt(1+t(vj[seqPred])%*%Psi22.1%*%vj[seqPred]))*Psi22.1%*%vj[seqPred]
    ypred[dataPred$ind==indj] <- ypredj
    xpred[dataPred$ind==indj,] <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
    timepred[dataPred$ind==indj] <- dataPredj$time
  }
  xpred = as.data.frame(xpred)
  colnames(xpred) = colnames(xPlus1)
  if (all(xpred[,1]==1)) xpred=xpred[-1]
  if (is.null(timeVar)) dfout = data.frame(groupVar=dataPred$ind,xpred,ypred)
  else dfout = data.frame(groupVar=dataPred$ind,time=timepred,xpred,ypred)
  dfout
}

#Information matrix for SMSN-LMM and SMSN-LMM-AR(p) with E(bi)=0
IPhi <- function(w,di,Ai,distr,nu) {
  if (distr=="sn") intval <- exp(-.5*di)*pnorm(Ai)
  else if (distr=="st") intval<- 2^w*nu^(nu/2)*gamma(w+nu/2)*
      pt(Ai*sqrt((2*w+nu)/(di+nu)),2*w+nu)/((di+nu)^(w+nu/2)*gamma(nu/2))
  else if (distr=="ss") {
    U <- runif(5000)
    V <- pgamma(1,w+nu, di/2)*U
    S <- qgamma(V,w+nu, di/2)
    intval<- nu*2^(nu+w)*gamma(nu+w)*pgamma(1,nu+w,di/2)*mean(pnorm(S^(1/2)*Ai))/(di^(nu+w))
  }
  else if (distr=="scn") intval <- nu[1]*nu[2]^w*exp(-.5*nu[2]*di)*pnorm(nu[2]^.5*Ai)+
      (1-nu[1])*exp(-.5*di)*pnorm(Ai)
  return(intval)
}
Iphi <- function(w,di,Ai,distr,nu) {
  if (distr=="sn") intval <- exp(-.5*di)*dnorm(Ai)
  else if (distr=="st") intval<- 2^w*nu^(nu/2)*gamma(w+nu/2)/
      (sqrt(2*pi)*gamma(nu/2)*(di+Ai^2+nu)^(w+nu/2))
  else if (distr=="ss") intval<- nu*2^(nu+w)*gamma(nu+w)*pgamma(1,nu+w,(di+Ai^2)/2)/
      (sqrt(2*pi)*(di+Ai^2)^(nu+w))
  else if (distr=="scn") intval <- nu[1]*nu[2]^w*exp(-.5*nu[2]*di)*dnorm(nu[2]^.5*Ai)+
      (1-nu[1])*exp(-.5*di)*dnorm(Ai)
  return(intval)
}

F.r <- function(r,q1){
  Fmat. <- matrix(0,ncol=q1,nrow=q1)
  Fmat.[upper.tri(Fmat.,diag=T)][r] = 1
  Fmat.[lower.tri(Fmat.)] = t(Fmat.)[lower.tri(Fmat.)]
  Fmat.
}

autocovsAR <- function(phi,n) {
  p <- length(phi)
  if (n==1) Rn <- 1
  else Rn<- ARMAacf(ar=phi, ma=0, lag.max = n-1)
  rhos <- ARMAacf(ar=phi, ma=0, lag.max = p)[-1]
  return(Rn/(1-sum(rhos*phi)))
}

scorei <- function(jseq,y,x,z,beta1,sigmae,D1,lambda,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[jseq]
  p= ncol(x)
  q1=ncol(z)
  q2 = q1*(q1+1)/2
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  ni = length(y1)
  Fmat = matrix.sqrt(D1)
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-Fmat%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  Sigma <- sigmae*diag(ni)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  sPsi <- solve(Psi)
  di<-as.numeric(t(y1-med)%*%sPsi%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ai<-as.numeric(mutj/sqrt(Mtj2))
  sFmat = solve(Fmat)
  Lambda = solve(solve(D1)+ t(z1)%*%z1/sigmae)
  F.lista <- lapply(1:q2,F.r,q1=q1)
  #theta = c(beta1,sigmae,dd,lambda,nu) - para independente
  indpar = c(rep("beta",p),"sigma",rep("dd",q2),rep("lambda",q1))
  lpar = length(indpar)
  ##### derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)
  dlogdpsi[indpar=="sigma"] =traceM(sPsi)
  for (i in 1:q2) dlogdpsi[indpar=="dd"][i] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                    Fmat%*%F.lista[[i]])%*%t(z1))

  ##### derivadas de Ai
  dAi = numeric(lpar)
  ai = as.numeric((1+t(lambda)%*%sFmat%*%Lambda%*%sFmat%*%lambda)^.5)
  bi = as.numeric((1+t(lambda)%*%lambda)^.5)
  Bi = as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)
  dAi[indpar=="beta"] = -1/ai*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda
  dAi[indpar=="lambda"] = 1/ai*Fmat%*%t(z1)%*%sPsi%*%(y1-x1%*%beta1- 2*c.*z1%*%Deltab)-
    1/ai^2*Ai*sFmat%*%Lambda%*%sFmat%*%lambda + c.*Bi/ai/(bi^3)*lambda
  dAi[indpar=="sigma"] = -1/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae^2)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda
  for (i in 1:q2) dAi[indpar=="dd"][i] = 1/ai*(t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med) -
                                                 c./bi*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda)+
    1/ai^2*Ai/2*t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda

  ##### derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)
  ddi[indpar=="lambda"] = -2*c./bi*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="sigma"] = -t(y1-med)%*%sPsi%*%sPsi%*%(y1-med)
  for (i in 1:q2) ddi[indpar=="dd"][i] = -2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med) -
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)

  ##### derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  sihat = -.5*dlogdpsi+1/ki*dki
  sihat

}


Infmatrix <- function(y,x,z,ind,beta1,sigmae,D1,lambda,distr,nu,diagD,skewind){
  N <-length(y)
  p= ncol(x)
  q1=ncol(z)
  q2 = q1*(q1+1)/2
  #indpar = c(rep("beta",p),"sigma",rep("dd",q2),rep("lambda",q1))
  score_list=tapply(1:N,ind,scorei,y=y, x=x, z=z, beta1=beta1, sigmae=sigmae,D1=D1,lambda=lambda,distr=distr,nu=nu)
  lpar <-p+q2+q1+1 
  indplus <- rep(1,lpar)
  indplus[(p+q2+2):lpar] <- skewind
  if (diagD) indplus[(p+2):(p+q2+1)] <- diag(q1)[upper.tri(D1, diag = T)]
  mi_list = lapply(score_list,function(tt) {xm = matrix(tt[indplus==1],ncol=1);xm%*%t(xm)})
  infmat <- Reduce("+",mi_list)
  if (abs(det(infmat))<1e-5) infmat= infmat+1e-20*diag(nrow(infmat))
  sqrt(diag(solve(infmat)))
}

#########################################################################
#ar(p)
#########################################################################
scoreARi <- function(jseq,y,x,z,time,beta1,sigmae,phiAR,D1,lambda,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[jseq]
  t1 = time[jseq]
  p= ncol(x);q1=ncol(z);pAR=length(phiAR)
  q2 = q1*(q1+1)/2
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  ni = length(y1)
  Fmat = matrix.sqrt(D1)
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-Fmat%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  MniAR <- CovARp(phi = phiAR,t1)
  sMniAR<-solve(MniAR)
  Sigma <- sigmae*MniAR
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  sPsi <- solve(Psi)
  di<-as.numeric(t(y1-med)%*%sPsi%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ai<-as.numeric(mutj/sqrt(Mtj2))
  sFmat = solve(Fmat)
  Lambda = solve(solve(D1)+ t(z1)%*%solve(Sigma)%*%z1)
  F.lista <- lapply(1:q2,F.r,q1=q1)
  #theta = c(beta1,sigmae,phi,dd,lambda,nu) - para AR(p)
  indpar = c(rep("beta",p),"sigma",rep("phi",pAR),rep("dd",q2),rep("lambda",q1))
  lpar = length(indpar)
  ##### derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)
  dlogdpsi[indpar=="sigma"] =traceM(sPsi%*%MniAR)
  for (i in 1:q2) dlogdpsi[indpar=="dd"][i] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                    Fmat%*%F.lista[[i]])%*%t(z1))
  #jacobAR <- jacobian(Mnp,phiAR,n=ni) #matrix(jacobAR[,1],ncol=ni)
  jacobARautocovs <- matrix(jacobian(autocovsAR,phiAR,n=max(t1))[t1,],ncol=pAR) #toeplitz(jacobARautocovs[,1])
  for (i in 1:pAR) dlogdpsi[indpar=="phi"][i] = sigmae*traceM(sPsi%*%toeplitz(jacobARautocovs[,i]))

  ##### derivadas de Ai para diferente de nu
  dAi = numeric(lpar)
  ai = as.numeric((1+t(lambda)%*%sFmat%*%Lambda%*%sFmat%*%lambda)^.5)
  bi = as.numeric((1+t(lambda)%*%lambda)^.5)
  Bi = as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)
  dAi[indpar=="beta"] = -1/ai*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda
  dAi[indpar=="lambda"] = 1/ai*Fmat%*%t(z1)%*%sPsi%*%(y1-x1%*%beta1- 2*c.*z1%*%Deltab)-
    1/ai^2*Ai*sFmat%*%Lambda%*%sFmat%*%lambda + c.*Bi/ai/(bi^3)*lambda
  dAi[indpar=="sigma"] = -1/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae^2)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda

  for (i in 1:q2) dAi[indpar=="dd"][i] = 1/ai*(t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 c./bi*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda)+
    1/ai^2*Ai/2*t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                       Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:pAR) dAi[indpar=="phi"][i] = -sigmae/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[,i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda

  ##### derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)
  ddi[indpar=="lambda"] = -2*c./bi*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="sigma"] = -t(y1-med)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)
  for (i in 1:q2) ddi[indpar=="dd"][i] =-2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)
  for (i in 1:pAR) ddi[indpar=="phi"][i] = -sigmae*t(y1-med)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%(y1-med)

  ##### derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = numeric(lpar)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  sihat = -.5*dlogdpsi+1/ki*dki
  sihat

}


InfmatrixAR <- function(y,x,z,time,ind,beta1,sigmae,phiAR,D1,lambda,distr,nu,
                        diagD,skewind){
  N <-length(y)
  p= ncol(x)
  q1=ncol(z)
  q2 = q1*(q1+1)/2
  score_list=tapply(1:N,ind,scoreARi,y=y, x=x, z=z,time=time, beta1=beta1, sigmae=sigmae,phiAR=phiAR,D1=D1,lambda=lambda,distr=distr,nu=nu)
  #indpar = c(rep("beta",p),"sigma",rep("phi",pAR),rep("dd",q2),rep("lambda",q1))
  lpar <-p+1+length(phiAR)+q2+q1
  indplus <- rep(1,lpar)
  indplus[(p+1+length(phiAR)+q2+1):lpar] <- skewind
  if (diagD) indplus[(p+2+length(phiAR)):(p+q2+1+length(phiAR))] <- diag(q1)[upper.tri(D1, diag = T)]
  mi_list = lapply(score_list,function(tt) {xm = matrix(tt[indplus==1],ncol=1);xm%*%t(xm)})
  infmat <- Reduce("+",mi_list)
  if (abs(det(infmat))<1e-5) infmat= infmat+1e-10*diag(nrow(infmat))
  sqrt(diag(solve(infmat)))
}

#########################################################################
#CS
#########################################################################
scoreCSi <- function(jseq,y,x,z,beta1,sigmae,phiCS,D1,lambda,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[jseq]
  p= ncol(x);q1=ncol(z)
  q2 = q1*(q1+1)/2
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  ni = length(y1)
  Fmat = matrix.sqrt(D1)
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-Fmat%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  Covmat <- CovCS(phiCS,ni)
  sCovmat<-solve(Covmat)
  Sigma <- sigmae*Covmat
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  sPsi <- solve(Psi)
  di<-as.numeric(t(y1-med)%*%sPsi%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ai<-as.numeric(mutj/sqrt(Mtj2))
  sFmat = solve(Fmat)
  Lambda = solve(solve(D1)+ t(z1)%*%solve(Sigma)%*%z1)
  F.lista <- lapply(1:q2,F.r,q1=q1)
  #theta = c(beta1,sigmae,phi,dd,lambda,nu) - para AR(p)
  indpar = c(rep("beta",p),"sigma","phi",rep("dd",q2),rep("lambda",q1))
  lpar = length(indpar)
  ##### derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)
  dlogdpsi[indpar=="sigma"] =traceM(sPsi%*%Covmat)
  for (i in 1:q2) dlogdpsi[indpar=="dd"][i] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                    Fmat%*%F.lista[[i]])%*%t(z1))
  derCov <- matrix(1,nrow=ni,ncol=ni)-diag(ni)
  dlogdpsi[indpar=="phi"] = sigmae*traceM(sPsi%*%derCov)

  ##### derivadas de Ai para diferente de nu
  dAi = numeric(lpar)
  ai = as.numeric((1+t(lambda)%*%sFmat%*%Lambda%*%sFmat%*%lambda)^.5)
  bi = as.numeric((1+t(lambda)%*%lambda)^.5)
  Bi = as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)
  dAi[indpar=="beta"] = -1/ai*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda
  dAi[indpar=="lambda"] = 1/ai*Fmat%*%t(z1)%*%sPsi%*%(y1-x1%*%beta1- 2*c.*z1%*%Deltab)-
    1/ai^2*Ai*sFmat%*%Lambda%*%sFmat%*%lambda + c.*Bi/ai/(bi^3)*lambda
  dAi[indpar=="sigma"] = -1/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae^2)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  for (i in 1:q2) dAi[indpar=="dd"][i] = 1/ai*(t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 c./bi*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda)+
    1/ai^2*Ai/2*t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                       Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  dAi[indpar=="phi"] = -sigmae/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%derCov%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%derCov%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  ##### derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)
  ddi[indpar=="lambda"] = -2*c./bi*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="sigma"] = -t(y1-med)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)
  for (i in 1:q2) ddi[indpar=="dd"][i] =-2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="phi"] = -sigmae*t(y1-med)%*%sPsi%*%derCov%*%sPsi%*%(y1-med)

  ##### derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = numeric(lpar)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  sihat = -.5*dlogdpsi+1/ki*dki
  sihat

}


InfmatrixCS <- function(y,x,z,ind,beta1,sigmae,phiCS,D1,lambda,distr,nu,diagD,skewind){
  N <-length(y)
  p= ncol(x)
  q1=ncol(z)
  q2 = q1*(q1+1)/2
  score_list=tapply(1:N,ind,scoreCSi,y=y, x=x, z=z, beta1=beta1, sigmae=sigmae,phiCS=phiCS,D1=D1,lambda=lambda,distr=distr,nu=nu)
  #indpar = c(rep("beta",p),"sigma","phi",rep("dd",q2),rep("lambda",q1))
  lpar <-p+2+q2+q1
  indplus <- rep(1,lpar)
  indplus[(p+2+q2+1):lpar] <- skewind
  if (diagD) indplus[(p+3):(p+q2+2)] <- diag(q1)[upper.tri(D1, diag = T)]
  mi_list = lapply(score_list,function(tt) {xm = matrix(tt[indplus==1],ncol=1);xm%*%t(xm)})
  infmat <- Reduce("+",mi_list)
  if (abs(det(infmat))<2e-5) infmat= infmat+1e-10*diag(nrow(infmat))
  sqrt(diag(solve(infmat)))
}

#########################################################################
#DEC
#########################################################################
dphiCovDEC<-function(phi,theta=1,ti) {
  ni <- length(ti)
  Rn <- matrix(0,nrow=ni,ncol=ni)
  if (ni==1) Rn <- 0
  else {
    for (i in 1:(ni-1)) for (j in (i+1):ni) Rn[i,j] <- abs(ti[i]-ti[j])^theta*phi^(abs(ti[i]-ti[j])^theta-1)
    Rn[lower.tri(Rn)] <-  t(Rn)[lower.tri(Rn)]
  }
  return(Rn)
}
dthetaCovDEC<-function(phi,theta,ti) {
  ni <- length(ti)
  Rn <- matrix(0,nrow=ni,ncol=ni)
  if (ni==1) Rn <- 0
  else {
    for (i in 1:(ni-1)) for (j in (i+1):ni) Rn[i,j] <- abs(ti[i]-ti[j])^theta*
        phi^(abs(ti[i]-ti[j])^theta)*log(phi)*log(abs(ti[i]-ti[j]))
    Rn[lower.tri(Rn)] <-  t(Rn)[lower.tri(Rn)]
  }
  return(Rn)
}

scoreDECi <- function(jseq,y,x,z,time,beta1,sigmae,phiDEC,thetaDEC,D1,lambda,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[jseq]
  t1 = time[jseq]
  p= ncol(x);q1=ncol(z)
  q2 = q1*(q1+1)/2
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  ni = length(y1)
  Fmat = matrix.sqrt(D1)
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-Fmat%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  Covmat <- CovDEC(phiDEC,thetaDEC,t1)
  sCovmat<-solve(Covmat)
  Sigma <- sigmae*Covmat
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  sPsi <- solve(Psi)
  di<-as.numeric(t(y1-med)%*%sPsi%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ai<-as.numeric(mutj/sqrt(Mtj2))
  sFmat = solve(Fmat)
  Lambda = solve(solve(D1)+ t(z1)%*%solve(Sigma)%*%z1)
  F.lista <- lapply(1:q2,F.r,q1=q1)
  #theta = c(beta1,sigmae,phi,dd,lambda,nu) - para AR(p)
  indpar = c(rep("beta",p),"sigma","phi","theta",rep("dd",q2),rep("lambda",q1))
  lpar = length(indpar)
  ##### derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)
  dlogdpsi[indpar=="sigma"] =traceM(sPsi%*%Covmat)
  for (i in 1:q2) dlogdpsi[indpar=="dd"][i] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                    Fmat%*%F.lista[[i]])%*%t(z1))
  dphiDEC <- dphiCovDEC(phiDEC,thetaDEC,t1)
  dthetaDEC <- dthetaCovDEC(phiDEC,thetaDEC,t1)
  dlogdpsi[indpar=="phi"] = sigmae*traceM(sPsi%*%dphiDEC)
  dlogdpsi[indpar=="theta"] = sigmae*traceM(sPsi%*%dthetaDEC)

  ##### derivadas de Ai para diferente de nu
  dAi = numeric(lpar)
  ai = as.numeric((1+t(lambda)%*%sFmat%*%Lambda%*%sFmat%*%lambda)^.5)
  bi = as.numeric((1+t(lambda)%*%lambda)^.5)
  Bi = as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)
  dAi[indpar=="beta"] = -1/ai*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda
  dAi[indpar=="lambda"] = 1/ai*Fmat%*%t(z1)%*%sPsi%*%(y1-x1%*%beta1- 2*c.*z1%*%Deltab)-
    1/ai^2*Ai*sFmat%*%Lambda%*%sFmat%*%lambda + c.*Bi/ai/(bi^3)*lambda
  dAi[indpar=="sigma"] = -1/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae^2)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  for (i in 1:q2) dAi[indpar=="dd"][i] = 1/ai*(t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 c./bi*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda)+
    1/ai^2*Ai/2*t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                       Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  dAi[indpar=="phi"] = -sigmae/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda
  dAi[indpar=="theta"] = -sigmae/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  ##### derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)
  ddi[indpar=="lambda"] = -2*c./bi*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="sigma"] = -t(y1-med)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)
  for (i in 1:q2) ddi[indpar=="dd"][i] =-2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="phi"] = -sigmae*t(y1-med)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)
  ddi[indpar=="theta"] = -sigmae*t(y1-med)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)

  ##### derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = numeric(lpar)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  sihat = -.5*dlogdpsi+1/ki*dki
  sihat

}


InfmatrixDEC <- function(y,x,z,time,ind,beta1,sigmae,phiDEC,thetaDEC,D1,lambda,distr,
                         nu,diagD,skewind){
  N <-length(y)
  p= ncol(x)
  q1=ncol(z)
  q2 = q1*(q1+1)/2
  score_list=tapply(1:N,ind,scoreDECi,y=y, x=x, z=z,time=time, beta1=beta1, sigmae=sigmae,
                    phiDEC=phiDEC,thetaDEC=thetaDEC,D1=D1,lambda=lambda,distr=distr,nu=nu)
  #indpar = c(rep("beta",p),"sigma","phi","theta",rep("dd",q2),rep("lambda",q1))
  lpar <-p+3+q2+q1
  indplus <- rep(1,lpar)
  indplus[(p+3+q2+1):lpar] <- skewind
  if (diagD) indplus[(p+4):(p+q2+3)] <- diag(q1)[upper.tri(D1, diag = T)]
  mi_list = lapply(score_list,function(tt) {xm = matrix(tt[indplus==1],ncol=1);xm%*%t(xm)})
  infmat <- Reduce("+",mi_list)
  if (abs(det(infmat))<1e-5) infmat= infmat+1e-10*diag(nrow(infmat))
  sqrt(diag(solve(infmat)))
}

#########################################################################
#CAR(1)
#########################################################################
scoreCAR1i <- function(jseq,y,x,z,time,beta1,sigmae,phiDEC,D1,lambda,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[jseq]
  t1 = time[jseq]
  p= ncol(x);q1=ncol(z)
  q2 = q1*(q1+1)/2
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  ni = length(y1)
  Fmat = matrix.sqrt(D1)
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-Fmat%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  Covmat <- CovDEC(phiDEC,1,t1)
  sCovmat<-solve(Covmat)
  Sigma <- sigmae*Covmat
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  sPsi <- solve(Psi)
  di<-as.numeric(t(y1-med)%*%sPsi%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ai<-as.numeric(mutj/sqrt(Mtj2))
  sFmat = solve(Fmat)
  Lambda = solve(solve(D1)+ t(z1)%*%solve(Sigma)%*%z1)
  F.lista <- lapply(1:q2,F.r,q1=q1)
  #theta = c(beta1,sigmae,phi,dd,lambda,nu) - para AR(p)
  indpar = c(rep("beta",p),"sigma","phi",rep("dd",q2),rep("lambda",q1))
  lpar = length(indpar)
  ##### derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)
  dlogdpsi[indpar=="sigma"] =traceM(sPsi%*%Covmat)
  for (i in 1:q2) dlogdpsi[indpar=="dd"][i] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                    Fmat%*%F.lista[[i]])%*%t(z1))
  dphiDEC <- dphiCovDEC(phiDEC,1,t1)
  dlogdpsi[indpar=="phi"] = sigmae*traceM(sPsi%*%dphiDEC)

  ##### derivadas de Ai para diferente de nu
  dAi = numeric(lpar)
  ai = as.numeric((1+t(lambda)%*%sFmat%*%Lambda%*%sFmat%*%lambda)^.5)
  bi = as.numeric((1+t(lambda)%*%lambda)^.5)
  Bi = as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)
  dAi[indpar=="beta"] = -1/ai*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda
  dAi[indpar=="lambda"] = 1/ai*Fmat%*%t(z1)%*%sPsi%*%(y1-x1%*%beta1- 2*c.*z1%*%Deltab)-
    1/ai^2*Ai*sFmat%*%Lambda%*%sFmat%*%lambda + c.*Bi/ai/(bi^3)*lambda
  dAi[indpar=="sigma"] = -1/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae^2)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  for (i in 1:q2) dAi[indpar=="dd"][i] = 1/ai*(t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 c./bi*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda)+
    1/ai^2*Ai/2*t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                       Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  dAi[indpar=="phi"] = -sigmae/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  ##### derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)
  ddi[indpar=="lambda"] = -2*c./bi*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="sigma"] = -t(y1-med)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)
  for (i in 1:q2) ddi[indpar=="dd"][i] =-2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="phi"] = -sigmae*t(y1-med)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)

  ##### derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = numeric(lpar)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  sihat = -.5*dlogdpsi+1/ki*dki
  sihat

}


InfmatrixCAR1 <- function(y,x,z,time,ind,beta1,sigmae,phiDEC,D1,lambda,distr,nu,
                          diagD,skewind){
  N <-length(y)
  p= ncol(x)
  q1=ncol(z)
  q2 = q1*(q1+1)/2
  score_list=tapply(1:N,ind,scoreCAR1i,y=y, x=x, z=z,time=time, beta1=beta1, sigmae=sigmae,
                    phiDEC=phiDEC,D1=D1,lambda=lambda,distr=distr,nu=nu)
  #indpar = c(rep("beta",p),"sigma","phi",rep("dd",q2),rep("lambda",q1))
  lpar <-p+2+q2+q1
  indplus <- rep(1,lpar)
  indplus[(p+2+q2+1):lpar] <- skewind
  if (diagD) indplus[(p+3):(p+q2+2)] <- diag(q1)[upper.tri(D1, diag = T)]
  mi_list = lapply(score_list,function(tt) {xm = matrix(tt[indplus==1],ncol=1);xm%*%t(xm)})
  infmat <- Reduce("+",mi_list)
  if (abs(det(infmat))<1e-5) infmat= infmat+1e-10*diag(nrow(infmat))
  sqrt(diag(solve(infmat)))
}

