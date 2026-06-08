# Sandwich Variance Estimator --------------------------------------------------
# first and second derivatives of likelihood

# Auxiliary  --------------------------------------------------------------------------

gerar_smsn = function(jvec,x,z,sigma2,Dsqrti,beta1,lambda,distr,nu,
                      ind,time,depStruct,phi) {
  if (distr=="sn"||distr=="norm") {ui=1; c.=-sqrt(2/pi)}
  if (distr=="st"||distr=="t") {ui=rgamma(1,nu/2,nu/2); c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss"||distr=="ssl"||distr=="sl") {ui=rbeta(1,nu,1); c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn"||distr=="cn") {ui=ifelse(runif(1)<nu[1],nu[2],1);
  c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jvec,  ],ncol=p)
  z1=matrix(z[jvec,  ],ncol=q1)
  nj = nrow(x1)
  Sig = errorVar(time[jvec],sigma2=sigma2,depStruct=depStruct,phi=phi)#sigma2*diag(nj)
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%(lambda)))
  Delta = Dsqrti%*%delta
  Gammab = Dsqrti%*%Dsqrti - Delta%*%t(Delta)
  Beta = matrix(beta1,ncol=1)
  ti = c.+abs(rnorm(1,0,ui^-.5))
  bi = t(rmvnorm(1,Delta*ti,sigma=ui^(-1)*Gammab))
  Yi = t(rmvnorm(1,x1%*%Beta+z1%*%bi,sigma=ui^(-1)*Sig))
  if (all(x1[,1]==1)) x1 = x1[,-1]
  if (all(z1[,1]==1)) z1 = z1[,-1]
  return(data.frame(y=Yi,ind=ind[jvec],time=time[jvec]))
}

F.r <- function(r,q1){
  Fmat. <- matrix(0,ncol=q1,nrow=q1)
  Fmat.[upper.tri(Fmat.,diag=T)][r] = 1
  Fmat.[lower.tri(Fmat.)] = t(Fmat.)[lower.tri(Fmat.)]
  Fmat.
}

autocovsAR2 <- function(phi,n1,j) {
  p <- length(phi)
  if (n1==1) Rn <- 1
  else Rn<- ARMAacf(ar=phi, ma=0, lag.max = n1-1)[j]
  rhos <- ARMAacf(ar=phi, ma=0, lag.max = p)[-1]

  return(Rn/(1-sum(rhos*phi)))
}

selectHessian <- function(k, l, n1, lista){
  aux <- 0
  for(s in 1:n1){
    aux[s] <- lista[[s]][k,l]
  }
  return(aux)
}

# DEC - derivadas segundas de Ri em relacao a phi_1 e phi_2
dphiphiCovDEC <- function(phi, theta, ti){
  ni <- length(ti)
  Rn <- matrix(0, nrow=ni, ncol=ni)
  if (ni == 1) Rn <- 0
  else {
    for (i in 1:(ni-1)) for (j in (i+1):ni) Rn[i,j] <- abs(ti[i]-ti[j])^theta*(abs(ti[i]-ti[j])^theta - 1)*phi^(abs(ti[i]-ti[j])^theta-2)
    Rn[lower.tri(Rn)] <-  t(Rn)[lower.tri(Rn)]
  }
  return(Rn)
}

dthetathetaCovDEC <- function(phi, theta, ti){
  ni <- length(ti)
  Rn <- matrix(0,nrow=ni,ncol=ni)
  if (ni == 1) Rn <- 0
  else {
    for (i in 1:(ni-1)) for (j in (i+1):ni) Rn[i,j] <- abs(ti[i]-ti[j])^theta*log(abs(ti[i]-ti[j]))^2*phi^(abs(ti[i]-ti[j])^theta)*log(phi)*(abs(ti[i]-ti[j])^theta*log(phi) + 1)
    Rn[lower.tri(Rn)] <-  t(Rn)[lower.tri(Rn)]
  }
  return(Rn)
}

dphithetaCovDEC <- function(phi, theta, ti){
  ni <- length(ti)
  Rn <- matrix(0, nrow=ni, ncol=ni)
  if (ni==1) Rn <- 0
  else {
    for (i in 1:(ni-1)) for (j in (i+1):ni) Rn[i,j] <- abs(ti[i]-ti[j])^theta*phi^(abs(ti[i]-ti[j])^theta-1)*log(abs(ti[i]-ti[j]))*(abs(ti[i]-ti[j])^theta*log(phi) + 1)
    Rn[lower.tri(Rn)] <-  t(Rn)[lower.tri(Rn)]
  }
  return(Rn)
}

# SMSN - AR(p) -----------------------------------------------------------------
derivatesARi <- function(jseq,y,x,z,time,beta1,sigmae,phiAR,D1,lambda,distr,nu){
  if (distr=="sn"|distr=="norm") c.=-sqrt(2/pi)
  if (distr=="st"|distr=="t") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ssl"|distr=="sl") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn"|distr=="cn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[jseq]
  t1 = time[jseq]
  p= ncol(x);q1=ncol(z);pAR=length(phiAR)
  q2 = q1*(q1+1)/2
  if(q1 != 0){ q3 = 1} else{q3 = 0}
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
  # theta = c(beta1,sigmae,phi,dd,lambda,nu) - para AR(p)
  indpar = c(rep("beta",p),"sigma",rep("phi",pAR),rep("dd",q2),rep("lambda",q1))
  lpar = length(indpar)
  indpar2 = matrix(0, nrow = (p+1+pAR+q2+q1), ncol = (p+1+pAR+q2+q1))

  ##### Primeiras derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)
  dlogdpsi[indpar=="sigma"] =traceM(sPsi%*%MniAR)
  for (i in 1:q2) dlogdpsi[indpar=="dd"][i] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                    Fmat%*%F.lista[[i]])%*%t(z1))
  #jacobAR <- jacobian(Mnp,phiAR,n=ni) #matrix(jacobAR[,1],ncol=ni)
  jacobARautocovs <- matrix(jacobian(autocovsAR,phiAR,n=max(t1))[t1,],ncol=pAR) #toeplitz(jacobARautocovs[,1])
  for (i in 1:pAR) dlogdpsi[indpar=="phi"][i] = sigmae*traceM(sPsi%*%toeplitz(jacobARautocovs[,i]))

  ##### Primeiras derivadas de Ai para diferente de nu
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

  ##### Primeiras derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)
  ddi[indpar=="lambda"] = -2*c./bi*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="sigma"] = -t(y1-med)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)
  for (i in 1:q2) ddi[indpar=="dd"][i] =-2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)
  for (i in 1:pAR) ddi[indpar=="phi"][i] = -sigmae*t(y1-med)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%(y1-med)

  ##### Primeiras derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = numeric(lpar)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  ##### Segundas Derivadas de log(det(psi))
  ddlogdpsi = indpar2
  # dsigmae.dsigmae
  ddlogdpsi[p+1, p+1] = traceM(-sPsi%*%MniAR%*%sPsi%*%MniAR)
  # dsigmae.dphi
  for (i in 1:pAR) ddlogdpsi[p+1,p+1+i]= ddlogdpsi[p+1+i,p+1] = traceM(sPsi%*%toeplitz(jacobARautocovs[,i]) -
                                                                         sigmae*sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%MniAR)
  # dsigmae.dalphar
  for (i in 1:q2) ddlogdpsi[p+1,p+1+pAR+i] = ddlogdpsi[p+1+pAR+i,p+1] = traceM(-sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                                               Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%MniAR)
  # dalpha.dalpha
  for (i in 1:q2) for (j in 1:q2) ddlogdpsi[p+1+pAR+i,p+1+pAR+j] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%F.lista[[j]] + F.lista[[j]]%*%F.lista[[i]])%*%t(z1)-
                                                                            sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%
                                                                            sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1))
  # dalpha.dphi
  for (i in 1:q2) for (j in 1:pAR) ddlogdpsi[p+1+pAR+i, p+1+j] = ddlogdpsi[p+1+j, p+1+pAR+i] = - traceM(sigmae*sPsi%*%toeplitz(jacobARautocovs[,j])%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                                                                                                                     Fmat%*%F.lista[[i]])%*%t(z1))
  # dphi.dphi
  n1 <- length(y1)
  hessARautocovs <- list() # Derivadas dos phis em cada lag
  for(k in 1:n1){ hessARautocovs[[k]] <- hessian(autocovsAR2, phiAR, n=max(t1), j = k)}

  for (i in 1:pAR) for (j in 1:pAR) ddlogdpsi[p+1+i,p+1+j] = sigmae*traceM(sPsi%*%toeplitz(selectHessian(i, j, n1, hessARautocovs)) -
                                                                             sigmae*sPsi%*%toeplitz(jacobARautocovs[, j])%*%sPsi%*%toeplitz(jacobARautocovs[, i]))
  ##### Segundas Derivadas de Ai
  ddAi <- indpar2
  # dbeta.dbeta = 0

  # dbeta.dsigmae
  ddAi[1:p, p+1] <- (1/ai)*t(x1)%*%sPsi%*%MniAR%*%sPsi%*%z1%*%Fmat%*%lambda+
    (1/(2*ai^3*sigmae^2))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda
  ddAi[p+1, 1:p] <- t(ddAi[1:p, p+1])
  # dbeta.dphi
  for (i in 1:pAR) ddAi[1:p, p+1+i] <- (sigmae/ai)*t(x1)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%z1%*%Fmat%*%lambda +
    (1/(2*ai^3*sigmae))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[,i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda
  for (i in 1:pAR) ddAi[p+1+i, 1:p] <- ddAi[1:p, p+1+i]

  #dbeta.dalphar
  for (i in 1:q2) ddAi[1:p, p+1+pAR+i] <- (1/ai)*t(x1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda-
    (1/ai)*t(x1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda - (1/(2*ai^3))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                                               Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat + sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:q2) ddAi[p+1+pAR+i,1:p] <- ddAi[1:p, p+1+pAR+i]
  # dbeta.dlambda
  ddAi[1:p, (p+1+pAR+q2+q3):(p+1+pAR+q2+q1)] <- (-1/ai)*t(x1)%*%sPsi%*%z1%*%Fmat +(1/ai^3)*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat
  ddAi[(p+1+pAR+q2+q3):(p+1+pAR+q2+q1), 1:p] <- ddAi[1:p, (p+1+pAR+q2+q3):(p+1+pAR+q2+q1)]

  # dsigma.dsigma
  ddAi[p+1,p+1] <- (2/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)+
    (1/(2*sigmae^2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="sigma"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda -
    (1/(sigmae^4*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(sigmae^3*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda +
    (1/(2*sigmae^4*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda

  #dsigma.dphi
  for (i in 1:pAR) ddAi[p+1,p+1+i] <- (-1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%(y1-med)+
    (2*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)+
    (1/(2*ai^3*sigmae))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[,i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="phi"][i]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(sigmae^3*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[,i])%*%sMniAR%*%z1%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[,i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^3*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[,i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda
  for (i in 1:pAR) ddAi[p+1+i,p+1] <- ddAi[p+1,p+1+i]

  #dsigma.dalphar
  for (i in 1:q2) ddAi[p+1, p+1+pAR+i] <- (-1/ai)*t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%z1%*%F.lista[[i]]%*%delta-
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                    Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+ sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="dd"][i]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%F.lista[[i]]%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%(sFmat%*%F.lista[[i]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[i]]%*%sFmat)%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%(sFmat%*%F.lista[[i]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[i]]%*%sFmat)%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%F.lista[[i]]%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^4))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                                                 Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+ sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:q2) ddAi[p+1+pAR+i, p+1] <- ddAi[p+1, p+1+pAR+i]

  #dsigmae.dlambda
  ddAi[p+1, (p+1+pAR+q2+q3):(p+1+pAR+q2+q1)] <- (-1/ai)*t(Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med))+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))+
    (1/ai^3)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat -
    (1/(2*sigmae^2*ai^2))*(t(dAi[indpar=="lambda"])*as.numeric(t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda)+
                             2*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat)+
    (1/(sigmae^2*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat

  ddAi[(p+1+pAR+q2+q3):(p+1+pAR+q2+q1), p+1] <- ddAi[p+1,(p+1+pAR+q2+q3):(p+1+pAR+q2+q1)] #ddAi[p+1+pAR+q2+q1, p+1] <- ddAi[p+1, p+1+pAR+q2+q1]
  #dphi.dphi
  for (i in 1:pAR) for (j in 1:pAR) ddAi[p+1+i, p+1+j] <- (-sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%
    (-sigmae*sPsi%*%toeplitz(jacobARautocovs[, j])%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi +
       sPsi%*%toeplitz(selectHessian(i, j, n1, hessARautocovs))%*%sPsi -
       sigmae*sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%toeplitz(jacobARautocovs[, j])%*%sPsi)%*%(y1-med)+
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, j])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*dAi[indpar=="phi"][j]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, j])%*%sMniAR%*%z1%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%(-sMniAR%*%toeplitz(jacobARautocovs[, j])%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR+
                                                                     sMniAR%*%toeplitz(selectHessian(i, j, n1, hessARautocovs))%*%sMniAR -
                                                                     sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%toeplitz(jacobARautocovs[, j])%*%sMniAR)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, j])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*ai^4*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, j])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda

  #dphi.dalphar
  for (i in 1:pAR) for(j in 1:q2) ddAi[p+1+i, p+1+pAR+j] <- (-sigmae/ai)*t(lambda)%*%F.lista[[j]]%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%(y1-med)+
    (sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[j]] + F.lista[[j]]%*%Fmat)%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%(y1-med)+
    (sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%z1%*%(Fmat%*%F.lista[[j]] + F.lista[[j]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c.*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%z1%*%F.lista[[j]]%*%delta-
    (sigmae/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[j]]%*%sFmat%*%Lambda+
                                                                                                                                  Lambda%*%sFmat%*%F.lista[[j]]-
                                                                                                                                  Lambda%*%sFmat%*%(F.lista[[j]]%*%sFmat+
                                                                                                                                                      sFmat%*%F.lista[[j]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*dAi[indpar=="dd"][j]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%F.lista[[j]]%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%(sFmat%*%F.lista[[j]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[j]]%*%sFmat)%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%(sFmat%*%F.lista[[j]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[j]]%*%sFmat)%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%F.lista[[j]]%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%(F.lista[[j]]%*%sFmat%*%Lambda+
                                                                                                                                                                   Lambda%*%sFmat%*%F.lista[[j]]-
                                                                                                                                                                   Lambda%*%sFmat%*%(F.lista[[j]]%*%sFmat+
                                                                                                                                                                                       sFmat%*%F.lista[[j]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:pAR) for(j in 1:q2) ddAi[p+1+pAR+ j, p+1+i] <- ddAi[p+1+i, p+1+pAR+j]

  #dphi.dlambda
  for (i in 1:pAR) ddAi[p+1+i, (p+1+pAR+q2+q3):(p+1+pAR+q2+q1)] <- (-sigmae/ai)*t(Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%(y1-med))+
    (c.*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))+
    (sigmae/ai^3)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat-
    (1/(2*sigmae*ai^2))*t(dAi[indpar=="lambda"])*as.numeric(t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda)-
    (1/(sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat+
    (1/(sigmae*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat

  for (i in 1:pAR) ddAi[(p+1+pAR+q2+q3):(p+1+pAR+q2+q1),p+1+i] <- ddAi[p+1+i, (p+1+pAR+q2+q3):(p+1+pAR+q2+q1)]

  #dalpha.dalpha
  for(r in 1:q2) for(s in 1:q2) ddAi[p+1+pAR+r, p+1+pAR+s] <- (-1/ai)*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]]+F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)-
    (c./ai)*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[s]]%*%delta+
    (1/(2*ai^3))*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                             Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (c./ai)*t(delta)%*%F.lista[[s]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda+
    (c./ai)*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda-
    (c./(2*ai^3))*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                 Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/ai)*t(lambda)%*%F.lista[[s]]%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)-
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[r]]%*%F.lista[[s]] + F.lista[[s]]%*%F.lista[[r]])%*%t(z1)%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[s]]%*%delta-
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                                                                       Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda+
    (1/(2*ai^2))*dAi[indpar=="dd"][s]%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                                               Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%F.lista[[s]]%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                                                    Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%F.lista[[s]]%*%sFmat%*%lambda+
    (1/(2*ai^4))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                                                                                      Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda+
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(-F.lista[[r]]%*%sFmat%*%F.lista[[s]]%*%sFmat%*%Lambda+
                                             F.lista[[r]]%*%sFmat%*%Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda+
                                             Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%F.lista[[s]]%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%F.lista[[s]]%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]]%*%sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%F.lista[[s]]%*%sFmat%*%Lambda-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda)%*%sFmat%*%lambda
  # dalpha.dlambda
  for(i in 1:q2) ddAi[p+1+pAR+i, (p+1+pAR+q2+q3):(p+1+pAR+q2+q1)] <- (1/ai)*t(F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med))-
    (c./ai)*t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))-
    (1/ai^3)*t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat-
    (c./ai)*t(t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda)-
    (c./ai)*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]+ #talvez t
    (c./ai^3)*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat-
    (1/ai)*t(Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med))+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))+
    (1/ai^3)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat+
    (1/(2*ai^2))*t(dAi[indpar=="lambda"])*as.numeric(t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                            Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda)+
    (1/(2*ai^2))*Ai%*%t(sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                   Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda)+
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                             Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat-
    (1/(ai^4))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                           Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat

  for(i in 1:q2) ddAi[(p+1+pAR+q2+q3):(p+1+pAR+q2+q1),p+1+pAR+i] <- ddAi[p+1+pAR+i, (p+1+pAR+q2+q3):(p+1+pAR+q2+q1)]

  # dlambda.dlambda
  ddAi[(p+1+pAR+q2+q3):(p+1+pAR+q2+q1), (p+1+pAR+q2+q3):(p+1+pAR+q2+q1)] <- (-2*c./ai)*Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))-
    (1/ai^3)*Fmat%*%t(z1)%*%sPsi%*%(y1-x1%*%beta1 - 2*c.*z1%*%Deltab)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat+
    as.numeric(2*c./(ai*(1+t(lambda)%*%lambda)^1.5))*Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)+
    as.numeric(c./(ai*(1+t(lambda)%*%lambda)))*as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)*t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))-
    as.numeric(c./(ai^2*(1+t(lambda)%*%lambda)^2.5))*as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)*lambda%*%((1/ai)*t(lambda)%*%sFmat%*%Lambda%*%sFmat*as.numeric((1+t(lambda)%*%lambda))+2*ai*t(lambda))-
    (1/ai^2)*dAi[indpar=="lambda"]%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat - (1/ai^2)*Ai*sFmat%*%Lambda%*%sFmat + (2/ai^4)*Ai*sFmat%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat

  ##### Segundas Derivadas de di
  dddi <- indpar2
  # dbeta.dbeta
  dddi[1:p,1:p] <- 2*t(x1)%*%sPsi%*%x1
  # dbeta.dsigmae
  dddi[1:p, p+1] <- dddi[p+1, 1:p] <- 2*t(x1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)
  # dbeta.dphi
  for(i in 1:pAR) dddi[1:p, p+1+i]  <- dddi[p+1+i, 1:p] <- 2*sigmae*t(x1)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%(y1-med)
  # dbeta.dalphar
  for(i in 1:q2) dddi[1:p, p+1+pAR+i] <- dddi[p+1+pAR+i,1:p] <- 2*t(x1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)+ 2*c.*t(x1)%*%sPsi%*%z1%*%F.lista[[i]]%*%delta
  # dbeta.dlambda
  dddi[1:p, (p+1+pAR+q2+q3):(p+1+pAR+q2+q1)] <- dddi[(p+1+pAR+q2+q3):(p+1+pAR+q2+q1),1:p] <- 2*c.*t(x1)%*%sPsi%*%z1%*%Fmat%*%t(1/((1+ as.numeric(t(lambda)%*%lambda))^.5)*diag(q1) - lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)
  # dsigma.dsigma
  dddi[p+1,p+1] <- 2*t(y1-med)%*%sPsi%*%MniAR%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)
  # dsigma.dphi
  for(i in 1:pAR) dddi[p+1,p+1+i] <- dddi[p+1+i,p+1] <-t(y1-med)%*%(sigmae*sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%MniAR%*%sPsi - sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi + sigmae*sPsi%*%MniAR%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi)%*%(y1-med)
  # dsigmae.dalpha
  for(i in 1:q2) dddi[p+1,p+1+pAR+i] <- dddi[p+1+pAR+i, p+1] <- c.*t(z1%*%F.lista[[i]]%*%delta)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med) + t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)+
    t(y1-med)%*%sPsi%*%MniAR%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med) + c.*t(y1-med)%*%sPsi%*%MniAR%*%sPsi%*%z1%*%F.lista[[i]]%*%delta
  # dsigmae.dlambda
  dddi[p+1,(p+1+pAR+q2+q3):(p+1+pAR+q2+q1)] <- dddi[(p+1+pAR+q2+q3):(p+1+pAR+q2+q1), p+1]<- (2*c./(1+as.numeric(t(lambda)%*%lambda))^.5)*(Fmat - delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)
  # dphi.dphi
  for(i in 1:pAR) for(j in 1:pAR) dddi[p+1+i,p+1+j] <- dddi[p+1+j,p+1+i] <- sigmae*t(y1-med)%*%(sigmae*sPsi%*%toeplitz(jacobARautocovs[,j])%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi - sPsi%*%toeplitz(selectHessian(i, j, n1, hessARautocovs))%*%sPsi + sigmae*sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%toeplitz(jacobARautocovs[,j])%*%sPsi)%*%(y1-med)
  # dphi.dalpha
  for(i in 1:pAR) for(j in 1:q2) dddi[p+1+i, p+1+pAR+j] <- dddi[p+1+pAR+j, p+1+i] <- 2*c.*sigmae*t(y1-med)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%z1%*%F.lista[[j]]%*%delta+
    2*sigmae%*%t(y1-med)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat + Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)
  # dphi. dlambda
  for(i in 1:pAR) dddi[p+1+i, (p+1+pAR+q2+q3):(p+1+pAR+q2+q1)] <- dddi[(p+1+pAR+q2+q3):(p+1+pAR+q2+q1), p+1+i] <- 2*c.*sigmae*t(y1-med)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%z1%*%Fmat%*%t(1/(1+as.numeric(t(lambda)%*%lambda))^.5 - lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)
  # dalpha.dalpha
  for(i in 1:q2) for(j in 1:q2) dddi[p+1+pAR+i, p+1+pAR+j] <- dddi[p+1+pAR+j, p+1+pAR+j] <- 2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)+
    2*c.^2*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[j]]%*%delta+
    2*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)+
    2*c.*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%F.lista[[j]]%*%delta-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%F.lista[[j]] + F.lista[[j]]%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)
  #dalpha.dlambda
  for(i in 1:q2) dddi[p+1+pAR+i, (p+1+pAR+q2+q3):(p+1+pAR+q2+q1)] <- dddi[(p+1+pAR+q2+q3):(p+1+pAR+q2+q1), p+1+pAR+i] <- (-2*c.)*t(F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med))%*%(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)+
    2*c.^2*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)+
    2*c.*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)

  #dlambda.dlambda
  dddi[(p+1+pAR+q2+q3):(p+1+pAR+q2+q1),(p+1+pAR+q2+q3):(p+1+pAR+q2+q1)] <- (2*c./(1+as.numeric(t(lambda)%*%lambda))^.5)*(t(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)%*%(Deltab) +
                                                                                                                           Fmat%*%(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)%*%delta)%*%t(t(z1)%*%sPsi%*%(y1-med))+
    (2*c.^2/(1+as.numeric(t(lambda)%*%lambda))^.5)*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)+
    (2*c./(1+as.numeric(t(lambda)%*%lambda)))*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%(y1-med)%*%t(delta)

  ##### Segundas Derivadas de ki
  ddki <- indpar2
  ntheta <- p+1+pAR+q1+q2
  for(i in 1:ntheta) for(j in 1:ntheta) ddki[i,j] = Iphi((ni+1)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddAi[i,j] -
    .5*Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*(ddi[j]+2*Ai*dAi[j])*dAi[i]-
    .5*IPhi((ni+2)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dddi[i,j]-
    .5*(Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dAi[j]-.5*IPhi((ni+4)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddi[j])*ddi[i]

  # Resultados
  sihat = -.5*dlogdpsi+1/ki*dki # primeira derivada
  hessianHat = indpar2
  for(i in 1:ntheta) for(j in 1:ntheta) hessianHat[i,j] = -.5*ddlogdpsi[i,j] + (1/ki)*ddki[i,j] - (1/ki^2)*dki[i]*dki[j] # Segunda derivada

  return(list(sihat = sihat, hessianHat = hessianHat))
}

derBetasARi <- function(jseq,y,x,z,time,beta1,sigmae,phiAR,D1,lambda,distr,nu){
  if (distr=="sn"|distr=="norm") c.=-sqrt(2/pi)
  if (distr=="st"|distr=="t") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ssl"|distr=="sl") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn"|distr=="cn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
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
  # theta = c(beta1,sigmae,phi,dd,lambda,nu) - para AR(p)
  indpar = c(rep("beta",p))
  lpar = length(indpar)
  indpar2 = matrix(0, nrow = (p), ncol = (p))

  ##### Primeiras derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)

  ##### Primeiras derivadas de Ai para diferente de nu
  dAi = numeric(lpar)
  ai = as.numeric((1+t(lambda)%*%sFmat%*%Lambda%*%sFmat%*%lambda)^.5)
  bi = as.numeric((1+t(lambda)%*%lambda)^.5)
  Bi = as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)
  dAi[indpar=="beta"] = -1/ai*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda

  ##### Primeiras derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)

  ##### Primeiras derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = numeric(lpar)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  ##### Segundas Derivadas de Ai
  ddAi <- indpar2
  # dbeta.dbeta = 0
  ##### Segundas Derivadas de di
  dddi <- indpar2
  # dbeta.dbeta
  dddi[1:p,1:p] <- 2*t(x1)%*%sPsi%*%x1
  ##### Segundas Derivadas de ki
  ddki <- indpar2
  ntheta <- p
  for(i in 1:ntheta) for(j in 1:ntheta) ddki[i,j] = Iphi((ni+1)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddAi[i,j] -
    .5*Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*(ddi[j]+2*Ai*dAi[j])*dAi[i]-
    .5*IPhi((ni+2)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dddi[i,j]-
    .5*(Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dAi[j]-.5*IPhi((ni+4)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddi[j])*ddi[i]

  # Resultados
  sihat = -.5*dlogdpsi+1/ki*dki # primeira derivada
  hessianHat = indpar2
  for(i in 1:ntheta) for(j in 1:ntheta) hessianHat[i,j] = (1/ki)*ddki[i,j] - (1/ki^2)*dki[i]*dki[j] # Segunda derivada

  return(list(sihat = sihat, hessianHat = hessianHat))
}

expectBetasARi <- function(jseq,y,x,z,time,beta1,sigmae,phiAR,D1,lambda,distr,nu){
  y1 = y[jseq]
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
  med<-x1%*%beta1 #+ c.*z1%*%Deltab
  MniAR<-CovARp(phi = phiAR,t1)
  sMniAR<-solve(MniAR)
  Sigma <- sigmae*MniAR
  Psi<-(z1)%*%(Gammab)%*%t(z1)+Sigma
  sPsi <- solve(Psi)
  di<-as.numeric(t(y1-med)%*%sPsi%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ai<-as.numeric(mutj/sqrt(Mtj2))
  sFmat = solve(Fmat)
  Lambda = solve(solve(D1)+ t(z1)%*%solve(Sigma)%*%z1)
  F.lista <- lapply(1:q2,F.r,q1=q1)

  if (distr=="sn"){
    e1 = e4 = e7 = 1# =matrix(1, nrow = (p), ncol = (p))
    e2 = e5 = Psi + med%*%t(med)
    e3 = e6 = med
  }
  if (distr=="st"){
    auxc <- function(nnu,r,p) (gamma((p+nnu)/2)*gamma((nnu+2*r)/2))/(nnu^r*gamma(nnu/2)*gamma((p+nnu+2*r)/2))
    e1 = (ni+nu)*auxc(nu,1,ni)
    e2 = (ni+nu+2)*(ni+nu)*auxc(nu,2,ni)*(med%*%t(med) + (nu/(nu+2))*Psi)
    e3 = (ni+nu+2)*(ni+nu)*auxc(nu,2,ni)*med
    e4 = (ni+nu+2)*(ni+nu)*auxc(nu,2,ni)
    e5 = (ni+nu)^2*auxc(nu,2,ni)*(med%*%t(med) + (nu/(nu+2))*Psi)
    e6 = (ni+nu)^2*auxc(nu,2,ni)*med
    e7 = (ni+nu)^2*auxc(nu,2,ni)
  }
  #if (distr=="ssl"){}
  #if (distr=="scn"){}

  ## E1 = E(d2log(det(Psi))) = 0
  ## E2 = E(1/Ki*d2Ki)
  E2i = -e1*t(x1)%*%sPsi%*%x1 + t(x1)%*%sPsi%*%e2%*%sPsi%*%x1 -
    2*t(x1)%*%sPsi%*%e3%*%t(med)%*%sPsi%*%x1 + e4*t(x1)%*%sPsi%*%med%*%t(med)%*%sPsi%*%x1
  ## E3 E(1/Ki^2*dKi/dbeta*dKi/dbeta)
  E3i = t(x1)%*%sPsi%*%e5%*%sPsi%*%x1 - 2*t(x1)%*%sPsi%*%e6%*%t(med)%*%sPsi%*%x1+
    e7*t(x1)%*%sPsi%*%med%*%t(med)%*%sPsi%*%x1
  ## E4 = E(d2log(det(Psi))/dbeta*d2log(det(Psi))/dbeta) = 0
  ## E5 = E(d2log(det(Psi))/dbeta*1/Ki*dKi/dbeta) = 0

  return(list(E2i = E2i, E3i = E3i))
}


# SMSN - UNC -------------------------------------------------------------------
derivatesUNC <- function(jseq,y,x,z,time,beta1,sigmae,D1,lambda,distr,nu){
  if (distr=="sn"|distr=="norm") c.=-sqrt(2/pi)
  if (distr=="st"|distr=="t") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ssl"|distr=="sl") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn"|distr=="cn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[jseq]
  t1 = time[jseq]
  p= ncol(x);q1=ncol(z)
  q2 = q1*(q1+1)/2
  if(q1 != 0){ q3 = 1} else{q3 = 0}
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
  Lambda = solve(solve(D1)+ t(z1)%*%solve(Sigma)%*%z1)
  F.lista <- lapply(1:q2,F.r,q1=q1)
  indpar = c(rep("beta",p),"sigma",rep("dd",q2),rep("lambda",q1))
  lpar = length(indpar)
  indpar2 = matrix(0, nrow = (p+1+q2+q1), ncol = (p+1+q2+q1))

  ##### Primeiras derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)
  dlogdpsi[indpar=="sigma"] =traceM(sPsi)
  for (i in 1:q2) dlogdpsi[indpar=="dd"][i] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                    Fmat%*%F.lista[[i]])%*%t(z1))
  ##### Primeiras derivadas de Ai para diferente de nu
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
                                                 t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 c./bi*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda)+
    1/ai^2*Ai/2*t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                       Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda

  ##### Primeiras derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)
  ddi[indpar=="lambda"] = -2*c./bi*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="sigma"] = -t(y1-med)%*%sPsi%*%sPsi%*%(y1-med)
  for (i in 1:q2) ddi[indpar=="dd"][i] =-2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)

  ##### Primeiras derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = numeric(lpar)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  ##### Segundas Derivadas de log(det(psi))
  ddlogdpsi = indpar2
  # dsigmae.dsigmae
  ddlogdpsi[p+1, p+1] = traceM(-sPsi%*%sPsi)
  # dsigmae.dalphar
  for (i in 1:q2) ddlogdpsi[p+1,p+1+i] = ddlogdpsi[p+1+i,p+1] = traceM(-sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                                       Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi)
  # dalpha.dalpha
  for (i in 1:q2) for (j in 1:q2) ddlogdpsi[p+1+i,p+1+j] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%F.lista[[j]] + F.lista[[j]]%*%F.lista[[i]])%*%t(z1)-
                                                                    sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%
                                                                    sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1))
  ##### Segundas Derivadas de Ai
  ddAi <- indpar2
  # dbeta.dbeta = 0

  # dbeta.dsigmae
  ddAi[1:p, p+1] <- (1/ai)*t(x1)%*%sPsi%*%sPsi%*%z1%*%Fmat%*%lambda+
    (1/(2*ai^3*sigmae^2))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda
  ddAi[p+1, 1:p] <- t(ddAi[1:p, p+1])
  #dbeta.dalphar
  for (i in 1:q2) ddAi[1:p, p+1+i] <- (1/ai)*t(x1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda-
    (1/ai)*t(x1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda - (1/(2*ai^3))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                                               Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat + sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:q2) ddAi[p+1+i,1:p] <- ddAi[1:p, p+1+i]
  # dbeta.dlambda
  ddAi[1:p, (p+1+q2+q3):(p+1+q2+q1)] <- (-1/ai)*t(x1)%*%sPsi%*%z1%*%Fmat +(1/ai^3)*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat
  ddAi[(p+1+q2+q3):(p+1+q2+q1), 1:p] <- ddAi[1:p, (p+1+q2+q3):(p+1+q2+q1)]

  # dsigma.dsigma
  ddAi[p+1,p+1] <- (2/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%sPsi%*%(y1-med)+
    (1/(2*sigmae^2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="sigma"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda -
    (1/(sigmae^4*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(sigmae^3*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda +
    (1/(2*sigmae^4*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda

  #dsigma.dalphar
  for (i in 1:q2) ddAi[p+1, p+1+i] <- (-1/ai)*t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%z1%*%F.lista[[i]]%*%delta-
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                            Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+ sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="dd"][i]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%F.lista[[i]]%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%(sFmat%*%F.lista[[i]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[i]]%*%sFmat)%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%(sFmat%*%F.lista[[i]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[i]]%*%sFmat)%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%F.lista[[i]]%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^4))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                                        Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+ sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:q2) ddAi[p+1+i, p+1] <- ddAi[p+1, p+1+i]

  #dsigmae.dlambda
  ddAi[p+1, (p+1+q2+q3):(p+1+q2+q1)] <- (-1/ai)*t(Fmat%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med))+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))+
    (1/ai^3)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat -
    (1/(2*sigmae^2*ai^2))*(t(dAi[indpar=="lambda"])*as.numeric(t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda)+
                             2*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat)+
    (1/(sigmae^2*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat

  ddAi[(p+1+q2+q3):(p+1+q2+q1), p+1] <- ddAi[p+1,(p+1+q2+q3):(p+1+q2+q1)] #ddAi[p+1+q2+q1, p+1] <- ddAi[p+1, p+1+q2+q1]

  #dalpha.dalpha
  for(r in 1:q2) for(s in 1:q2) ddAi[p+1+r, p+1+s] <- (-1/ai)*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]]+F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)-
    (c./ai)*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[s]]%*%delta+
    (1/(2*ai^3))*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                             Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (c./ai)*t(delta)%*%F.lista[[s]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda+
    (c./ai)*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda-
    (c./(2*ai^3))*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                 Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/ai)*t(lambda)%*%F.lista[[s]]%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)-
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[r]]%*%F.lista[[s]] + F.lista[[s]]%*%F.lista[[r]])%*%t(z1)%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[s]]%*%delta-
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                                                                       Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda+
    (1/(2*ai^2))*dAi[indpar=="dd"][s]%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                                               Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%F.lista[[s]]%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                                                    Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%F.lista[[s]]%*%sFmat%*%lambda+
    (1/(2*ai^4))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                                                                                      Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda+
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(-F.lista[[r]]%*%sFmat%*%F.lista[[s]]%*%sFmat%*%Lambda+
                                             F.lista[[r]]%*%sFmat%*%Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda+
                                             Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%F.lista[[s]]%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%F.lista[[s]]%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]]%*%sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%F.lista[[s]]%*%sFmat%*%Lambda-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda)%*%sFmat%*%lambda
  # dalpha.dlambda
  for(i in 1:q2) ddAi[p+1+i, (p+1+q2+q3):(p+1+q2+q1)] <- (1/ai)*t(F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med))-
    (c./ai)*t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))-
    (1/ai^3)*t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat-
    (c./ai)*t(t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda)-
    (c./ai)*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]+ #talvez t
    (c./ai^3)*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat-
    (1/ai)*t(Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med))+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))+
    (1/ai^3)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat+
    (1/(2*ai^2))*t(dAi[indpar=="lambda"])*as.numeric(t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                            Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda)+
    (1/(2*ai^2))*Ai%*%t(sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                   Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda)+
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                             Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat-
    (1/(ai^4))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                           Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat

  for(i in 1:q2) ddAi[(p+1+q2+q3):(p+1+q2+q1),p+1+i] <- ddAi[p+1+i, (p+1+q2+q3):(p+1+q2+q1)]

  # dlambda.dlambda
  ddAi[(p+1+q2+q3):(p+1+q2+q1), (p+1+q2+q3):(p+1+q2+q1)] <- (-2*c./ai)*Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))-
    (1/ai^3)*Fmat%*%t(z1)%*%sPsi%*%(y1-x1%*%beta1 - 2*c.*z1%*%Deltab)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat+
    as.numeric(2*c./(ai*(1+t(lambda)%*%lambda)^1.5))*Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)+
    as.numeric(c./(ai*(1+t(lambda)%*%lambda)))*as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)*t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))-
    as.numeric(c./(ai^2*(1+t(lambda)%*%lambda)^2.5))*as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)*lambda%*%((1/ai)*t(lambda)%*%sFmat%*%Lambda%*%sFmat*as.numeric((1+t(lambda)%*%lambda))+2*ai*t(lambda))-
    (1/ai^2)*dAi[indpar=="lambda"]%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat - (1/ai^2)*Ai*sFmat%*%Lambda%*%sFmat + (2/ai^4)*Ai*sFmat%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat

  ##### Segundas Derivadas de di
  dddi <- indpar2
  # dbeta.dbeta
  dddi[1:p,1:p] <- 2*t(x1)%*%sPsi%*%x1
  # dbeta.dsigmae
  dddi[1:p, p+1] <- dddi[p+1, 1:p] <- 2*t(x1)%*%sPsi%*%sPsi%*%(y1-med)
  # dbeta.dalphar
  for(i in 1:q2) dddi[1:p, p+1+i] <- dddi[p+1+i,1:p] <- 2*t(x1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)+ 2*c.*t(x1)%*%sPsi%*%z1%*%F.lista[[i]]%*%delta
  # dbeta.dlambda
  dddi[1:p, (p+1+q2+q3):(p+1+q2+q1)] <- dddi[(p+1+q2+q3):(p+1+q2+q1),1:p] <- 2*c.*t(x1)%*%sPsi%*%z1%*%Fmat%*%t(1/((1+ as.numeric(t(lambda)%*%lambda))^.5)*diag(q1) - lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)
  # dsigma.dsigma
  dddi[p+1,p+1] <- 2*t(y1-med)%*%sPsi%*%sPsi%*%sPsi%*%(y1-med)
  # dsigmae.dalpha
  for(i in 1:q2) dddi[p+1,p+1+i] <- dddi[p+1+i, p+1] <- c.*t(z1%*%F.lista[[i]]%*%delta)%*%sPsi%*%sPsi%*%(y1-med) + t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)+
    t(y1-med)%*%sPsi%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med) + c.*t(y1-med)%*%sPsi%*%sPsi%*%z1%*%F.lista[[i]]%*%delta
  # dsigmae.dlambda
  dddi[p+1,(p+1+q2+q3):(p+1+q2+q1)] <- dddi[(p+1+q2+q3):(p+1+q2+q1), p+1]<- (2*c./(1+as.numeric(t(lambda)%*%lambda))^.5)*(Fmat - delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)
  # dalpha.dalpha
  for(i in 1:q2) for(j in 1:q2) dddi[p+1+i, p+1+j] <- dddi[p+1+j, p+1+j] <- 2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)+
    2*c.^2*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[j]]%*%delta+
    2*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)+
    2*c.*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%F.lista[[j]]%*%delta-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%F.lista[[j]] + F.lista[[j]]%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)
  #dalpha.dlambda
  for(i in 1:q2) dddi[p+1+i, (p+1+q2+q3):(p+1+q2+q1)] <- dddi[(p+1+q2+q3):(p+1+q2+q1), p+1+i] <- (-2*c.)*t(F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med))%*%(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)+
    2*c.^2*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)+
    2*c.*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)
  #dlambda.dlambda
  dddi[(p+1+q2+q3):(p+1+q2+q1),(p+1+q2+q3):(p+1+q2+q1)] <- (2*c./(1+as.numeric(t(lambda)%*%lambda))^.5)*(t(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)%*%(Deltab) +
                                                                                                           Fmat%*%(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)%*%delta)%*%t(t(z1)%*%sPsi%*%(y1-med))+
    (2*c.^2/(1+as.numeric(t(lambda)%*%lambda))^.5)*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)+
    (2*c./(1+as.numeric(t(lambda)%*%lambda)))*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%(y1-med)%*%t(delta)

  ##### Segundas Derivadas de ki
  ddki <- indpar2
  ntheta <- p+1+q1+q2
  for(i in 1:ntheta) for(j in 1:ntheta) ddki[i,j] = Iphi((ni+1)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddAi[i,j] -
    .5*Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*(ddi[j]+2*Ai*dAi[j])*dAi[i]-
    .5*IPhi((ni+2)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dddi[i,j]-
    .5*(Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dAi[j]-.5*IPhi((ni+4)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddi[j])*ddi[i]

  # Resultados
  sihat = -.5*dlogdpsi+1/ki*dki # primeira derivada
  hessianHat = indpar2
  for(i in 1:ntheta) for(j in 1:ntheta) hessianHat[i,j] = -.5*ddlogdpsi[i,j] + (1/ki)*ddki[i,j] - (1/ki^2)*dki[i]*dki[j] # Segunda derivada

  return(list(sihat = sihat, hessianHat = hessianHat))
}

derBetasUNC <- function(jseq,y,x,z,time,beta1,sigmae,D1,lambda,distr,nu){
  if (distr=="sn"|distr=="norm") c.=-sqrt(2/pi)
  if (distr=="st"|distr=="t") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ssl"|distr=="sl") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn"|distr=="cn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
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
  Sigma <- sigmae*diag(ni)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  sPsi <- solve(Psi)
  di<-as.numeric(t(y1-med)%*%sPsi%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ai<-as.numeric(mutj/sqrt(Mtj2))
  sFmat = solve(Fmat)
  Lambda = solve(solve(D1)+ t(z1)%*%solve(Sigma)%*%z1)
  F.lista <- lapply(1:q2,F.r,q1=q1)
  # theta = c(beta1,sigmae,phi,dd,lambda,nu) - para AR(p)
  indpar = c(rep("beta",p))
  lpar = length(indpar)
  indpar2 = matrix(0, nrow = (p), ncol = (p))

  ##### Primeiras derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)

  ##### Primeiras derivadas de Ai para diferente de nu
  dAi = numeric(lpar)
  ai = as.numeric((1+t(lambda)%*%sFmat%*%Lambda%*%sFmat%*%lambda)^.5)
  bi = as.numeric((1+t(lambda)%*%lambda)^.5)
  Bi = as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)
  dAi[indpar=="beta"] = -1/ai*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda

  ##### Primeiras derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)

  ##### Primeiras derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = numeric(lpar)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  ##### Segundas Derivadas de Ai
  ddAi <- indpar2
  # dbeta.dbeta = 0
  ##### Segundas Derivadas de di
  dddi <- indpar2
  # dbeta.dbeta
  dddi[1:p,1:p] <- 2*t(x1)%*%sPsi%*%x1
  ##### Segundas Derivadas de ki
  ddki <- indpar2
  ntheta <- p
  for(i in 1:ntheta) for(j in 1:ntheta) ddki[i,j] = Iphi((ni+1)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddAi[i,j] -
    .5*Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*(ddi[j]+2*Ai*dAi[j])*dAi[i]-
    .5*IPhi((ni+2)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dddi[i,j]-
    .5*(Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dAi[j]-.5*IPhi((ni+4)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddi[j])*ddi[i]

  # Resultados
  sihat = -.5*dlogdpsi+1/ki*dki # primeira derivada
  hessianHat = indpar2
  for(i in 1:ntheta) for(j in 1:ntheta) hessianHat[i,j] = (1/ki)*ddki[i,j] - (1/ki^2)*dki[i]*dki[j] # Segunda derivada

  return(list(sihat = sihat, hessianHat = hessianHat))
}

expectBetasUNC <- function(jseq,y,x,z,time,beta1,sigmae,D1,lambda,distr,nu){
  y1 = y[jseq]
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
  med<-x1%*%beta1 #+ c.*z1%*%Deltab
  Sigma <- sigmae*diag(ni)
  Psi<-(z1)%*%(Gammab)%*%t(z1)+Sigma
  sPsi <- solve(Psi)
  di<-as.numeric(t(y1-med)%*%sPsi%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ai<-as.numeric(mutj/sqrt(Mtj2))
  sFmat = solve(Fmat)
  Lambda = solve(solve(D1)+ t(z1)%*%solve(Sigma)%*%z1)
  F.lista <- lapply(1:q2,F.r,q1=q1)

  if (distr=="sn"){
    e1 = e4 = e7 = 1# =matrix(1, nrow = (p), ncol = (p))
    e2 = e5 = Psi + med%*%t(med)
    e3 = e6 = med
  }
  if (distr=="st"){
    auxc <- function(nnu,r,p) (gamma((p+nnu)/2)*gamma((nnu+2*r)/2))/(nnu^r*gamma(nnu/2)*gamma((p+nnu+2*r)/2))
    e1 = (ni+nu)*auxc(nu,1,ni)
    e2 = (ni+nu+2)*(ni+nu)*auxc(nu,2,ni)*(med%*%t(med) + (nu/(nu+2))*Psi)
    e3 = (ni+nu+2)*(ni+nu)*auxc(nu,2,ni)*med
    e4 = (ni+nu+2)*(ni+nu)*auxc(nu,2,ni)
    e5 = (ni+nu)^2*auxc(nu,2,ni)*(med%*%t(med) + (nu/(nu+2))*Psi)
    e6 = (ni+nu)^2*auxc(nu,2,ni)*med
    e7 = (ni+nu)^2*auxc(nu,2,ni)
  }
  if (distr=="ssl"){}
  if (distr=="scn"){}

  ## E1 = E(d2log(det(Psi))) = 0
  ## E2 = E(1/Ki*d2Ki)
  E2i = -e1*t(x1)%*%sPsi%*%x1 + t(x1)%*%sPsi%*%e2%*%sPsi%*%x1 -
    2*t(x1)%*%sPsi%*%e3%*%t(med)%*%sPsi%*%x1 + e4*t(x1)%*%sPsi%*%med%*%t(med)%*%sPsi%*%x1
  ## E3 E(1/Ki^2*dKi/dbeta*dKi/dbeta)
  E3i = t(x1)%*%sPsi%*%e5%*%sPsi%*%x1 - 2*t(x1)%*%sPsi%*%e6%*%t(med)%*%sPsi%*%x1+
    e7*t(x1)%*%sPsi%*%med%*%t(med)%*%sPsi%*%x1
  ## E4 = E(d2log(det(Psi))/dbeta*d2log(det(Psi))/dbeta) = 0
  ## E5 = E(d2log(det(Psi))/dbeta*1/Ki*dKi/dbeta) = 0

  return(list(E2i = E2i, E3i = E3i))
}

# SMSN - DEC -------------------------------------------------------------------
derivatesDEC <- function(jseq,y,x,z,time,beta1,sigmae,phiDEC,thetaDEC,D1,lambda,distr,nu){
  if (distr=="sn"|distr=="norm") c.=-sqrt(2/pi)
  if (distr=="st"|distr=="t") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ssl"|distr=="sl") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn"|distr=="cn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[jseq]
  t1 = time[jseq]
  p= ncol(x);q1=ncol(z); #pAR= 2 #length(phiAR)
  q2 = q1*(q1+1)/2
  if(q1 != 0){ q3 = 1} else{q3 = 0}
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  ni = length(y1)
  Fmat = matrix.sqrt(D1)
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-Fmat%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  Covmat <- CovDEC(phiDEC,thetaDEC,t1)  #Covmat <- CovDEC(phiDEC, thetaDEC, t1)
  sCovmat<-solve(Covmat)   #sCovmat<-solve(Covmat)
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
  # indpar = c(rep("beta",p),"sigma",rep("phi",pAR),rep("dd",q2),rep("lambda",q1))
  indpar = c(rep("beta",p),"sigma","phi","theta",rep("dd",q2),rep("lambda",q1))
  lpar = length(indpar)
  indpar2 = matrix(0, nrow = (p+1+2+q2+q1), ncol = (p+1+2+q2+q1))

  ##### Primeiras derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)
  dlogdpsi[indpar=="sigma"] =traceM(sPsi%*%Covmat)
  for (i in 1:q2) dlogdpsi[indpar=="dd"][i] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                    Fmat%*%F.lista[[i]])%*%t(z1))
  dphiDEC <- dphiCovDEC(phiDEC, thetaDEC, t1)
  dthetaDEC <- dthetaCovDEC(phiDEC, thetaDEC, t1)
  dphiphiDEC <- dphiphiCovDEC(phiDEC, thetaDEC, t1)
  dthetathetaDEC <- dthetathetaCovDEC(phiDEC, thetaDEC, t1)
  dphithetaDEC <- dphithetaCovDEC(phiDEC, thetaDEC, t1)

  ##### Primeiras derivadas de Ai para diferente de nu
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

  ##### Primeiras derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)
  ddi[indpar=="lambda"] = -2*c./bi*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="sigma"] = -t(y1-med)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)
  for (i in 1:q2) ddi[indpar=="dd"][i] =-2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="phi"] = -sigmae*t(y1-med)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)
  ddi[indpar=="theta"] = -sigmae*t(y1-med)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)

  ##### Primeiras derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = numeric(lpar)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  ##### Segundas Derivadas de log(det(psi))
  ddlogdpsi = indpar2
  # dsigmae.dsigmae
  ddlogdpsi[p+1, p+1] = traceM(-sPsi%*%Covmat%*%sPsi%*%Covmat)
  # dsigmae.dphi
  ddlogdpsi[p+1,p+1+1]= ddlogdpsi[p+1+1,p+1] = traceM(sPsi%*%dphiDEC - sigmae*sPsi%*%dphiDEC%*%sPsi%*%Covmat)
  ddlogdpsi[p+2,p+2+2]= ddlogdpsi[p+2+2,p+2] = traceM(sPsi%*%dthetaDEC - sigmae*sPsi%*%dthetaDEC%*%sPsi%*%Covmat)

  # dsigmae.dalphar
  for (i in 1:q2) ddlogdpsi[p+1,p+1+2+i] = ddlogdpsi[p+1+2+i,p+1] = traceM(-sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                                           Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%Covmat)
  # dalpha.dalpha
  for (i in 1:q2) for (j in 1:q2) ddlogdpsi[p+1+2+i,p+1+2+j] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%F.lista[[j]] + F.lista[[j]]%*%F.lista[[i]])%*%t(z1)-
                                                                        sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%
                                                                        sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1))
  # dalpha.dphi
  for (i in 1:q2) ddlogdpsi[p+1+1+i, p+1+1] = ddlogdpsi[p+1+1, p+1+1+i] = -traceM(sigmae*sPsi%*%dphiDEC%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+ Fmat%*%F.lista[[i]])%*%t(z1))
  for (i in 1:q2) ddlogdpsi[p+1+2+i, p+1+2] = ddlogdpsi[p+1+2, p+1+2+i] = -traceM(sigmae*sPsi%*%dthetaDEC%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1))

  # dphi.dphi
  ddlogdpsi[p+1+1,p+1+1] = sigmae*traceM(sPsi%*%dphiphiDEC - sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC)
  ddlogdpsi[p+1+2,p+1+2] = sigmae*traceM(sPsi%*%dthetathetaDEC - sigmae*sPsi%*%dthetaDEC%*%sPsi%*%dthetaDEC)
  ddlogdpsi[p+1+1,p+1+2] = ddlogdpsi[p+1+2,p+1+1] = sigmae*traceM(sPsi%*%dphithetaDEC - sigmae*sPsi%*%dphiDEC%*%sPsi%*%dthetaDEC)

  ##### Segundas Derivadas de Ai
  ddAi <- indpar2
  # dbeta.dbeta = 0

  # dbeta.dsigmae
  ddAi[1:p, p+1] <- (1/ai)*t(x1)%*%sPsi%*%Covmat%*%sPsi%*%z1%*%Fmat%*%lambda+
    (1/(2*ai^3*sigmae^2))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda
  ddAi[p+1, 1:p] <- t(ddAi[1:p, p+1])
  # dbeta.dphi
  ddAi[1:p, p+1+1] <- (sigmae/ai)*t(x1)%*%sPsi%*%dphiDEC%*%sPsi%*%z1%*%Fmat%*%lambda +
    (1/(2*ai^3*sigmae))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda
  ddAi[1:p, p+1+2] <- (sigmae/ai)*t(x1)%*%sPsi%*%dthetaDEC%*%sPsi%*%z1%*%Fmat%*%lambda +
    (1/(2*ai^3*sigmae))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  ddAi[p+1+1, 1:p] <- ddAi[1:p, p+1+1]
  ddAi[p+1+2, 1:p] <- ddAi[1:p, p+1+2]

  #dbeta.dalphar
  for (i in 1:q2) ddAi[1:p, p+1+2+i] <- (1/ai)*t(x1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda-
    (1/ai)*t(x1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda - (1/(2*ai^3))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                                               Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat + sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:q2) ddAi[p+1+2+i,1:p] <- ddAi[1:p, p+1+2+i]
  # dbeta.dlambda
  ddAi[1:p, (p+1+2+q2+q3):(p+1+2+q2+q1)] <- (-1/ai)*t(x1)%*%sPsi%*%z1%*%Fmat +(1/ai^3)*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat
  ddAi[(p+1+2+q2+q3):(p+1+2+q2+q1), 1:p] <- ddAi[1:p, (p+1+2+q2+q3):(p+1+2+q2+q1)]

  # dsigma.dsigma
  ddAi[p+1,p+1] <- (2/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)+
    (1/(2*sigmae^2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="sigma"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda -
    (1/(sigmae^4*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(sigmae^3*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda +
    (1/(2*sigmae^4*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  #dsigma.dphi
  ddAi[p+1,p+1+1] <- (-1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)+
    (2*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)+
    (1/(2*ai^3*sigmae))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="phi"][i]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(sigmae^3*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^3*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda
  ddAi[p+1+1,p+1] <- ddAi[p+1,p+1+1]

  ddAi[p+1,p+1+2] <- (-1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)+
    (2*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)+
    (1/(2*ai^3*sigmae))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="theta"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(sigmae^3*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^3*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda
  ddAi[p+1+2,p+1] <- ddAi[p+1,p+1+2]

  #dsigma.dalphar
  for (i in 1:q2) ddAi[p+1, p+1+2+i] <- (-1/ai)*t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%z1%*%F.lista[[i]]%*%delta-
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                     Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+ sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="dd"][i]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%F.lista[[i]]%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%(sFmat%*%F.lista[[i]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[i]]%*%sFmat)%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%(sFmat%*%F.lista[[i]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[i]]%*%sFmat)%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%F.lista[[i]]%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^4))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                                                  Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+ sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:q2) ddAi[p+1+2+i, p+1] <- ddAi[p+1, p+1+2+i]

  #dsigmae.dlambda
  ddAi[p+1, (p+1+2+q2+q3):(p+1+2+q2+q1)] <- (-1/ai)*t(Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med))+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))+
    (1/ai^3)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat -
    (1/(2*sigmae^2*ai^2))*(t(dAi[indpar=="lambda"])*as.numeric(t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda)+
                             2*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat)+
    (1/(sigmae^2*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat

  ddAi[(p+1+2+q2+q3):(p+1+2+q2+q1), p+1] <- ddAi[p+1,(p+1+2+q2+q3):(p+1+2+q2+q1)] #ddAi[p+1+pAR+q2+q1, p+1] <- ddAi[p+1, p+1+pAR+q2+q1]

  #dphi.dphi
  ddAi[p+1+1, p+1+1] <- (-sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%
    (-sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC%*%sPsi + sPsi%*%dphiphiDEC%*%sPsi -
       sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC%*%sPsi)%*%(y1-med)+ (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*dAi[indpar=="phi"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%
    (-sCovmat%*%dphiDEC%*%sCovmat%*%dphiDEC%*%sCovmat+ sCovmat%*%dphiphiDEC%*%sCovmat - sCovmat%*%dphiDEC%*%sCovmat%*%dphiDEC%*%sCovmat)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*ai^4*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  ddAi[p+1+2, p+1+2] <- (-sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%
    (-sigmae*sPsi%*%dthetaDEC%*%sPsi%*%dthetaDEC%*%sPsi + sPsi%*%dthetathetaDEC%*%sPsi -
       sigmae*sPsi%*%dthetaDEC%*%sPsi%*%dthetaDEC%*%sPsi)%*%(y1-med)+ (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*dAi[indpar=="theta"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%
    (-sCovmat%*%dthetaDEC%*%sCovmat%*%dthetaDEC%*%sCovmat+ sCovmat%*%dthetathetaDEC%*%sCovmat - sCovmat%*%dthetaDEC%*%sCovmat%*%dthetaDEC%*%sCovmat)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*ai^4*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  ddAi[p+1+1, p+1+2] <- (-sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%
    (-sigmae*sPsi%*%dthetaDEC%*%sPsi%*%dphiDEC%*%sPsi + sPsi%*%dphithetaDEC%*%sPsi - sigmae*sPsi%*%dphiDEC%*%sPsi%*%dthetaDEC%*%sPsi)%*%(y1-med)+
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*dAi[indpar=="theta"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%
    (-sCovmat%*%dthetaDEC%*%sCovmat%*%dphiDEC%*%sCovmat+sCovmat%*%dphithetaDEC%*%sCovmat - sCovmat%*%dphiDEC%*%sCovmat%*%dthetaDEC%*%sCovmat)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*ai^4*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  ddAi[p+1+2, p+1+1] = ddAi[p+1+1, p+1+2]
  #dphi.dalphar
  for(j in 1:q2) ddAi[p+1+1, p+1+2+j] <- (-sigmae/ai)*t(lambda)%*%F.lista[[j]]%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)+
    (sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[j]] + F.lista[[j]]%*%Fmat)%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)+
    (sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%z1%*%(Fmat%*%F.lista[[j]] + F.lista[[j]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c.*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%z1%*%F.lista[[j]]%*%delta-
    (sigmae/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[j]]%*%sFmat%*%Lambda+
                                                                                                           Lambda%*%sFmat%*%F.lista[[j]]-
                                                                                                           Lambda%*%sFmat%*%(F.lista[[j]]%*%sFmat+
                                                                                                                               sFmat%*%F.lista[[j]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*dAi[indpar=="dd"][j]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%F.lista[[j]]%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%(sFmat%*%F.lista[[j]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[j]]%*%sFmat)%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%(sFmat%*%F.lista[[j]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[j]]%*%sFmat)%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%F.lista[[j]]%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%(F.lista[[j]]%*%sFmat%*%Lambda+
                                                                                                                                              Lambda%*%sFmat%*%F.lista[[j]]-
                                                                                                                                              Lambda%*%sFmat%*%(F.lista[[j]]%*%sFmat+
                                                                                                                                                                  sFmat%*%F.lista[[j]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for(j in 1:q2) ddAi[p+1+2, p+1+2+j] <- (-sigmae/ai)*t(lambda)%*%F.lista[[j]]%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)+
    (sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[j]] + F.lista[[j]]%*%Fmat)%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)+
    (sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%z1%*%(Fmat%*%F.lista[[j]] + F.lista[[j]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c.*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%z1%*%F.lista[[j]]%*%delta-
    (sigmae/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[j]]%*%sFmat%*%Lambda+
                                                                                                             Lambda%*%sFmat%*%F.lista[[j]]-
                                                                                                             Lambda%*%sFmat%*%(F.lista[[j]]%*%sFmat+
                                                                                                                                 sFmat%*%F.lista[[j]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*dAi[indpar=="dd"][j]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%F.lista[[j]]%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%(sFmat%*%F.lista[[j]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[j]]%*%sFmat)%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%(sFmat%*%F.lista[[j]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[j]]%*%sFmat)%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%F.lista[[j]]%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%(F.lista[[j]]%*%sFmat%*%Lambda+
                                                                                                                                                Lambda%*%sFmat%*%F.lista[[j]]-
                                                                                                                                                Lambda%*%sFmat%*%(F.lista[[j]]%*%sFmat+
                                                                                                                                                                    sFmat%*%F.lista[[j]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for(j in 1:q2) ddAi[p+1+2+ j, p+1+1] <- ddAi[p+1+i, p+1+2+1]
  for(j in 1:q2) ddAi[p+1+2+ j, p+1+2] <- ddAi[p+1+i, p+1+2+2]

  #dphi.dlambda
  ddAi[p+1+1, (p+1+2+q2+q3):(p+1+2+q2+q1)] <- (-sigmae/ai)*t(Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med))+
    (c.*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))+
    (sigmae/ai^3)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat-
    (1/(2*sigmae*ai^2))*t(dAi[indpar=="lambda"])*as.numeric(t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda)-
    (1/(sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat+
    (1/(sigmae*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat

  ddAi[p+1+2, (p+1+2+q2+q3):(p+1+2+q2+q1)] <- (-sigmae/ai)*t(Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med))+
    (c.*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))+
    (sigmae/ai^3)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat-
    (1/(2*sigmae*ai^2))*t(dAi[indpar=="lambda"])*as.numeric(t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda)-
    (1/(sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat+
    (1/(sigmae*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat


  ddAi[(p+1+2+q2+q3):(p+1+2+q2+q1),p+1+1] <- ddAi[p+1+1, (p+1+2+q2+q3):(p+1+1+q2+q1)]
  ddAi[(p+1+2+q2+q3):(p+1+2+q2+q1),p+1+2] <- ddAi[p+1+2, (p+1+2+q2+q3):(p+1+2+q2+q1)]

  #dalpha.dalpha
  for(r in 1:q2) for(s in 1:q2) ddAi[p+1+2+r, p+1+2+s] <- (-1/ai)*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]]+F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)-
    (c./ai)*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[s]]%*%delta+
    (1/(2*ai^3))*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                             Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (c./ai)*t(delta)%*%F.lista[[s]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda+
    (c./ai)*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda-
    (c./(2*ai^3))*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                 Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/ai)*t(lambda)%*%F.lista[[s]]%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)-
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[r]]%*%F.lista[[s]] + F.lista[[s]]%*%F.lista[[r]])%*%t(z1)%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[s]]%*%delta-
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                                                                       Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda+
    (1/(2*ai^2))*dAi[indpar=="dd"][s]%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                                               Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%F.lista[[s]]%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                                                    Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%F.lista[[s]]%*%sFmat%*%lambda+
    (1/(2*ai^4))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                                                                                      Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda+
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(-F.lista[[r]]%*%sFmat%*%F.lista[[s]]%*%sFmat%*%Lambda+
                                             F.lista[[r]]%*%sFmat%*%Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda+
                                             Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%F.lista[[s]]%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%F.lista[[s]]%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]]%*%sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%F.lista[[s]]%*%sFmat%*%Lambda-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda)%*%sFmat%*%lambda
  # dalpha.dlambda
  for(i in 1:q2) ddAi[p+1+2+i, (p+1+2+q2+q3):(p+1+2+q2+q1)] <- (1/ai)*t(F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med))-
    (c./ai)*t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))-
    (1/ai^3)*t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat-
    (c./ai)*t(t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda)-
    (c./ai)*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]+ #talvez t
    (c./ai^3)*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat-
    (1/ai)*t(Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med))+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))+
    (1/ai^3)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat+
    (1/(2*ai^2))*t(dAi[indpar=="lambda"])*as.numeric(t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                            Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda)+
    (1/(2*ai^2))*Ai%*%t(sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                   Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda)+
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                             Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat-
    (1/(ai^4))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                           Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat

  for(i in 1:q2) ddAi[(p+1+2+q2+q3):(p+1+2+q2+q1),p+1+2+i] <- ddAi[p+1+2+i, (p+1+2+q2+q3):(p+1+2+q2+q1)]

  # dlambda.dlambda
  ddAi[(p+1+2+q2+q3):(p+1+2+q2+q1), (p+1+2+q2+q3):(p+1+2+q2+q1)] <- (-2*c./ai)*Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))-
    (1/ai^3)*Fmat%*%t(z1)%*%sPsi%*%(y1-x1%*%beta1 - 2*c.*z1%*%Deltab)%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat+
    as.numeric(2*c./(ai*(1+t(lambda)%*%lambda)^1.5))*Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)+
    as.numeric(c./(ai*(1+t(lambda)%*%lambda)))*as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)*t(as.numeric(1/((1+t(lambda)%*%lambda)^.5))*diag(q1) - (lambda%*%t(lambda))/as.numeric((1+t(lambda)%*%lambda)^1.5))-
    as.numeric(c./(ai^2*(1+t(lambda)%*%lambda)^2.5))*as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)*lambda%*%((1/ai)*t(lambda)%*%sFmat%*%Lambda%*%sFmat*as.numeric((1+t(lambda)%*%lambda))+2*ai*t(lambda))-
    (1/ai^2)*dAi[indpar=="lambda"]%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat - (1/ai^2)*Ai*sFmat%*%Lambda%*%sFmat + (2/ai^4)*Ai*sFmat%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%sFmat

  ##### Segundas Derivadas de di
  dddi <- indpar2
  # dbeta.dbeta
  dddi[1:p,1:p] <- 2*t(x1)%*%sPsi%*%x1
  # dbeta.dsigmae
  dddi[1:p, p+1] <- dddi[p+1, 1:p] <- 2*t(x1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)
  # dbeta.dphi
  dddi[1:p, p+1+1]  <- dddi[p+1+1, 1:p] <- 2*sigmae*t(x1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)
  dddi[1:p, p+1+2]  <- dddi[p+1+2, 1:p] <- 2*sigmae*t(x1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)
  # dbeta.dalphar
  for(i in 1:q2) dddi[1:p, p+1+2+i] <- dddi[p+1+2+i,1:p] <- 2*t(x1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)+ 2*c.*t(x1)%*%sPsi%*%z1%*%F.lista[[i]]%*%delta
  # dbeta.dlambda
  dddi[1:p, (p+1+2+q2+q3):(p+1+2+q2+q1)] <- dddi[(p+1+2+q2+q3):(p+1+2+q2+q1),1:p] <- 2*c.*t(x1)%*%sPsi%*%z1%*%Fmat%*%t(1/((1+ as.numeric(t(lambda)%*%lambda))^.5)*diag(q1) - lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)
  # dsigma.dsigma
  dddi[p+1,p+1] <- 2*t(y1-med)%*%sPsi%*%Covmat%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)
  # dsigma.dphi
  dddi[p+1,p+1+1] <- dddi[p+1+1,p+1] <-t(y1-med)%*%(sigmae*sPsi%*%dphiDEC%*%sPsi%*%Covmat%*%sPsi - sPsi%*%dphiDEC%*%sPsi + sigmae*sPsi%*%Covmat%*%sPsi%*%dphiDEC%*%sPsi)%*%(y1-med)
  dddi[p+1,p+1+2] <- dddi[p+1+2,p+1] <-t(y1-med)%*%(sigmae*sPsi%*%dthetaDEC%*%sPsi%*%Covmat%*%sPsi - sPsi%*%dthetaDEC%*%sPsi + sigmae*sPsi%*%Covmat%*%sPsi%*%dthetaDEC%*%sPsi)%*%(y1-med)
  # dsigmae.dalpha
  for(i in 1:q2) dddi[p+1,p+1+2+i] <- dddi[p+1+2+i, p+1] <- c.*t(z1%*%F.lista[[i]]%*%delta)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med) + t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)+
    t(y1-med)%*%sPsi%*%Covmat%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med) + c.*t(y1-med)%*%sPsi%*%Covmat%*%sPsi%*%z1%*%F.lista[[i]]%*%delta
  # dsigmae.dlambda
  dddi[p+1,(p+1+2+q2+q3):(p+1+2+q2+q1)] <- dddi[(p+1+2+q2+q3):(p+1+2+q2+q1), p+1]<- (2*c./(1+as.numeric(t(lambda)%*%lambda))^.5)*(Fmat - delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)
  # dphi.dphi
  dddi[p+1+1,p+1+1] <- sigmae*t(y1-med)%*%(sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC%*%sPsi - sPsi%*%dphiphiDEC%*%sPsi + sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC%*%sPsi)%*%(y1-med)
  dddi[p+1+2,p+1+2] <- sigmae*t(y1-med)%*%(sigmae*sPsi%*%dthetaDEC%*%sPsi%*%dthetaDEC%*%sPsi - sPsi%*%dthetathetaDEC%*%sPsi + sigmae*sPsi%*%dthetaDEC%*%sPsi%*%dthetaDEC%*%sPsi)%*%(y1-med)
  dddi[p+1+1,p+1+2] <- dddi[p+1+2,p+1+1] <- sigmae*t(y1-med)%*%(sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC%*%sPsi - sPsi%*%dphithetaDEC%*%sPsi + sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC%*%sPsi)%*%(y1-med)

  # dphi.dalpha
  for(j in 1:q2) dddi[p+1+1, p+1+2+j] <- dddi[p+1+2+j, p+1+1] <- 2*c.*sigmae*t(y1-med)%*%sPsi%*%dphiDEC%*%sPsi%*%z1%*%F.lista[[j]]%*%delta+
    2*sigmae%*%t(y1-med)%*%sPsi%*%dphiDEC%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat + Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)
  for(j in 1:q2) dddi[p+1+2, p+1+2+j] <- dddi[p+1+2+j, p+1+2] <- 2*c.*sigmae*t(y1-med)%*%sPsi%*%dthetaDEC%*%sPsi%*%z1%*%F.lista[[j]]%*%delta+
    2*sigmae%*%t(y1-med)%*%sPsi%*%dthetaDEC%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat + Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)
  # dphi. dlambda
  dddi[p+1+1, (p+1+2+q2+q3):(p+1+2+q2+q1)] <- dddi[(p+1+2+q2+q3):(p+1+2+q2+q1), p+1+1] <- 2*c.*sigmae*t(y1-med)%*%sPsi%*%dphiDEC%*%sPsi%*%z1%*%Fmat%*%t(1/(1+as.numeric(t(lambda)%*%lambda))^.5 - lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)
  dddi[p+1+2, (p+1+2+q2+q3):(p+1+2+q2+q1)] <- dddi[(p+1+2+q2+q3):(p+1+2+q2+q1), p+1+2] <- 2*c.*sigmae*t(y1-med)%*%sPsi%*%dthetaDEC%*%sPsi%*%z1%*%Fmat%*%t(1/(1+as.numeric(t(lambda)%*%lambda))^.5 - lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)

  # dalpha.dalpha
  for(i in 1:q2) for(j in 1:q2) dddi[p+1+2+i, p+1+2+j] <- dddi[p+1+2+j, p+1+2+j] <- 2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)+
    2*c.^2*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[j]]%*%delta+
    2*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)+
    2*c.*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%F.lista[[j]]%*%delta-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%F.lista[[j]] + F.lista[[j]]%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)

  #dalpha.dlambda
  for(i in 1:q2) dddi[p+1+2+i, (p+1+2+q2+q3):(p+1+2+q2+q1)] <- dddi[(p+1+2+q2+q3):(p+1+2+q2+q1), p+1+2+i] <- (-2*c.)*t(F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med))%*%(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)+
    2*c.^2*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)+
    2*c.*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)

  #dlambda.dlambda
  dddi[(p+1+2+q2+q3):(p+1+2+q2+q1),(p+1+2+q2+q3):(p+1+2+q2+q1)] <- (2*c./(1+as.numeric(t(lambda)%*%lambda))^.5)*(t(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)%*%(Deltab) +
                                                                                                                   Fmat%*%(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)%*%delta)%*%t(t(z1)%*%sPsi%*%(y1-med))+
    (2*c.^2/(1+as.numeric(t(lambda)%*%lambda))^.5)*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%t(1/(1+as.numeric(t(lambda)%*%lambda))^.5*diag(q1)-lambda%*%t(lambda)/(1+as.numeric(t(lambda)%*%lambda))^1.5)+
    (2*c./(1+as.numeric(t(lambda)%*%lambda)))*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%(y1-med)%*%t(delta)

  ##### Segundas Derivadas de ki
  ddki <- indpar2
  ntheta <- p+1+2+q1+q2
  for(i in 1:ntheta) for(j in 1:ntheta) ddki[i,j] = Iphi((ni+1)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddAi[i,j] -
    .5*Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*(ddi[j]+2*Ai*dAi[j])*dAi[i]-
    .5*IPhi((ni+2)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dddi[i,j]-
    .5*(Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dAi[j]-.5*IPhi((ni+4)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddi[j])*ddi[i]

  # Resultados
  sihat = -.5*dlogdpsi+1/ki*dki # primeira derivada
  hessianHat = indpar2
  for(i in 1:ntheta) for(j in 1:ntheta) hessianHat[i,j] = -.5*ddlogdpsi[i,j] + (1/ki)*ddki[i,j] - (1/ki^2)*dki[i]*dki[j] # Segunda derivada

  return(list(sihat = sihat, hessianHat = hessianHat))
}

derBetasDEC <- function(jseq,y,x,z,time,beta1,sigmae,phiDEC,thetaDEC, D1,lambda,distr,nu){
  if (distr=="sn"|distr=="norm") c.=-sqrt(2/pi)
  if (distr=="st"|distr=="t") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ssl"|distr=="sl") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn"|distr=="cn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
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
  Covmat <- CovDEC(phiDEC, thetaDEC, t1)
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
  # theta = c(beta1,sigmae,phi,dd,lambda,nu) - para AR(p)
  indpar = c(rep("beta",p))
  lpar = length(indpar)
  indpar2 = matrix(0, nrow = (p), ncol = (p))

  ##### Primeiras derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)

  ##### Primeiras derivadas de Ai para diferente de nu
  dAi = numeric(lpar)
  ai = as.numeric((1+t(lambda)%*%sFmat%*%Lambda%*%sFmat%*%lambda)^.5)
  bi = as.numeric((1+t(lambda)%*%lambda)^.5)
  Bi = as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)
  dAi[indpar=="beta"] = -1/ai*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda

  ##### Primeiras derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)

  ##### Primeiras derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = numeric(lpar)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  ##### Segundas Derivadas de Ai
  ddAi <- indpar2
  # dbeta.dbeta = 0
  ##### Segundas Derivadas de di
  dddi <- indpar2
  # dbeta.dbeta
  dddi[1:p,1:p] <- 2*t(x1)%*%sPsi%*%x1
  ##### Segundas Derivadas de ki
  ddki <- indpar2
  ntheta <- p
  for(i in 1:ntheta) for(j in 1:ntheta) ddki[i,j] = Iphi((ni+1)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddAi[i,j] -
    .5*Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*(ddi[j]+2*Ai*dAi[j])*dAi[i]-
    .5*IPhi((ni+2)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dddi[i,j]-
    .5*(Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dAi[j]-.5*IPhi((ni+4)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddi[j])*ddi[i]

  # Resultados
  sihat = -.5*dlogdpsi+1/ki*dki # primeira derivada
  hessianHat = indpar2
  for(i in 1:ntheta) for(j in 1:ntheta) hessianHat[i,j] = (1/ki)*ddki[i,j] - (1/ki^2)*dki[i]*dki[j] # Segunda derivada

  return(list(sihat = sihat, hessianHat = hessianHat))
}

expectBetasDEC <- function(jseq,y,x,z,time,beta1,sigmae,phiDEC,thetaDEC,D1,lambda,distr,nu){
  y1 = y[jseq]
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
  med<-x1%*%beta1 #+ c.*z1%*%Deltab
  Covmat<-CovDEC(phiDEC, thetaDEC, t1)
  sCovmat<-solve(Covmat)
  Sigma <- sigmae*Covmat
  Psi<-(z1)%*%(Gammab)%*%t(z1)+Sigma
  sPsi <- solve(Psi)
  di<-as.numeric(t(y1-med)%*%sPsi%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ai<-as.numeric(mutj/sqrt(Mtj2))
  sFmat = solve(Fmat)
  Lambda = solve(solve(D1)+ t(z1)%*%solve(Sigma)%*%z1)
  F.lista <- lapply(1:q2,F.r,q1=q1)

  if (distr=="sn"){
    e1 = e4 = e7 = 1# =matrix(1, nrow = (p), ncol = (p))
    e2 = e5 = Psi + med%*%t(med)
    e3 = e6 = med
  }
  if (distr=="st"){
    auxc <- function(nnu,r,p) (gamma((p+nnu)/2)*gamma((nnu+2*r)/2))/(nnu^r*gamma(nnu/2)*gamma((p+nnu+2*r)/2))
    e1 = (ni+nu)*auxc(nu,1,ni)
    e2 = (ni+nu+2)*(ni+nu)*auxc(nu,2,ni)*(med%*%t(med) + (nu/(nu+2))*Psi)
    e3 = (ni+nu+2)*(ni+nu)*auxc(nu,2,ni)*med
    e4 = (ni+nu+2)*(ni+nu)*auxc(nu,2,ni)
    e5 = (ni+nu)^2*auxc(nu,2,ni)*(med%*%t(med) + (nu/(nu+2))*Psi)
    e6 = (ni+nu)^2*auxc(nu,2,ni)*med
    e7 = (ni+nu)^2*auxc(nu,2,ni)
  }
  #if (distr=="ssl"){}
  #if (distr=="scn"){}

  ## E1 = E(d2log(det(Psi))) = 0
  ## E2 = E(1/Ki*d2Ki)
  E2i = -e1*t(x1)%*%sPsi%*%x1 + t(x1)%*%sPsi%*%e2%*%sPsi%*%x1 -
    2*t(x1)%*%sPsi%*%e3%*%t(med)%*%sPsi%*%x1 + e4*t(x1)%*%sPsi%*%med%*%t(med)%*%sPsi%*%x1
  ## E3 E(1/Ki^2*dKi/dbeta*dKi/dbeta)
  E3i = t(x1)%*%sPsi%*%e5%*%sPsi%*%x1 - 2*t(x1)%*%sPsi%*%e6%*%t(med)%*%sPsi%*%x1+
    e7*t(x1)%*%sPsi%*%med%*%t(med)%*%sPsi%*%x1
  ## E4 = E(d2log(det(Psi))/dbeta*d2log(det(Psi))/dbeta) = 0
  ## E5 = E(d2log(det(Psi))/dbeta*1/Ki*dKi/dbeta) = 0

  return(list(E2i = E2i, E3i = E3i))
}

# SMN - AR(p) ------------------------------------------------------------------
derivatesARis <- function(jseq,y,x,z,time,beta1,sigmae,phiAR,D1,lambda,distr,nu){
  if (distr=="sn"|distr=="norm") c.=-sqrt(2/pi)
  if (distr=="st"|distr=="t") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ssl"|distr=="sl") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn"|distr=="cn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[jseq]
  t1 = time[jseq]
  p= ncol(x);q1=ncol(z);pAR=length(phiAR)
  q2 = q1*(q1+1)/2
  if(q1 != 0){ q3 = 1} else{q3 = 0}
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
  # theta = c(beta1,sigmae,phi,dd,lambda,nu) - para AR(p)
  indpar = c(rep("beta",p),"sigma",rep("phi",pAR),rep("dd",q2))
  lpar = length(indpar)
  indpar2 = matrix(0, nrow = (p+1+pAR+q2), ncol = (p+1+pAR+q2))

  ##### Primeiras derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)
  dlogdpsi[indpar=="sigma"] =traceM(sPsi%*%MniAR)
  for (i in 1:q2) dlogdpsi[indpar=="dd"][i] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                    Fmat%*%F.lista[[i]])%*%t(z1))
  #jacobAR <- jacobian(Mnp,phiAR,n=ni) #matrix(jacobAR[,1],ncol=ni)
  jacobARautocovs <- matrix(jacobian(autocovsAR,phiAR,n=max(t1))[t1,],ncol=pAR) #toeplitz(jacobARautocovs[,1])
  for (i in 1:pAR) dlogdpsi[indpar=="phi"][i] = sigmae*traceM(sPsi%*%toeplitz(jacobARautocovs[,i]))

  ##### Primeiras derivadas de Ai para diferente de nu
  dAi = numeric(lpar)
  ai = as.numeric((1+t(lambda)%*%sFmat%*%Lambda%*%sFmat%*%lambda)^.5)
  bi = as.numeric((1+t(lambda)%*%lambda)^.5)
  Bi = as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)
  dAi[indpar=="beta"] = -1/ai*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda
  dAi[indpar=="sigma"] = -1/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae^2)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda

  for (i in 1:q2) dAi[indpar=="dd"][i] = 1/ai*(t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 c./bi*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda)+
    1/ai^2*Ai/2*t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                       Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:pAR) dAi[indpar=="phi"][i] = -sigmae/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[,i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda

  ##### Primeiras derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)
  ddi[indpar=="sigma"] = -t(y1-med)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)
  for (i in 1:q2) ddi[indpar=="dd"][i] =-2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)
  for (i in 1:pAR) ddi[indpar=="phi"][i] = -sigmae*t(y1-med)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%(y1-med)

  ##### Primeiras derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = numeric(lpar)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  ##### Segundas Derivadas de log(det(psi))
  ddlogdpsi = indpar2
  # dsigmae.dsigmae
  ddlogdpsi[p+1, p+1] = traceM(-sPsi%*%MniAR%*%sPsi%*%MniAR)
  # dsigmae.dphi
  for (i in 1:pAR) ddlogdpsi[p+1,p+1+i]= ddlogdpsi[p+1+i,p+1] = traceM(sPsi%*%toeplitz(jacobARautocovs[,i]) -
                                                                         sigmae*sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%MniAR)
  # dsigmae.dalphar
  for (i in 1:q2) ddlogdpsi[p+1,p+1+pAR+i] = ddlogdpsi[p+1+pAR+i,p+1] = traceM(-sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                                               Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%MniAR)
  # dalpha.dalpha
  for (i in 1:q2) for (j in 1:q2) ddlogdpsi[p+1+pAR+i,p+1+pAR+j] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%F.lista[[j]] + F.lista[[j]]%*%F.lista[[i]])%*%t(z1)-
                                                                            sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%
                                                                            sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1))
  # dalpha.dphi
  for (i in 1:q2) for (j in 1:pAR) ddlogdpsi[p+1+pAR+i, p+1+j] = ddlogdpsi[p+1+j, p+1+pAR+i] = - traceM(sigmae*sPsi%*%toeplitz(jacobARautocovs[,j])%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                                                                                                                     Fmat%*%F.lista[[i]])%*%t(z1))
  # dphi.dphi
  n1 <- length(y1)
  hessARautocovs <- list() # Derivadas dos phis em cada lag
  for(k in 1:n1){ hessARautocovs[[k]] <- hessian(autocovsAR2, phiAR, n=max(t1), j = k)}

  for (i in 1:pAR) for (j in 1:pAR) ddlogdpsi[p+1+i,p+1+j] = sigmae*traceM(sPsi%*%toeplitz(selectHessian(i, j, n1, hessARautocovs)) -
                                                                             sigmae*sPsi%*%toeplitz(jacobARautocovs[, j])%*%sPsi%*%toeplitz(jacobARautocovs[, i]))
  ##### Segundas Derivadas de Ai
  ddAi <- indpar2
  # dbeta.dbeta = 0

  # dbeta.dsigmae
  ddAi[1:p, p+1] <- (1/ai)*t(x1)%*%sPsi%*%MniAR%*%sPsi%*%z1%*%Fmat%*%lambda+
    (1/(2*ai^3*sigmae^2))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda
  ddAi[p+1, 1:p] <- t(ddAi[1:p, p+1])
  # dbeta.dphi
  for (i in 1:pAR) ddAi[1:p, p+1+i] <- (sigmae/ai)*t(x1)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%z1%*%Fmat%*%lambda +
    (1/(2*ai^3*sigmae))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[,i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda
  for (i in 1:pAR) ddAi[p+1+i, 1:p] <- ddAi[1:p, p+1+i]

  #dbeta.dalphar
  for (i in 1:q2) ddAi[1:p, p+1+pAR+i] <- (1/ai)*t(x1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda-
    (1/ai)*t(x1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda - (1/(2*ai^3))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                                               Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat + sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:q2) ddAi[p+1+pAR+i,1:p] <- ddAi[1:p, p+1+pAR+i]

  # dsigma.dsigma
  ddAi[p+1,p+1] <- (2/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)+
    (1/(2*sigmae^2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="sigma"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda -
    (1/(sigmae^4*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(sigmae^3*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda +
    (1/(2*sigmae^4*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda

  #dsigma.dphi
  for (i in 1:pAR) ddAi[p+1,p+1+i] <- (-1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%(y1-med)+
    (2*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)+
    (1/(2*ai^3*sigmae))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[,i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="phi"][i]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(sigmae^3*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[,i])%*%sMniAR%*%z1%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[,i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^3*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[,i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda
  for (i in 1:pAR) ddAi[p+1+i,p+1] <- ddAi[p+1,p+1+i]

  #dsigma.dalphar
  for (i in 1:q2) ddAi[p+1, p+1+pAR+i] <- (-1/ai)*t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%z1%*%F.lista[[i]]%*%delta-
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                    Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+ sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="dd"][i]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%F.lista[[i]]%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%(sFmat%*%F.lista[[i]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[i]]%*%sFmat)%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%(sFmat%*%F.lista[[i]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[i]]%*%sFmat)%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%F.lista[[i]]%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^4))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                                                 Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+ sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:q2) ddAi[p+1+pAR+i, p+1] <- ddAi[p+1, p+1+pAR+i]

  #dphi.dphi
  for (i in 1:pAR) for (j in 1:pAR) ddAi[p+1+i, p+1+j] <- (-sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%
    (-sigmae*sPsi%*%toeplitz(jacobARautocovs[, j])%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi +
       sPsi%*%toeplitz(selectHessian(i, j, n1, hessARautocovs))%*%sPsi -
       sigmae*sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%toeplitz(jacobARautocovs[, j])%*%sPsi)%*%(y1-med)+
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, j])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*dAi[indpar=="phi"][j]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, j])%*%sMniAR%*%z1%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%(-sMniAR%*%toeplitz(jacobARautocovs[, j])%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR+
                                                                     sMniAR%*%toeplitz(selectHessian(i, j, n1, hessARautocovs))%*%sMniAR -
                                                                     sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%toeplitz(jacobARautocovs[, j])%*%sMniAR)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, j])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*ai^4*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, j])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda

  #dphi.dalphar
  for (i in 1:pAR) for(j in 1:q2) ddAi[p+1+i, p+1+pAR+j] <- (-sigmae/ai)*t(lambda)%*%F.lista[[j]]%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%(y1-med)+
    (sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[j]] + F.lista[[j]]%*%Fmat)%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%(y1-med)+
    (sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%z1%*%(Fmat%*%F.lista[[j]] + F.lista[[j]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c.*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%z1%*%F.lista[[j]]%*%delta-
    (sigmae/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%toeplitz(jacobARautocovs[, i])%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[j]]%*%sFmat%*%Lambda+
                                                                                                                                  Lambda%*%sFmat%*%F.lista[[j]]-
                                                                                                                                  Lambda%*%sFmat%*%(F.lista[[j]]%*%sFmat+
                                                                                                                                                      sFmat%*%F.lista[[j]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*dAi[indpar=="dd"][j]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%F.lista[[j]]%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%(sFmat%*%F.lista[[j]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[j]]%*%sFmat)%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%(sFmat%*%F.lista[[j]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[j]]%*%sFmat)%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%F.lista[[j]]%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sMniAR%*%toeplitz(jacobARautocovs[, i])%*%sMniAR%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%(F.lista[[j]]%*%sFmat%*%Lambda+
                                                                                                                                                                   Lambda%*%sFmat%*%F.lista[[j]]-
                                                                                                                                                                   Lambda%*%sFmat%*%(F.lista[[j]]%*%sFmat+
                                                                                                                                                                                       sFmat%*%F.lista[[j]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:pAR) for(j in 1:q2) ddAi[p+1+pAR+ j, p+1+i] <- ddAi[p+1+i, p+1+pAR+j]

  #dalpha.dalpha
  for(r in 1:q2) for(s in 1:q2) ddAi[p+1+pAR+r, p+1+pAR+s] <- (-1/ai)*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]]+F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)-
    (c./ai)*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[s]]%*%delta+
    (1/(2*ai^3))*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                             Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (c./ai)*t(delta)%*%F.lista[[s]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda+
    (c./ai)*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda-
    (c./(2*ai^3))*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                 Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/ai)*t(lambda)%*%F.lista[[s]]%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)-
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[r]]%*%F.lista[[s]] + F.lista[[s]]%*%F.lista[[r]])%*%t(z1)%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[s]]%*%delta-
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                                                                       Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda+
    (1/(2*ai^2))*dAi[indpar=="dd"][s]%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                                               Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%F.lista[[s]]%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                                                    Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%F.lista[[s]]%*%sFmat%*%lambda+
    (1/(2*ai^4))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                                                                                      Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda+
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(-F.lista[[r]]%*%sFmat%*%F.lista[[s]]%*%sFmat%*%Lambda+
                                             F.lista[[r]]%*%sFmat%*%Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda+
                                             Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%F.lista[[s]]%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%F.lista[[s]]%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]]%*%sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%F.lista[[s]]%*%sFmat%*%Lambda-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda)%*%sFmat%*%lambda

  ##### Segundas Derivadas de di
  dddi <- indpar2
  # dbeta.dbeta
  dddi[1:p,1:p] <- 2*t(x1)%*%sPsi%*%x1
  # dbeta.dsigmae
  dddi[1:p, p+1] <- dddi[p+1, 1:p] <- 2*t(x1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)
  # dbeta.dphi
  for(i in 1:pAR) dddi[1:p, p+1+i]  <- dddi[p+1+i, 1:p] <- 2*sigmae*t(x1)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%(y1-med)
  # dbeta.dalphar
  for(i in 1:q2) dddi[1:p, p+1+pAR+i] <- dddi[p+1+pAR+i,1:p] <- 2*t(x1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)+ 2*c.*t(x1)%*%sPsi%*%z1%*%F.lista[[i]]%*%delta
  # dsigma.dsigma
  dddi[p+1,p+1] <- 2*t(y1-med)%*%sPsi%*%MniAR%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)
  # dsigma.dphi
  for(i in 1:pAR) dddi[p+1,p+1+i] <- dddi[p+1+i,p+1] <-t(y1-med)%*%(sigmae*sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%MniAR%*%sPsi - sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi + sigmae*sPsi%*%MniAR%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi)%*%(y1-med)
  # dsigmae.dalpha
  for(i in 1:q2) dddi[p+1,p+1+pAR+i] <- dddi[p+1+pAR+i, p+1] <- c.*t(z1%*%F.lista[[i]]%*%delta)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med) + t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%MniAR%*%sPsi%*%(y1-med)+
    t(y1-med)%*%sPsi%*%MniAR%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med) + c.*t(y1-med)%*%sPsi%*%MniAR%*%sPsi%*%z1%*%F.lista[[i]]%*%delta
  # dphi.dphi
  for(i in 1:pAR) for(j in 1:pAR) dddi[p+1+i,p+1+j] <- dddi[p+1+j,p+1+i] <- sigmae*t(y1-med)%*%(sigmae*sPsi%*%toeplitz(jacobARautocovs[,j])%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi - sPsi%*%toeplitz(selectHessian(i, j, n1, hessARautocovs))%*%sPsi + sigmae*sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%toeplitz(jacobARautocovs[,j])%*%sPsi)%*%(y1-med)
  # dphi.dalpha
  for(i in 1:pAR) for(j in 1:q2) dddi[p+1+i, p+1+pAR+j] <- dddi[p+1+pAR+j, p+1+i] <- 2*c.*sigmae*t(y1-med)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%z1%*%F.lista[[j]]%*%delta+
    2*sigmae%*%t(y1-med)%*%sPsi%*%toeplitz(jacobARautocovs[,i])%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat + Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)
  # dalpha.dalpha
  for(i in 1:q2) for(j in 1:q2) dddi[p+1+pAR+i, p+1+pAR+j] <- dddi[p+1+pAR+j, p+1+pAR+j] <- 2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)+
    2*c.^2*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[j]]%*%delta+
    2*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)+
    2*c.*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%F.lista[[j]]%*%delta-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%F.lista[[j]] + F.lista[[j]]%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)

  ##### Segundas Derivadas de ki
  ddki <- indpar2
  ntheta <- p+1+pAR+q2
  for(i in 1:ntheta) for(j in 1:ntheta) ddki[i,j] = Iphi((ni+1)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddAi[i,j] -
    .5*Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*(ddi[j]+2*Ai*dAi[j])*dAi[i]-
    .5*IPhi((ni+2)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dddi[i,j]-
    .5*(Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dAi[j]-.5*IPhi((ni+4)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddi[j])*ddi[i]

  # Resultados
  sihat = -.5*dlogdpsi+1/ki*dki # primeira derivada
  hessianHat = indpar2
  for(i in 1:ntheta) for(j in 1:ntheta) hessianHat[i,j] = -.5*ddlogdpsi[i,j] + (1/ki)*ddki[i,j] - (1/ki^2)*dki[i]*dki[j] # Segunda derivada

  return(list(sihat = sihat, hessianHat = hessianHat))
}

# SMN -UNC ---------------------------------------------------------------------
derivatesUNCs <- function(jseq,y,x,z,time,beta1,sigmae,D1,lambda,distr,nu){
  if (distr=="sn"|distr=="norm") c.=-sqrt(2/pi)
  if (distr=="st"|distr=="t") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ssl"|distr=="sl") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn"|distr=="cn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[jseq]
  t1 = time[jseq]
  p= ncol(x);q1=ncol(z)
  q2 = q1*(q1+1)/2
  if(q1 != 0){ q3 = 1} else{q3 = 0}
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
  Lambda = solve(solve(D1)+ t(z1)%*%solve(Sigma)%*%z1)
  F.lista <- lapply(1:q2,F.r,q1=q1)
  indpar = c(rep("beta",p),"sigma",rep("dd",q2))
  lpar = length(indpar)
  indpar2 = matrix(0, nrow = (p+1+q2), ncol = (p+1+q2))

  ##### Primeiras derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)
  dlogdpsi[indpar=="sigma"] =traceM(sPsi)
  for (i in 1:q2) dlogdpsi[indpar=="dd"][i] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                    Fmat%*%F.lista[[i]])%*%t(z1))
  ##### Primeiras derivadas de Ai para diferente de nu
  dAi = numeric(lpar)
  ai = as.numeric((1+t(lambda)%*%sFmat%*%Lambda%*%sFmat%*%lambda)^.5)
  bi = as.numeric((1+t(lambda)%*%lambda)^.5)
  Bi = as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)
  dAi[indpar=="beta"] = -1/ai*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda
  dAi[indpar=="sigma"] = -1/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae^2)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda

  for (i in 1:q2) dAi[indpar=="dd"][i] = 1/ai*(t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 c./bi*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda)+
    1/ai^2*Ai/2*t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-
                                       Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda

  ##### Primeiras derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)
  ddi[indpar=="lambda"] = -2*c./bi*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="sigma"] = -t(y1-med)%*%sPsi%*%sPsi%*%(y1-med)
  for (i in 1:q2) ddi[indpar=="dd"][i] =-2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)

  ##### Primeiras derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = numeric(lpar)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  ##### Segundas Derivadas de log(det(psi))
  ddlogdpsi = indpar2
  # dsigmae.dsigmae
  ddlogdpsi[p+1, p+1] = traceM(-sPsi%*%sPsi)
  # dsigmae.dalphar
  for (i in 1:q2) ddlogdpsi[p+1,p+1+i] = ddlogdpsi[p+1+i,p+1] = traceM(-sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                                       Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi)
  # dalpha.dalpha
  for (i in 1:q2) for (j in 1:q2) ddlogdpsi[p+1+i,p+1+j] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%F.lista[[j]] + F.lista[[j]]%*%F.lista[[i]])%*%t(z1)-
                                                                    sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%
                                                                    sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1))
  ##### Segundas Derivadas de Ai
  ddAi <- indpar2
  # dbeta.dbeta = 0

  # dbeta.dsigmae
  ddAi[1:p, p+1] <- (1/ai)*t(x1)%*%sPsi%*%sPsi%*%z1%*%Fmat%*%lambda+
    (1/(2*ai^3*sigmae^2))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda
  ddAi[p+1, 1:p] <- t(ddAi[1:p, p+1])
  #dbeta.dalphar
  for (i in 1:q2) ddAi[1:p, p+1+i] <- (1/ai)*t(x1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda-
    (1/ai)*t(x1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda - (1/(2*ai^3))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                                               Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat + sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:q2) ddAi[p+1+i,1:p] <- ddAi[1:p, p+1+i]
  # dsigma.dsigma
  ddAi[p+1,p+1] <- (2/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%sPsi%*%(y1-med)+
    (1/(2*sigmae^2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="sigma"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda -
    (1/(sigmae^4*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(sigmae^3*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda +
    (1/(2*sigmae^4*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda

  #dsigma.dalphar
  for (i in 1:q2) ddAi[p+1, p+1+i] <- (-1/ai)*t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%z1%*%F.lista[[i]]%*%delta-
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                            Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+ sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="dd"][i]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%F.lista[[i]]%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%(sFmat%*%F.lista[[i]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[i]]%*%sFmat)%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%(sFmat%*%F.lista[[i]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[i]]%*%sFmat)%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%F.lista[[i]]%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^4))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                                        Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+ sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:q2) ddAi[p+1+i, p+1] <- ddAi[p+1, p+1+i]

  #dalpha.dalpha
  for(r in 1:q2) for(s in 1:q2) ddAi[p+1+r, p+1+s] <- (-1/ai)*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]]+F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)-
    (c./ai)*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[s]]%*%delta+
    (1/(2*ai^3))*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                             Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (c./ai)*t(delta)%*%F.lista[[s]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda+
    (c./ai)*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda-
    (c./(2*ai^3))*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                 Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/ai)*t(lambda)%*%F.lista[[s]]%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)-
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[r]]%*%F.lista[[s]] + F.lista[[s]]%*%F.lista[[r]])%*%t(z1)%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[s]]%*%delta-
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                                                                       Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda+
    (1/(2*ai^2))*dAi[indpar=="dd"][s]%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                                               Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%F.lista[[s]]%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                                                    Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%F.lista[[s]]%*%sFmat%*%lambda+
    (1/(2*ai^4))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                                                                                      Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda+
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(-F.lista[[r]]%*%sFmat%*%F.lista[[s]]%*%sFmat%*%Lambda+
                                             F.lista[[r]]%*%sFmat%*%Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda+
                                             Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%F.lista[[s]]%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%F.lista[[s]]%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]]%*%sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%F.lista[[s]]%*%sFmat%*%Lambda-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda)%*%sFmat%*%lambda
  ##### Segundas Derivadas de di
  dddi <- indpar2
  # dbeta.dbeta
  dddi[1:p,1:p] <- 2*t(x1)%*%sPsi%*%x1
  # dbeta.dsigmae
  dddi[1:p, p+1] <- dddi[p+1, 1:p] <- 2*t(x1)%*%sPsi%*%sPsi%*%(y1-med)
  # dbeta.dalphar
  for(i in 1:q2) dddi[1:p, p+1+i] <- dddi[p+1+i,1:p] <- 2*t(x1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)+ 2*c.*t(x1)%*%sPsi%*%z1%*%F.lista[[i]]%*%delta
  # dsigma.dsigma
  dddi[p+1,p+1] <- 2*t(y1-med)%*%sPsi%*%sPsi%*%sPsi%*%(y1-med)
  # dsigmae.dalpha
  for(i in 1:q2) dddi[p+1,p+1+i] <- dddi[p+1+i, p+1] <- c.*t(z1%*%F.lista[[i]]%*%delta)%*%sPsi%*%sPsi%*%(y1-med) + t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)+
    t(y1-med)%*%sPsi%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med) + c.*t(y1-med)%*%sPsi%*%sPsi%*%z1%*%F.lista[[i]]%*%delta
  # dalpha.dalpha
  for(i in 1:q2) for(j in 1:q2) dddi[p+1+i, p+1+j] <- dddi[p+1+j, p+1+j] <- 2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)+
    2*c.^2*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[j]]%*%delta+
    2*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)+
    2*c.*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%F.lista[[j]]%*%delta-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%F.lista[[j]] + F.lista[[j]]%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)
  ##### Segundas Derivadas de ki
  ddki <- indpar2
  ntheta <- p+1+q2
  for(i in 1:ntheta) for(j in 1:ntheta) ddki[i,j] = Iphi((ni+1)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddAi[i,j] -
    .5*Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*(ddi[j]+2*Ai*dAi[j])*dAi[i]-
    .5*IPhi((ni+2)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dddi[i,j]-
    .5*(Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dAi[j]-.5*IPhi((ni+4)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddi[j])*ddi[i]

  # Resultados
  sihat = -.5*dlogdpsi+1/ki*dki # primeira derivada
  hessianHat = indpar2
  for(i in 1:ntheta) for(j in 1:ntheta) hessianHat[i,j] = -.5*ddlogdpsi[i,j] + (1/ki)*ddki[i,j] - (1/ki^2)*dki[i]*dki[j] # Segunda derivada

  return(list(sihat = sihat, hessianHat = hessianHat))
}

# SMSN - DEC -------------------------------------------------------------------
derivatesDECs <- function(jseq,y,x,z,time,beta1,sigmae,phiDEC,thetaDEC,D1,lambda,distr,nu){
  if (distr=="sn"|distr=="norm") c.=-sqrt(2/pi)
  if (distr=="st"|distr=="t") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ssl"|distr=="sl") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn"|distr=="cn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[jseq]
  t1 = time[jseq]
  p= ncol(x);q1=ncol(z); #pAR= 2 #length(phiAR)
  q2 = q1*(q1+1)/2
  if(q1 != 0){ q3 = 1} else{q3 = 0}
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  ni = length(y1)
  Fmat = matrix.sqrt(D1)
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-Fmat%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  Covmat <- CovDEC(phiDEC,thetaDEC,t1)  #Covmat <- CovDEC(phiDEC, thetaDEC, t1)
  sCovmat<-solve(Covmat)   #sCovmat<-solve(Covmat)
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
  # indpar = c(rep("beta",p),"sigma",rep("phi",pAR),rep("dd",q2),rep("lambda",q1))
  indpar = c(rep("beta",p),"sigma","phi","theta",rep("dd",q2))
  lpar = length(indpar)
  indpar2 = matrix(0, nrow = (p+1+2+q2), ncol = (p+1+2+q2))

  ##### Primeiras derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)
  dlogdpsi[indpar=="sigma"] =traceM(sPsi%*%Covmat)
  for (i in 1:q2) dlogdpsi[indpar=="dd"][i] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                    Fmat%*%F.lista[[i]])%*%t(z1))
  dphiDEC <- dphiCovDEC(phiDEC, thetaDEC, t1)
  dthetaDEC <- dthetaCovDEC(phiDEC, thetaDEC, t1)
  dphiphiDEC <- dphiphiCovDEC(phiDEC, thetaDEC, t1)
  dthetathetaDEC <- dthetathetaCovDEC(phiDEC, thetaDEC, t1)
  dphithetaDEC <- dphithetaCovDEC(phiDEC, thetaDEC, t1)

  ##### Primeiras derivadas de Ai para diferente de nu
  dAi = numeric(lpar)
  ai = as.numeric((1+t(lambda)%*%sFmat%*%Lambda%*%sFmat%*%lambda)^.5)
  bi = as.numeric((1+t(lambda)%*%lambda)^.5)
  Bi = as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)
  dAi[indpar=="beta"] = -1/ai*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda
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

  ##### Primeiras derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)
  ddi[indpar=="sigma"] = -t(y1-med)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)
  for (i in 1:q2) ddi[indpar=="dd"][i] =-2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="phi"] = -sigmae*t(y1-med)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)
  ddi[indpar=="theta"] = -sigmae*t(y1-med)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)

  ##### Primeiras derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = numeric(lpar)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi

  ##### Segundas Derivadas de log(det(psi))
  ddlogdpsi = indpar2
  # dsigmae.dsigmae
  ddlogdpsi[p+1, p+1] = traceM(-sPsi%*%Covmat%*%sPsi%*%Covmat)
  # dsigmae.dphi
  ddlogdpsi[p+1,p+1+1]= ddlogdpsi[p+1+1,p+1] = traceM(sPsi%*%dphiDEC - sigmae*sPsi%*%dphiDEC%*%sPsi%*%Covmat)
  ddlogdpsi[p+2,p+2+2]= ddlogdpsi[p+2+2,p+2] = traceM(sPsi%*%dthetaDEC - sigmae*sPsi%*%dthetaDEC%*%sPsi%*%Covmat)

  # dsigmae.dalphar
  for (i in 1:q2) ddlogdpsi[p+1,p+1+2+i] = ddlogdpsi[p+1+2+i,p+1] = traceM(-sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                                           Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%Covmat)
  # dalpha.dalpha
  for (i in 1:q2) for (j in 1:q2) ddlogdpsi[p+1+2+i,p+1+2+j] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%F.lista[[j]] + F.lista[[j]]%*%F.lista[[i]])%*%t(z1)-
                                                                        sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%
                                                                        sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1))
  # dalpha.dphi
  for (i in 1:q2) ddlogdpsi[p+1+1+i, p+1+1] = ddlogdpsi[p+1+1, p+1+1+i] = -traceM(sigmae*sPsi%*%dphiDEC%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+ Fmat%*%F.lista[[i]])%*%t(z1))
  for (i in 1:q2) ddlogdpsi[p+1+2+i, p+1+2] = ddlogdpsi[p+1+2, p+1+2+i] = -traceM(sigmae*sPsi%*%dthetaDEC%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1))

  # dphi.dphi
  ddlogdpsi[p+1+1,p+1+1] = sigmae*traceM(sPsi%*%dphiphiDEC - sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC)
  ddlogdpsi[p+1+2,p+1+2] = sigmae*traceM(sPsi%*%dthetathetaDEC - sigmae*sPsi%*%dthetaDEC%*%sPsi%*%dthetaDEC)
  ddlogdpsi[p+1+1,p+1+2] = ddlogdpsi[p+1+2,p+1+1] = sigmae*traceM(sPsi%*%dphithetaDEC - sigmae*sPsi%*%dphiDEC%*%sPsi%*%dthetaDEC)

  ##### Segundas Derivadas de Ai
  ddAi <- indpar2
  # dbeta.dbeta = 0

  # dbeta.dsigmae
  ddAi[1:p, p+1] <- (1/ai)*t(x1)%*%sPsi%*%Covmat%*%sPsi%*%z1%*%Fmat%*%lambda+
    (1/(2*ai^3*sigmae^2))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda
  ddAi[p+1, 1:p] <- t(ddAi[1:p, p+1])
  # dbeta.dphi
  ddAi[1:p, p+1+1] <- (sigmae/ai)*t(x1)%*%sPsi%*%dphiDEC%*%sPsi%*%z1%*%Fmat%*%lambda +
    (1/(2*ai^3*sigmae))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda
  ddAi[1:p, p+1+2] <- (sigmae/ai)*t(x1)%*%sPsi%*%dthetaDEC%*%sPsi%*%z1%*%Fmat%*%lambda +
    (1/(2*ai^3*sigmae))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  ddAi[p+1+1, 1:p] <- ddAi[1:p, p+1+1]
  ddAi[p+1+2, 1:p] <- ddAi[1:p, p+1+2]

  #dbeta.dalphar
  for (i in 1:q2) ddAi[1:p, p+1+2+i] <- (1/ai)*t(x1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda-
    (1/ai)*t(x1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda - (1/(2*ai^3))*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                                               Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat + sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:q2) ddAi[p+1+2+i,1:p] <- ddAi[1:p, p+1+2+i]

  # dsigma.dsigma
  ddAi[p+1,p+1] <- (2/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)+
    (1/(2*sigmae^2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="sigma"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda -
    (1/(sigmae^4*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(sigmae^3*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda +
    (1/(2*sigmae^4*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  #dsigma.dphi
  ddAi[p+1,p+1+1] <- (-1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)+
    (2*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)+
    (1/(2*ai^3*sigmae))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="phi"][i]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(sigmae^3*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^3*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda
  ddAi[p+1+1,p+1] <- ddAi[p+1,p+1+1]

  ddAi[p+1,p+1+2] <- (-1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)+
    (2*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)+
    (1/(2*ai^3*sigmae))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="theta"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(sigmae^3*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^3*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda
  ddAi[p+1+2,p+1] <- ddAi[p+1,p+1+2]

  #dsigma.dalphar
  for (i in 1:q2) ddAi[p+1, p+1+2+i] <- (-1/ai)*t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%z1%*%(Fmat%*%F.lista[[i]] + F.lista[[i]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%z1%*%F.lista[[i]]%*%delta-
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                     Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+ sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*dAi[indpar=="dd"][i]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%F.lista[[i]]%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%(sFmat%*%F.lista[[i]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[i]]%*%sFmat)%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%(sFmat%*%F.lista[[i]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[i]]%*%sFmat)%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae^2*ai^2))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%F.lista[[i]]%*%sFmat%*%lambda-
    (1/(2*sigmae^2*ai^4))*Ai*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[i]]-
                                                                                                                                  Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+ sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for (i in 1:q2) ddAi[p+1+2+i, p+1] <- ddAi[p+1, p+1+2+i]

  #dphi.dphi
  ddAi[p+1+1, p+1+1] <- (-sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%
    (-sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC%*%sPsi + sPsi%*%dphiphiDEC%*%sPsi -
       sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC%*%sPsi)%*%(y1-med)+ (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*dAi[indpar=="phi"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%
    (-sCovmat%*%dphiDEC%*%sCovmat%*%dphiDEC%*%sCovmat+ sCovmat%*%dphiphiDEC%*%sCovmat - sCovmat%*%dphiDEC%*%sCovmat%*%dphiDEC%*%sCovmat)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*ai^4*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  ddAi[p+1+2, p+1+2] <- (-sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%
    (-sigmae*sPsi%*%dthetaDEC%*%sPsi%*%dthetaDEC%*%sPsi + sPsi%*%dthetathetaDEC%*%sPsi -
       sigmae*sPsi%*%dthetaDEC%*%sPsi%*%dthetaDEC%*%sPsi)%*%(y1-med)+ (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*dAi[indpar=="theta"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%
    (-sCovmat%*%dthetaDEC%*%sCovmat%*%dthetaDEC%*%sCovmat+ sCovmat%*%dthetathetaDEC%*%sCovmat - sCovmat%*%dthetaDEC%*%sCovmat%*%dthetaDEC%*%sCovmat)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*ai^4*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  ddAi[p+1+1, p+1+2] <- (-sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%
    (-sigmae*sPsi%*%dthetaDEC%*%sPsi%*%dphiDEC%*%sPsi + sPsi%*%dphithetaDEC%*%sPsi - sigmae*sPsi%*%dphiDEC%*%sPsi%*%dthetaDEC%*%sPsi)%*%(y1-med)+
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*dAi[indpar=="theta"]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%
    (-sCovmat%*%dthetaDEC%*%sCovmat%*%dphiDEC%*%sCovmat+sCovmat%*%dphithetaDEC%*%sCovmat - sCovmat%*%dphiDEC%*%sCovmat%*%dthetaDEC%*%sCovmat)%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*ai^2*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*ai^4*sigmae^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda

  ddAi[p+1+2, p+1+1] = ddAi[p+1+1, p+1+2]
  #dphi.dalphar
  for(j in 1:q2) ddAi[p+1+1, p+1+2+j] <- (-sigmae/ai)*t(lambda)%*%F.lista[[j]]%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)+
    (sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[j]] + F.lista[[j]]%*%Fmat)%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)+
    (sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%z1%*%(Fmat%*%F.lista[[j]] + F.lista[[j]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c.*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%z1%*%F.lista[[j]]%*%delta-
    (sigmae/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[j]]%*%sFmat%*%Lambda+
                                                                                                           Lambda%*%sFmat%*%F.lista[[j]]-
                                                                                                           Lambda%*%sFmat%*%(F.lista[[j]]%*%sFmat+
                                                                                                                               sFmat%*%F.lista[[j]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*dAi[indpar=="dd"][j]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%F.lista[[j]]%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%(sFmat%*%F.lista[[j]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[j]]%*%sFmat)%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%(sFmat%*%F.lista[[j]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[j]]%*%sFmat)%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%F.lista[[j]]%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dphiDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%(F.lista[[j]]%*%sFmat%*%Lambda+
                                                                                                                                              Lambda%*%sFmat%*%F.lista[[j]]-
                                                                                                                                              Lambda%*%sFmat%*%(F.lista[[j]]%*%sFmat+
                                                                                                                                                                  sFmat%*%F.lista[[j]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for(j in 1:q2) ddAi[p+1+2, p+1+2+j] <- (-sigmae/ai)*t(lambda)%*%F.lista[[j]]%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)+
    (sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[j]] + F.lista[[j]]%*%Fmat)%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)+
    (sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%z1%*%(Fmat%*%F.lista[[j]] + F.lista[[j]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c.*sigmae/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%z1%*%F.lista[[j]]%*%delta-
    (sigmae/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[j]]%*%sFmat%*%Lambda+
                                                                                                             Lambda%*%sFmat%*%F.lista[[j]]-
                                                                                                             Lambda%*%sFmat%*%(F.lista[[j]]%*%sFmat+
                                                                                                                                 sFmat%*%F.lista[[j]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*dAi[indpar=="dd"][j]%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%F.lista[[j]]%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%(sFmat%*%F.lista[[j]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[j]]%*%sFmat)%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%(sFmat%*%F.lista[[j]]%*%sFmat%*%sFmat + sFmat%*%sFmat%*%F.lista[[j]]%*%sFmat)%*%Lambda%*%sFmat%*%lambda+
    (1/(2*sigmae*ai^2))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%F.lista[[j]]%*%sFmat%*%lambda-
    (1/(2*sigmae*ai^4))*Ai%*%t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%sCovmat%*%dthetaDEC%*%sCovmat%*%z1%*%Lambda%*%sFmat%*%lambda%*%t(lambda)%*%(F.lista[[j]]%*%sFmat%*%Lambda+
                                                                                                                                                Lambda%*%sFmat%*%F.lista[[j]]-
                                                                                                                                                Lambda%*%sFmat%*%(F.lista[[j]]%*%sFmat+
                                                                                                                                                                    sFmat%*%F.lista[[j]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  for(j in 1:q2) ddAi[p+1+2+ j, p+1+1] <- ddAi[p+1+i, p+1+2+1]
  for(j in 1:q2) ddAi[p+1+2+ j, p+1+2] <- ddAi[p+1+i, p+1+2+2]

  #dalpha.dalpha
  for(r in 1:q2) for(s in 1:q2) ddAi[p+1+2+r, p+1+2+s] <- (-1/ai)*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]]+F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)-
    (c./ai)*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[s]]%*%delta+
    (1/(2*ai^3))*t(lambda)%*%F.lista[[r]]%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                             Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (c./ai)*t(delta)%*%F.lista[[s]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda+
    (c./ai)*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda-
    (c./(2*ai^3))*t(Deltab)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[r]]%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                 Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/ai)*t(lambda)%*%F.lista[[s]]%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)-
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[r]]%*%F.lista[[s]] + F.lista[[s]]%*%F.lista[[r]])%*%t(z1)%*%sPsi%*%(y1-med)+
    (1/ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[s]] + F.lista[[s]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)+
    (c./ai)*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%z1%*%F.lista[[s]]%*%delta-
    (1/(2*ai^3))*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(Fmat%*%F.lista[[r]] + F.lista[[r]]%*%Fmat)%*%t(z1)%*%sPsi%*%(y1-med)%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                                                                       Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda+
    (1/(2*ai^2))*dAi[indpar=="dd"][s]%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                                               Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%F.lista[[s]]%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                                                    Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda-
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%F.lista[[s]]%*%sFmat%*%lambda+
    (1/(2*ai^4))*Ai%*%t(lambda)%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat+sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda%*%t(lambda)%*%sFmat%*%(F.lista[[s]]%*%sFmat%*%Lambda + Lambda%*%sFmat%*%F.lista[[s]]-
                                                                                                                                                                      Lambda%*%sFmat%*%(F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda+
    (1/(2*ai^2))*Ai%*%t(lambda)%*%sFmat%*%(-F.lista[[r]]%*%sFmat%*%F.lista[[s]]%*%sFmat%*%Lambda+
                                             F.lista[[r]]%*%sFmat%*%Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda+
                                             Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%sFmat%*%F.lista[[s]]%*%sFmat%*%F.lista[[r]]-
                                             Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%F.lista[[s]]%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat%*%F.lista[[s]]%*%sFmat+sFmat%*%F.lista[[s]]%*%sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda+
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%F.lista[[s]]%*%sFmat%*%Lambda-
                                             Lambda%*%sFmat%*%(F.lista[[r]]%*%sFmat + sFmat%*%F.lista[[r]])%*%sFmat%*%Lambda%*%(sFmat%*%sFmat%*%F.lista[[s]]%*%sFmat + sFmat%*%F.lista[[s]]%*%sFmat%*%sFmat)%*%Lambda)%*%sFmat%*%lambda
  ##### Segundas Derivadas de di
  dddi <- indpar2
  # dbeta.dbeta
  dddi[1:p,1:p] <- 2*t(x1)%*%sPsi%*%x1
  # dbeta.dsigmae
  dddi[1:p, p+1] <- dddi[p+1, 1:p] <- 2*t(x1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)
  # dbeta.dphi
  dddi[1:p, p+1+1]  <- dddi[p+1+1, 1:p] <- 2*sigmae*t(x1)%*%sPsi%*%dphiDEC%*%sPsi%*%(y1-med)
  dddi[1:p, p+1+2]  <- dddi[p+1+2, 1:p] <- 2*sigmae*t(x1)%*%sPsi%*%dthetaDEC%*%sPsi%*%(y1-med)
  # dbeta.dalphar
  for(i in 1:q2) dddi[1:p, p+1+2+i] <- dddi[p+1+2+i,1:p] <- 2*t(x1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)+ 2*c.*t(x1)%*%sPsi%*%z1%*%F.lista[[i]]%*%delta
  # dsigma.dsigma
  dddi[p+1,p+1] <- 2*t(y1-med)%*%sPsi%*%Covmat%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)
  # dsigma.dphi
  dddi[p+1,p+1+1] <- dddi[p+1+1,p+1] <-t(y1-med)%*%(sigmae*sPsi%*%dphiDEC%*%sPsi%*%Covmat%*%sPsi - sPsi%*%dphiDEC%*%sPsi + sigmae*sPsi%*%Covmat%*%sPsi%*%dphiDEC%*%sPsi)%*%(y1-med)
  dddi[p+1,p+1+2] <- dddi[p+1+2,p+1] <-t(y1-med)%*%(sigmae*sPsi%*%dthetaDEC%*%sPsi%*%Covmat%*%sPsi - sPsi%*%dthetaDEC%*%sPsi + sigmae*sPsi%*%Covmat%*%sPsi%*%dthetaDEC%*%sPsi)%*%(y1-med)
  # dsigmae.dalpha
  for(i in 1:q2) dddi[p+1,p+1+2+i] <- dddi[p+1+2+i, p+1] <- c.*t(z1%*%F.lista[[i]]%*%delta)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med) + t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%Covmat%*%sPsi%*%(y1-med)+
    t(y1-med)%*%sPsi%*%Covmat%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat + Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med) + c.*t(y1-med)%*%sPsi%*%Covmat%*%sPsi%*%z1%*%F.lista[[i]]%*%delta
  # dphi.dphi
  dddi[p+1+1,p+1+1] <- sigmae*t(y1-med)%*%(sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC%*%sPsi - sPsi%*%dphiphiDEC%*%sPsi + sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC%*%sPsi)%*%(y1-med)
  dddi[p+1+2,p+1+2] <- sigmae*t(y1-med)%*%(sigmae*sPsi%*%dthetaDEC%*%sPsi%*%dthetaDEC%*%sPsi - sPsi%*%dthetathetaDEC%*%sPsi + sigmae*sPsi%*%dthetaDEC%*%sPsi%*%dthetaDEC%*%sPsi)%*%(y1-med)
  dddi[p+1+1,p+1+2] <- dddi[p+1+2,p+1+1] <- sigmae*t(y1-med)%*%(sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC%*%sPsi - sPsi%*%dphithetaDEC%*%sPsi + sigmae*sPsi%*%dphiDEC%*%sPsi%*%dphiDEC%*%sPsi)%*%(y1-med)

  # dphi.dalpha
  for(j in 1:q2) dddi[p+1+1, p+1+2+j] <- dddi[p+1+2+j, p+1+1] <- 2*c.*sigmae*t(y1-med)%*%sPsi%*%dphiDEC%*%sPsi%*%z1%*%F.lista[[j]]%*%delta+
    2*sigmae%*%t(y1-med)%*%sPsi%*%dphiDEC%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat + Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)
  for(j in 1:q2) dddi[p+1+2, p+1+2+j] <- dddi[p+1+2+j, p+1+2] <- 2*c.*sigmae*t(y1-med)%*%sPsi%*%dthetaDEC%*%sPsi%*%z1%*%F.lista[[j]]%*%delta+
    2*sigmae%*%t(y1-med)%*%sPsi%*%dthetaDEC%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat + Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)
  # dalpha.dalpha
  for(i in 1:q2) for(j in 1:q2) dddi[p+1+2+i, p+1+2+j] <- dddi[p+1+2+j, p+1+2+j] <- 2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)+
    2*c.^2*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%z1%*%F.lista[[j]]%*%delta+
    2*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[j]]%*%Fmat+Fmat%*%F.lista[[j]])%*%t(z1)%*%sPsi%*%(y1-med)+
    2*c.*t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%z1%*%F.lista[[j]]%*%delta-
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%F.lista[[j]] + F.lista[[j]]%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)

  ##### Segundas Derivadas de ki ---- CHECK
  ddki <- indpar2
  ntheta <- p+1+2+q2
  for(i in 1:ntheta) for(j in 1:ntheta) ddki[i,j] = Iphi((ni+1)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddAi[i,j] -
    .5*Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*(ddi[j]+2*Ai*dAi[j])*dAi[i]-
    .5*IPhi((ni+2)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dddi[i,j]-
    .5*(Iphi((ni+3)/2, di = di, Ai=Ai, distr=distr, nu=nu)*dAi[j]-.5*IPhi((ni+4)/2, di = di, Ai=Ai, distr=distr, nu=nu)*ddi[j])*ddi[i]

  # Resultados
  sihat = -.5*dlogdpsi+1/ki*dki # primeira derivada
  hessianHat = indpar2
  for(i in 1:ntheta) for(j in 1:ntheta) hessianHat[i,j] = -.5*ddlogdpsi[i,j] + (1/ki)*ddki[i,j] - (1/ki^2)*dki[i]*dki[j] # Segunda derivada

  return(list(sihat = sihat, hessianHat = hessianHat))
}
