#fixpt functions for SMSN
fixpt.skewAR <- function(thetav,y,x,z,time,ind,distr,pAR,lb,lu,parallelphi,parallelnu,
                         diagD,skewind){
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2;N<-length(y)
  m <- n_distinct(ind)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  Gammab <- Dmatrix(thetav[(p+2):(p+1+q2)])
  Deltab<-thetav[(p+2+q2):(p+1+q2+q1)]
  piAR<-thetav[(p+2+q2+q1):(p+1+q2+q1+pAR)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+1+q2+q1+pAR))]
  #
  res_emj = revert_list(tapply(1:N,ind,emjAR,y=y, x=x, z=z,time=time, beta1=beta1, Gammab=Gammab,
                               Deltab=Deltab, sigmae=sigmae,piAR=piAR, distr=distr,nu=nu))
  sum1 = Reduce("+",res_emj$sum1)
  sum2 = Reduce("+",res_emj$sum2)
  sum3 = sum(unlist(res_emj$sum3))
  sum4 = Reduce("+",res_emj$sum4)
  sum5 = Reduce("+",res_emj$sum5)
  ut2j = unlist(res_emj$ut2j,use.names = F)
  #
  beta1<-solve(sum1)%*%sum2
  sigmae<-as.numeric(sum3)/N
  Gammab<-sum4/m
  Deltab<-sum5/sum(ut2j)
  #
  D1<-Gammab+Deltab%*%t(Deltab);sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))
  #special cases
  if (diagD||(sum(skewind)<q1)) {
    lambda<-lambda*skewind
    if (diagD) D1 <-diag(diag(D1))
    delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
    Deltab<-matrix.sqrt(D1)%*%delta
    Gammab<-D1-Deltab%*%t(Deltab)
  }
  #
  if (parallelphi) {
    piAR<- optimParallel(piAR,lcAR,gr = NULL,method = "L-BFGS-B", lower =rep(-.9999,pAR),
                 upper = rep(.9999,pAR),control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                 y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  } else {
    piAR<- optim(piAR,lcAR,gr = NULL,method = "L-BFGS-B", lower =rep(-.9999,pAR),
           upper = rep(.9999,pAR),control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
           y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  }
  #
  logvero1<-function(nu){logveroAR(y, x, z,time, ind, beta1, sigmae,estphit(piAR), D1, lambda, distr, nu)}
  
  if (distr=="sn"){ nu<-NULL} else {
    if (parallelnu) {
      nu <- optimParallel(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,
                  control = list(fnscale=-1))$par
    } else {
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,
                  control = list(fnscale=-1))$par 
    }
  }
  #
  c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,piAR,nu)
}
fixpt.skewUNC <- function(thetav,y,x,z,ind,distr,lb,lu,parallelnu,diagD,skewind){
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2;N<-length(y)
  m <- n_distinct(ind)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  Gammab <- Dmatrix(thetav[(p+2):(p+1+q2)])
  Deltab<-thetav[(p+2+q2):(p+1+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+1+q2+q1))]
  #
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
  sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))
  #special cases
  if (diagD||(sum(skewind)<q1)) {
    lambda<-lambda*skewind
    if (diagD) D1 <-diag(diag(D1))
    delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
    Deltab<-matrix.sqrt(D1)%*%delta
    Gammab<-D1-Deltab%*%t(Deltab)
  }
  #
  logvero1<-function(nu){logvero(y, x, z, ind, beta1, sigmae, D1, lambda, distr, nu)}
  
  if (distr=="sn"){ nu<-NULL} else {
    if (parallelnu) {
      nu <- optimParallel(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, 
                          upper = lu,control = list(fnscale=-1))$par
    } else {
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,
                  control = list(fnscale=-1))$par
    }
  }
  #
  c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,nu)
}
fixpt.skewCS <- function(thetav,y,x,z,ind,distr,lb,lu,parallelphi,parallelnu,diagD,skewind){
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2;N<-length(y)
  m <- n_distinct(ind)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  Gammab <- Dmatrix(thetav[(p+2):(p+1+q2)])
  Deltab<-thetav[(p+2+q2):(p+1+q2+q1)]
  phiCS<-thetav[(p+2+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+2+q2+q1))]
  #
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
  D1<-Gammab+Deltab%*%t(Deltab);sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))
  #special cases
  if (diagD||(sum(skewind)<q1)) {
    lambda<-lambda*skewind
    if (diagD) D1 <-diag(diag(D1))
    delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
    Deltab<-matrix.sqrt(D1)%*%delta
    Gammab<-D1-Deltab%*%t(Deltab)
  }
  #
  if (parallelphi) {
    phiCS <- optimParallel(phiCS,lcCS,gr = NULL,method = "L-BFGS-B", lower =0,
                   upper = .9999,control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  } else{
    phiCS <- optim(phiCS,lcCS,gr = NULL,method = "L-BFGS-B", lower =0,
                   upper = .9999,control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  }
  #
  logvero1<-function(nu){logveroCS(y, x, z, ind, beta1, sigmae,phiCS, D1, lambda, distr, nu)}
  
  if (distr=="sn"){ nu<-NULL} else {
    if (parallelnu) {
      nu <- optimParallel(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, 
                          upper = lu,control = list(fnscale=-1))$par
    } else{
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,
                  control = list(fnscale=-1))$par
    }
  }
  #
  c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,phiCS,nu)
}
fixpt.skewDEC <- function(thetav,y,x,z,time,ind,distr,lb,lu,luDEC,parallelphi,
                          parallelnu,diagD,skewind){
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2;N<-length(y)
  m <- n_distinct(ind)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  Gammab <- Dmatrix(thetav[(p+2):(p+1+q2)])
  Deltab<-thetav[(p+2+q2):(p+1+q2+q1)]
  phiDEC<-thetav[(p+2+q2+q1)]
  thetaDEC<-thetav[(p+3+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+3+q2+q1))]
  #
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
  D1<-Gammab+Deltab%*%t(Deltab);sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))
  #special cases
  if (diagD||(sum(skewind)<q1)) {
    lambda<-lambda*skewind
    if (diagD) D1 <-diag(diag(D1))
    delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
    Deltab<-matrix.sqrt(D1)%*%delta
    Gammab<-D1-Deltab%*%t(Deltab)
  }
  #
  if (parallelphi) {
    parDEC<- optimParallel(c(phiDEC,thetaDEC),lcDEC,gr = NULL,method = "L-BFGS-B", lower =rep(0.0001,2),
                   upper = c(.9999,luDEC),control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  } else{
    parDEC<- optim(c(phiDEC,thetaDEC),lcDEC,gr = NULL,method = "L-BFGS-B", lower =rep(0.0001,2),
                   upper = c(.9999,luDEC),control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  }
  phiDEC<-parDEC[1]
  thetaDEC<-parDEC[2]
  #
  logvero1<-function(nu){logveroDEC(y, x, z,time, ind, beta1, sigmae,phiDEC,thetaDEC, D1, lambda, distr, nu)}
  
  if (distr=="sn"){ nu<-NULL} else {
    if (parallelnu) {
      nu <- optimParallel(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par
    } else{
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par 
    }
    }
  #
  c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,phiDEC,thetaDEC,nu)
}
fixpt.skewCAR1 <- function(thetav,y,x,z,time,ind,distr,lb,lu,parallelphi,parallelnu,
                           diagD,skewind){
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2;N<-length(y)
  m <- n_distinct(ind)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  Gammab <- Dmatrix(thetav[(p+2):(p+1+q2)])
  Deltab<-thetav[(p+2+q2):(p+1+q2+q1)]
  phiDEC<-thetav[(p+2+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+2+q2+q1))]
  #
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
  D1<-Gammab+Deltab%*%t(Deltab);sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))
  #special cases
  if (diagD||(sum(skewind)<q1)) {
    lambda<-lambda*skewind
    if (diagD) D1 <-diag(diag(D1))
    delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
    Deltab<-matrix.sqrt(D1)%*%delta
    Gammab<-D1-Deltab%*%t(Deltab)
  }
  #
  if (parallelphi){
    phiDEC<- optimParallel(phiDEC,lcCAR1,gr = NULL,method = "L-BFGS-B", lower =0.0001,
                   upper = .9999,control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  } else{
    phiDEC<- optim(phiDEC,lcCAR1,gr = NULL,method = "L-BFGS-B", lower =0.0001,
                   upper = .9999,control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  }
  #
  logvero1<-function(nu){logveroCAR1(y, x, z,time, ind, beta1, sigmae,phiDEC, D1, lambda, distr, nu)}
  
  if (distr=="sn"){ nu<-NULL} else {
    if (parallelnu) {
      nu <- optimParallel(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, 
                          upper = lu,control = list(fnscale=-1))$par
    } else{
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,
                  control = list(fnscale=-1))$par
    }
  }
  #
  c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,phiDEC,nu)
}

#fixpt functions for SMN
fixpt.AR <- function(thetav,y,x,z,time,ind,distr,pAR,lb,lu,parallelphi,parallelnu,diagD){
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2;N<-length(y)
  m <- n_distinct(ind)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  D1 <- Dmatrix(thetav[(p+2):(p+1+q2)])
  piAR<-thetav[(p+2+q2):(p+1+q2+pAR)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+1+q2+pAR))]
  #
  res_emj = revert_list(tapply(1:N,ind,emjARs,y=y, x=x, z=z,time=time, beta1=beta1, D1=D1,
                               sigmae=sigmae,piAR=piAR, distr=distr,nu=nu))
  sum1 = Reduce("+",res_emj$sum1)
  sum2 = Reduce("+",res_emj$sum2)
  sum3 = sum(unlist(res_emj$sum3))
  sum4 = Reduce("+",res_emj$sum4)
  
  beta1<-solve(sum1)%*%sum2
  sigmae<-as.numeric(sum3)/N
  D1<-sum4/m
  if (diagD) D1 <-diag(diag(D1))
  #
  if (parallelphi){
    piAR<- optimParallel(piAR,lcAR,gr = NULL,method = "L-BFGS-B", lower =rep(-.9999,pAR),
                 upper = rep(.9999,pAR),control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                 y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  } else{
    piAR<- optim(piAR,lcAR,gr = NULL,method = "L-BFGS-B", lower =rep(-.9999,pAR),
                 upper = rep(.9999,pAR),control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                 y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  }
  #
  logvero1<-function(nu){logveroARs(y = y,x = x, z = z,time = time,ind =  ind,
                                    beta1 =  beta1,sigmae = sigmae,phiAR = estphit(piAR),
                                    D1 =  D1, distr, nu)}
  
  if (distr=="sn"){ nu<-NULL} else {
    if (parallelnu) {
      nu <- optimParallel(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, 
                          upper = lu,control = list(fnscale=-1))$par
    } else{
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,
                  control = list(fnscale=-1))$par
    }
  }
  #
  c(beta1,sigmae,D1[upper.tri(D1, diag = T)],piAR,nu)
}
fixpt.UNC <- function(thetav,y,x,z,ind,distr,lb,lu,parallelnu,diagD){
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2;N<-length(y)
  m <- n_distinct(ind)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  D1 <- Dmatrix(thetav[(p+2):(p+1+q2)])
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+1+q2))]
  #
  res_emj = revert_list(tapply(1:N,ind,emjs,y=y, x=x, z=z, beta1=beta1, D1=D1,
                               sigmae=sigmae, distr=distr,nu=nu))
  sum1 = Reduce("+",res_emj$sum1)
  sum2 = Reduce("+",res_emj$sum2)
  sum3 = sum(unlist(res_emj$sum3))
  sum4 = Reduce("+",res_emj$sum4)
  uj = unlist(res_emj$uj,use.names = F)
  #if (calcbi) bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
  
  beta1<-solve(sum1)%*%sum2
  sigmae<-as.numeric(sum3)/N
  D1<-sum4/m
  if (diagD) D1 <-diag(diag(D1))
  logvero1<-function(nu){logveros(y, x, z, ind, beta1, sigmae, D1, distr, nu)}
  
  if (distr=="sn"){ nu<-NULL} else {
    if (parallelnu) {
      nu <- optimParallel(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, 
                          upper = lu,control = list(fnscale=-1))$par
    } else{
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,
                  control = list(fnscale=-1))$par
    }
  }
  #
  c(beta1,sigmae,D1[upper.tri(D1, diag = T)],nu)
}
fixpt.CS <- function(thetav,y,x,z,ind,distr,lb,lu,parallelphi,parallelnu,diagD){
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2;N<-length(y)
  m <- n_distinct(ind)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  D1 <- Dmatrix(thetav[(p+2):(p+1+q2)])
  phiCS<-thetav[(p+2+q2)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+2+q2))]
  #
  res_emj = revert_list(tapply(1:N,ind,emjCSs,y=y, x=x, z=z, beta1=beta1, D1=D1,
                               sigmae=sigmae,phiCS=phiCS, distr=distr,nu=nu))
  sum1 = Reduce("+",res_emj$sum1)
  sum2 = Reduce("+",res_emj$sum2)
  sum3 = sum(unlist(res_emj$sum3))
  sum4 = Reduce("+",res_emj$sum4)
  
  #if (calcbi) bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
  beta1<-solve(sum1)%*%sum2
  sigmae<-as.numeric(sum3)/N
  D1<-sum4/m
  if (diagD) D1 <-diag(diag(D1))
  #
  if (parallelphi) {
    phiCS <- optimParallel(phiCS,lcCS,gr = NULL,method = "L-BFGS-B", lower =0,
                   upper = .9999,control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  } else {
    phiCS <- optim(phiCS,lcCS,gr = NULL,method = "L-BFGS-B", lower =0,
                   upper = .9999,control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  }
  #
  logvero1<-function(nu){logveroCSs(y, x, z, ind, beta1, sigmae,phiCS, D1, distr, nu)}
  
  if (distr=="sn"){ nu<-NULL} else {
    if (parallelnu) {
      nu <- optimParallel(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,
                  control = list(fnscale=-1))$par
    } else {
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,
                  control = list(fnscale=-1))$par
    }
  }
  #
  c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phiCS,nu)
}
fixpt.DEC <- function(thetav,y,x,z,time,ind,distr,lb,lu,luDEC,parallelphi,parallelnu,diagD){
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2;N<-length(y)
  m <- n_distinct(ind)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  D1 <- Dmatrix(thetav[(p+2):(p+1+q2)])
  phiDEC<-thetav[(p+2+q2)]
  thetaDEC<-thetav[(p+3+q2)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+3+q2))]
  #
  res_emj = revert_list(tapply(1:N,ind,emjDECs,y=y, x=x, z=z,time=time, beta1=beta1, D1=D1,
                               sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC, distr=distr,nu=nu))
  sum1 = Reduce("+",res_emj$sum1)
  sum2 = Reduce("+",res_emj$sum2)
  sum3 = sum(unlist(res_emj$sum3))
  sum4 = Reduce("+",res_emj$sum4)
  
  #if (calcbi) bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
  beta1<-solve(sum1)%*%sum2
  sigmae<-as.numeric(sum3)/N
  D1<-sum4/m
  if (diagD) D1 <-diag(diag(D1))
  #
  if (parallelphi) {
    parDEC<- optimParallel(c(phiDEC,thetaDEC),lcDEC,gr = NULL,method = "L-BFGS-B", lower =rep(0.0001,2),
                   upper = c(.9999,luDEC),control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  } else {
    parDEC<- optim(c(phiDEC,thetaDEC),lcDEC,gr = NULL,method = "L-BFGS-B", lower =rep(0.0001,2),
                   upper = c(.9999,luDEC),control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  }
  phiDEC<-parDEC[1]; thetaDEC<-parDEC[2]
  #
  logvero1<-function(nu){logveroDECs(y, x, z,time, ind, beta1, sigmae,phiDEC,thetaDEC, D1, distr, nu)}
  
  if (distr=="sn"){ nu<-NULL} else {
    if (parallelnu) {
      nu <- optimParallel(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, 
                  upper = lu,control = list(fnscale=-1))$par
    } else {
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, 
                  upper = lu,control = list(fnscale=-1))$par
    }
  }
  #
  c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phiDEC,thetaDEC,nu)
}
fixpt.CAR1 <- function(thetav,y,x,z,time,ind,distr,lb,lu,parallelphi,parallelnu,diagD){
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2;N<-length(y)
  m <- n_distinct(ind)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  D1 <- Dmatrix(thetav[(p+2):(p+1+q2)])
  phiDEC<-thetav[(p+2+q2)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+2+q2))]
  #
  res_emj = revert_list(tapply(1:N,ind,emjCAR1s,y=y, x=x, z=z,time=time, beta1=beta1, D1=D1,
                               sigmae=sigmae,phiDEC=phiDEC,distr=distr,nu=nu))
  sum1 = Reduce("+",res_emj$sum1)
  sum2 = Reduce("+",res_emj$sum2)
  sum3 = sum(unlist(res_emj$sum3))
  sum4 = Reduce("+",res_emj$sum4)
  
  #if (calcbi) bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
  
  beta1<-solve(sum1)%*%sum2
  sigmae<-as.numeric(sum3)/N
  D1<-sum4/m
  if (diagD) D1 <-diag(diag(D1))
  #
  if (parallelphi) {
    phiDEC<- optimParallel(phiDEC,lcCAR1,gr = NULL,method = "L-BFGS-B", lower =0.0001,
                   upper = .9999,control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  } else {
    phiDEC<- optim(phiDEC,lcCAR1,gr = NULL,method = "L-BFGS-B", lower =0.0001,
                   upper = .9999,control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
  }
  #
  logvero1<-function(nu){logveroCAR1s(y, x, z,time, ind, beta1, sigmae,phiDEC, D1, distr, nu)}
  
  if (distr=="sn"){ nu<-NULL} else {
    if (parallelnu) {
      nu <- optimParallel(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,
                  control = list(fnscale=-1))$par
    } else {
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,
                  control = list(fnscale=-1))$par
    }
  }
  #
  c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phiDEC,nu)
}
