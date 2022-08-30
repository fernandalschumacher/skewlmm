#SMSN
DAAREM.SkewAR<- function(formFixed,formRandom,data,groupVar,pAR,timeVar,
                        distr,beta1,sigmae,phiAR,D1,lambda,nu,lb,lu,
                        skewind,diagD,
                        precisao,informa,max.iter,showiter,showerroriter,
                        algorithm, #"EM" or "DAAREM"
                        parallelphi,parallelnu,ncores,#=detectCores()
                        control.daarem=list()){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  #varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]

  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if ((!is.null(phiAR)) && pAR!=length(phiAR)) stop("initial value from phi must be in agreement with pAR")
  if ((pAR%%1)!=0||pAR==0) stop("pAR must be integer greater than 1")

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)

  if (is.null(phiAR)) {
    lmeAR <- try(lme(formFixed,random=~1|ind,data=data,correlation=corARMA(p=pAR,q=0)),silent=T)
    if (class(lmeAR)=="try-error") piAR =as.numeric(pacf(y-x%*%beta1,lag.max=pAR,plot=F)$acf)
    else {
      phiAR <- capture.output(lmeAR$modelStruct$corStruct)[3]
      phiAR <- as.numeric(strsplit(phiAR, " ")[[1]])
      phiAR <- phiAR[!is.na(phiAR)]
      piAR <- tphitopi(phiAR)
    }
  } else piAR <- tphitopi(phiAR)
  if (any(piAR< -1 | piAR>1)) stop("invalid initial value from phi")

  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,piAR,nu)
  ##
  llji <- logveroARpi(y, x, z, time,ind, beta1, sigmae,piAR, D1, lambda, distr, nu)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,length(phiAR)*parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovARp","estphit","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovARp","estphit","logveroARpi","logveroAR","matrix.sqrt",
                          "dmvnorm","ljtAR","ljsAR","ljcnAR"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovARp","estphit","logveroARpi",
                          "logveroAR","matrix.sqrt","dmvnorm","ljtAR","ljsAR","ljcnAR"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.skewAR,objfn = objfn.skewAR,
                    y=y,x=x,z=z,time=time,ind=ind,distr=distr,pAR=pAR,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,skewind=skewind)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.skewAR,objfn = objfn.skewAR,
                    y=y,x=x,z=z,time=time,ind=ind,distr=distr,pAR=pAR,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,skewind=skewind)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  Gammab <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  Deltab<-EMout$par[(p+2+q2):(p+1+q2+q1)]
  piAR<-EMout$par[(p+2+q2+q1):(p+1+q2+q1+pAR)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-EMout$par[-(1:(p+1+q2+q1+pAR))]
  #
  D1<-Gammab+Deltab%*%t(Deltab)
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))*skewind
  zeta<-matrix.sqrt(sD1)%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjAR,y=y, x=x, z=z, time=time, beta1=beta1, Gammab=Gammab,
                             Deltab=Deltab, sigmae=sigmae,piAR=piAR,zeta=zeta, distr=distr,nu=nu,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calc_ui,y=y,time=time, x=x, z=z, beta1=beta1, Gammab=Gammab,
               Deltab=Deltab, sigmae=sigmae,phi=estphit(piAR),depStruct="ARp", distr=distr,nu=nu,simplify = TRUE)
  phiAR<-estphit(piAR)
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiAR,dd,lambda[skewind==1],nu)
  #theta <- c(beta1,sigmae,phiAR,dd,lambda,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",paste0("phiAR",1:length(piAR)),
                                   names_dd, paste0("lambda",(1:q1)[skewind==1]))
  else names(theta)<- c(colnames(x),"sigma2",paste0("phiAR",1:length(piAR)),
                        names_dd, paste0("lambda",(1:q1)[skewind==1]),
                        paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 phi=phiAR,dsqrt=dd,D=D1,lambda=as.numeric(lambda)),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixAR(y,x,z,time,ind,beta1,sigmae,phiAR,D1,lambda,
                             distr = distr,nu = nu,skewind = skewind,diagD=diagD),
                 silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      desvios[substr(names(theta),1,6) == 'lambda'] <- NA
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}

DAAREM.SkewUNC<- function(formFixed,formRandom,data,groupVar,
                         distr,beta1,sigmae,D1,lambda,nu,lb,lu,
                         skewind,diagD,
                         precisao,informa,max.iter,showiter,showerroriter,
                         algorithm, #"EM" or "DAAREM"
                         parallelnu,ncores,
                         control.daarem=list()){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  #
  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  #
  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,nu)
  ##
  llji <- logvero(y, x, z, ind, beta1, sigmae, D1, lambda, distr, nu)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu) {
    ncores <- min(ncores,1+2*length(nu))
    cl <- makeCluster(ncores) # set the number of processor cores
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    clusterExport(cl, c("logvero","matrix.sqrt","dmvnorm","ljt","ljs","ljcn",
                        "n_distinct"),
                    envir=environment())
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.skewUNC,objfn = objfn.skewUNC,
                    y=y,x=x,z=z,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelnu=parallelnu,
                    diagD=diagD,skewind=skewind)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.skewUNC,objfn = objfn.skewUNC,
                    y=y,x=x,z=z,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelnu=parallelnu,diagD=diagD,skewind=skewind)
  }

  if (parallelnu) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  Gammab <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  Deltab<-EMout$par[(p+2+q2):(p+1+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-EMout$par[-(1:(p+1+q2+q1))]
  #
  D1<-Gammab+Deltab%*%t(Deltab)
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))*skewind
  zeta<-matrix.sqrt(sD1)%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emj,y=y, x=x, z=z, beta1=beta1, Gammab=Gammab,
                             Deltab=Deltab, sigmae=sigmae, zeta=zeta, distr=distr,nu=nu,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calc_ui,y=y, x=x, z=z, beta1=beta1, Gammab=Gammab,Deltab=Deltab,
               sigmae=sigmae,depStruct="UNC", distr=distr,nu=nu,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,dd,lambda[skewind==1],nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",names_dd,
                                   paste0("lambda",(1:q1)[skewind==1]))
  else names(theta)<- c(colnames(x),"sigma2",names_dd,
                        paste0("lambda",(1:q1)[skewind==1]),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 dsqrt=dd,D=D1,lambda=as.numeric(lambda)),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(Infmatrix(y,x,z,ind,beta1,sigmae,D1,lambda,distr = distr,nu = nu,
                           skewind = skewind,diagD=diagD),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      desvios[substr(names(theta),1,6) == 'lambda'] <- NA
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}

DAAREM.SkewCS<- function(formFixed,formRandom,data,groupVar,
                         distr,beta1,sigmae,phiCS,D1,lambda,nu,lb,lu,
                         skewind,diagD,
                         precisao,informa,max.iter,showiter,showerroriter,
                         algorithm, #"EM" or "DAAREM"
                         parallelphi,parallelnu,ncores,
                         control.daarem=list()){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  #
  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  #
  if (!is.null(phiCS) && length(phiCS)!=1) stop ("initial value from phi must have length 1 or be NULL")
  if (!is.null(phiCS)) if (phiCS<=0 | phiCS>=1) stop("0<initialValue$phi<1 needed")
  #
  if (is.null(phiCS)) {
    phiCS <- abs(as.numeric(pacf(y-x%*%beta1,lag.max=1,plot=F)$acf))
  }
  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,phiCS,nu)
  ##
  llji <- logveroCS(y, x, z,ind, beta1, sigmae,phiCS, D1, lambda, distr, nu)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovCS","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovCS","logveroCS","matrix.sqrt",
                          "dmvnorm","ljtCS","ljsCS","ljcnCS"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovCS","logveroCS","matrix.sqrt",
                          "dmvnorm","ljtCS","ljsCS","ljcnCS"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.skewCS,objfn = objfn.skewCS,
                    y=y,x=x,z=z,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,skewind=skewind)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.skewCS,objfn = objfn.skewCS,
                    y=y,x=x,z=z,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,
                    diagD=diagD,skewind=skewind)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  Gammab <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  Deltab<-EMout$par[(p+2+q2):(p+1+q2+q1)]
  phiCS<-EMout$par[(p+2+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-EMout$par[-(1:(p+2+q2+q1))]
  #
  D1<-Gammab+Deltab%*%t(Deltab)
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))*skewind
  zeta<-matrix.sqrt(sD1)%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjCS,y=y, x=x, z=z, beta1=beta1, Gammab=Gammab,
                             Deltab=Deltab, sigmae=sigmae,phiCS=phiCS,zeta=zeta, distr=distr,
                             nu=nu,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calc_ui,y=y, x=x, z=z, beta1=beta1, Gammab=Gammab,Deltab=Deltab,
               sigmae=sigmae,depStruct="CS", phi=phiCS,distr=distr,nu=nu,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiCS,dd,lambda[skewind==1],nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phiCS",
                                   names_dd, paste0("lambda",(1:q1)[skewind==1]))
  else names(theta)<- c(colnames(x),"sigma2","phiCS",names_dd,
                        paste0("lambda",(1:q1)[skewind==1]),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae, phi=phiCS,
                                 dsqrt=dd,D=D1,lambda=as.numeric(lambda)),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixCS(y,x,z,ind,beta1,sigmae,phiCS,D1,lambda,distr = distr,
                             nu = nu,skewind = skewind,diagD=diagD),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      desvios[substr(names(theta),1,6) == 'lambda'] <- NA
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}

DAAREM.SkewDEC<- function(formFixed,formRandom,data,groupVar,timeVar,
                         distr,beta1,sigmae,parDEC,D1,lambda,nu,lb,lu,luDEC,
                         skewind,diagD,
                         precisao,informa,max.iter,showiter,showerroriter,
                         algorithm, #"EM" or "DAAREM"
                         parallelphi,parallelnu,ncores,
                         control.daarem=list()){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  #varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]

  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if (!is.null(parDEC)) {
    if (length(parDEC)!=2) stop ("initial value from phi should have length 2 or NULL")
    if (parDEC[1]<=0||parDEC[1]>=1) stop("invalid initial value from phi1")
    if (parDEC[2]<=0) stop("invalid initial value from phi2")
    if (parDEC[2]>= luDEC) stop("initial value from phi2 must be smaller than luDEC")
  }

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)

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
  #phiDEC=parDEC[1]
  #thetaDEC=parDEC[2]
  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,parDEC,nu)
  ##
  llji <- logveroDEC(y, x, z, time,ind, beta1, sigmae,parDEC[1],parDEC[2], D1, lambda, distr, nu)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,2*parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","logveroDEC","matrix.sqrt",
                          "dmvnorm","ljtDEC","ljsDEC","ljcnDEC"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovDEC","logveroDEC","matrix.sqrt",
                          "dmvnorm","ljtDEC","ljsDEC","ljcnDEC"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.skewDEC,objfn = objfn.skewDEC,
                    y=y,x=x,z=z,time=time,ind=ind,distr=distr,lb=lb,lu=lu,luDEC=luDEC,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,skewind=skewind)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.skewDEC,objfn = objfn.skewDEC,
                    y=y,x=x,z=z,time=time,ind=ind,distr=distr,lb=lb,lu=lu,luDEC=luDEC,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,
                    diagD=diagD,skewind=skewind)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  Gammab <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  Deltab<-EMout$par[(p+2+q2):(p+1+q2+q1)]
  phiDEC<-EMout$par[(p+2+q2+q1)]
  thetaDEC<-EMout$par[(p+3+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-EMout$par[-(1:(p+3+q2+q1))]
  #
  D1<-Gammab+Deltab%*%t(Deltab)
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))*skewind
  zeta<-matrix.sqrt(sD1)%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjDEC,y=y, x=x, z=z, time=time, beta1=beta1, Gammab=Gammab,
                             Deltab=Deltab, sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC,zeta=zeta, distr=distr,nu=nu,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calc_ui,y=y,time=time, x=x, z=z, beta1=beta1, Gammab=Gammab,
               Deltab=Deltab, sigmae=sigmae,phi=c(phiDEC,thetaDEC),depStruct="DEC", distr=distr,nu=nu,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiDEC,thetaDEC,dd,lambda[skewind==1],nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phi1DEC","phi2DEC",
                                   names_dd, paste0("lambda",(1:q1)[skewind==1]))
  else names(theta)<- c(colnames(x),"sigma2","phi1DEC","phi2DEC",
                        names_dd, paste0("lambda",(1:q1)[skewind==1]),
                        paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 phi=c(phiDEC,thetaDEC),dsqrt=dd,D=D1,lambda=as.numeric(lambda)),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixDEC(y,x,z,time,ind,beta1,sigmae,phiDEC,thetaDEC,D1,
                              lambda,distr = distr,nu = nu, skewind = skewind,
                              diagD=diagD),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      desvios[substr(names(theta),1,6) == 'lambda'] <- NA
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}

DAAREM.SkewCAR1 <- function(formFixed,formRandom,data,groupVar,timeVar,
                          distr,beta1,sigmae,phiCAR1,D1,lambda,nu,lb,lu,
                          skewind,diagD,
                          precisao,informa,max.iter,showiter,showerroriter,
                          algorithm, #"EM" or "DAAREM"
                          parallelphi,parallelnu,ncores,
                          control.daarem=list()){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  #varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]

  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if (!is.null(phiCAR1) && length(phiCAR1)!=1) stop("initial value from phi must have length 1 or be NULL")
  if (!is.null(phiCAR1)) if (phiCAR1>=1 || phiCAR1<=0) stop ("0<initialValue$phi<1 needed")

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)

  if (is.null(phiCAR1)) {
    lmeCAR = try(lme(formFixed,random=~1|ind,data=data,correlation=corCAR1(form = ~time)),silent=T)
    if (class(lmeCAR)=="try-error") phiDEC =abs(as.numeric(pacf(y-x%*%beta1,lag.max=1,plot=F)$acf))
    else {
      phiDEC = capture.output(lmeCAR$modelStruct$corStruct)[3]
      phiDEC = as.numeric(strsplit(phiDEC, " ")[[1]])
    }
  } else phiDEC <- phiCAR1

  teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,phiDEC,nu)
  ##
  llji <- logveroCAR1(y, x, z, time,ind, beta1, sigmae,phiDEC, D1, lambda, distr, nu)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","logveroCAR1","matrix.sqrt",
                          "dmvnorm","ljtCAR1","ljsCAR1","ljcnCAR1"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovDEC","logveroCAR1","matrix.sqrt",
                          "dmvnorm","ljtCAR1","ljsCAR1","ljcnCAR1"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.skewCAR1,objfn = objfn.skewCAR1,
                    y=y,x=x,z=z,time=time,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD,skewind=skewind)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.skewCAR1,objfn = objfn.skewCAR1,
                    y=y,x=x,z=z,time=time,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,
                    diagD=diagD,skewind=skewind)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  Gammab <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  Deltab<-EMout$par[(p+2+q2):(p+1+q2+q1)]
  phiDEC<-EMout$par[(p+2+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-EMout$par[-(1:(p+2+q2+q1))]
  #
  D1<-Gammab+Deltab%*%t(Deltab)
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  if ((t(Deltab)%*%sD1%*%Deltab)>=1) Deltab<-Deltab/as.numeric(sqrt(t(Deltab)%*%sD1%*%Deltab+1e-4))
  lambda<-matrix.sqrt(sD1)%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%sD1%*%Deltab))*skewind
  zeta<-matrix.sqrt(sD1)%*%lambda
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjDEC,y=y, x=x, z=z, time=time, beta1=beta1, Gammab=Gammab,
                             Deltab=Deltab, sigmae=sigmae,phiDEC=phiDEC,thetaDEC=1,zeta=zeta, distr=distr,nu=nu,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calc_ui,y=y,time=time, x=x, z=z, beta1=beta1, Gammab=Gammab,
               Deltab=Deltab, sigmae=sigmae,phi=c(phiDEC,1),depStruct="DEC", distr=distr,nu=nu,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiDEC,dd,lambda[skewind==1],nu)

  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phiCAR1",
                                   names_dd, paste0("lambda",(1:q1)[skewind==1]))
  else names(theta)<- c(colnames(x),"sigma2","phiCAR1",
                        names_dd, paste0("lambda",(1:q1)[skewind==1]),
                        paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 phi=phiDEC,dsqrt=dd,D=D1,lambda=as.numeric(lambda)),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixCAR1(y,x,z,time,ind,beta1,sigmae,phiDEC,D1,lambda,
                               distr = distr,nu = nu,skewind = skewind,
                               diagD=diagD),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      desvios[substr(names(theta),1,6) == 'lambda'] <- NA
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}

#SMN
DAAREM.AR<- function(formFixed,formRandom,data,groupVar,pAR,timeVar,
                     distr,beta1,sigmae,phiAR,D1,nu,lb,lu,diagD,
                     precisao,informa,max.iter,showiter,showerroriter,
                     algorithm, #"EM" or "DAAREM"
                     parallelphi,parallelnu,ncores,
                     control.daarem=list()){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  #varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]

  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if ((!is.null(phiAR)) && pAR!=length(phiAR)) stop("initial value from phi must be in agreement with pAR")
  if ((pAR%%1)!=0||pAR==0) stop("pAR must be integer greater than 1")

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

  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],piAR,nu)
  ##
  llji <- logveroARpis(y, x, z, time,ind, beta1, sigmae,piAR, D1, distr, nu)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,length(phiAR)*parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovARp","estphit","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovARp","estphit","logveroARpis","logveroARs",
                          "matrix.sqrt","dmvnorm","ljtARs","ljsARs","ljcnARs"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovARp","estphit","logveroARpis",
                          "logveroARs","matrix.sqrt","dmvnorm","ljtARs","ljsARs","ljcnARs"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.AR,objfn = objfn.AR,
                    y=y,x=x,z=z,time=time,ind=ind,distr=distr,pAR=pAR,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.AR,objfn = objfn.AR,
                    y=y,x=x,z=z,time=time,ind=ind,distr=distr,pAR=pAR,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,diagD=diagD)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  D1 <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  piAR<-EMout$par[(p+2+q2):(p+1+q2+pAR)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-EMout$par[-(1:(p+1+q2+pAR))]
  #
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjARs,y=y, x=x, z=z,time=time, beta1=beta1, D1=D1,
                             sigmae=sigmae,piAR=piAR, distr=distr,nu=nu,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calcs_ui,y=y,time=time, x=x, z=z, beta1=beta1, D1=D1,
               sigmae=sigmae,phi=estphit(piAR),depStruct="ARp", distr=distr,nu=nu,
               simplify = TRUE)
  #cat('\n',D1)
  phiAR<-estphit(piAR)
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiAR,dd,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",paste0("phiAR",1:length(piAR)),
                                   names_dd)
  else names(theta)<- c(colnames(x),"sigma2",paste0("phiAR",1:length(piAR)),
                        names_dd,paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 phi=phiAR,dsqrt=dd,D=D1),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixAR(y,x,z,time,ind,beta1,sigmae,phiAR,D1,lambda=rep(0,q1),
                             distr = distr,nu = nu,diagD=diagD,skewind=rep(0,q1)),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}

DAAREM.UNC<- function(formFixed,formRandom,data,groupVar,
                      distr,beta1,sigmae,D1,nu,lb,lu,diagD,
                      precisao,informa,max.iter,showiter,showerroriter,
                      algorithm, #"EM" or "DAAREM"
                      parallelnu,ncores,
                      control.daarem=list()){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  #
  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],nu)
  ##
  llji <- logveros(y, x, z, ind, beta1, sigmae, D1, distr, nu)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu) {
    ncores <- min(ncores,1+2*length(nu))
    cl <- makeCluster(ncores) # set the number of processor cores
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    clusterExport(cl, c("logveros","matrix.sqrt","dmvnorm","ljts","ljss","ljcns",
                        "n_distinct"),
                  envir=environment())
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.UNC,objfn = objfn.UNC,
                    y=y,x=x,z=z,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelnu=parallelnu,diagD=diagD)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.UNC,objfn = objfn.UNC,
                    y=y,x=x,z=z,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,
                    showerroriter = showerroriter,parallelnu=parallelnu,diagD=diagD)
  }

  if (parallelnu) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  D1 <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  if (distr=="sn") {
    nu <- NULL
  } else nu<-EMout$par[-(1:(p+1+q2))]
  #
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjs,y=y, x=x, z=z, beta1=beta1, D1=D1,
                             sigmae=sigmae, distr=distr,nu=nu,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calcs_ui,y=y, x=x, z=z, beta1=beta1, D1=D1,
               sigmae=sigmae,depStruct="UNC", distr=distr,nu=nu,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,dd,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",names_dd)
  else names(theta)<- c(colnames(x),"sigma2",names_dd,paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 dsqrt=dd,D=D1),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(Infmatrix(y,x,z,ind,beta1,sigmae,D1,lambda=rep(0,q1),distr = distr,
                           nu = nu,diagD=diagD,skewind=rep(0,q1)),
                 silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}

DAAREM.CS<- function(formFixed,formRandom,data,groupVar,
                     distr,beta1,sigmae,phiCS,D1,nu,lb,lu,diagD,
                     precisao,informa,max.iter,showiter,showerroriter,
                     algorithm, #"EM" or "DAAREM"
                     parallelphi,parallelnu,ncores,
                     control.daarem=list()){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  #
  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if (!is.null(phiCS) && length(phiCS)!=1) stop ("initial value from phi must have length 1 or be NULL")
  if (!is.null(phiCS)) if (phiCS<=0 | phiCS>=1) stop("0<initialValue$phi<1 needed")
  #
  if (is.null(phiCS)) {
    phiCS <- abs(as.numeric(pacf(y-x%*%beta1,lag.max=1,plot=F)$acf))
  }
  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phiCS,nu)
  ##
  llji <- logveroCSs(y, x, z,ind, beta1, sigmae,phiCS, D1, distr, nu)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovCS","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovCS","logveroCSs","matrix.sqrt",
                          "dmvnorm","ljtCSs","ljsCSs","ljcnCSs"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovCS","logveroCSs","matrix.sqrt",
                          "dmvnorm","ljtCSs","ljsCSs","ljcnCSs"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.CS,objfn = objfn.CS,
                    y=y,x=x,z=z,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.CS,objfn = objfn.CS,
                    y=y,x=x,z=z,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,diagD=diagD)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  D1 <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  phiCS<-EMout$par[(p+2+q2)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-EMout$par[-(1:(p+2+q2))]
  #
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjCSs,y=y, x=x, z=z, beta1=beta1, D1=D1,
                             sigmae=sigmae,phiCS=phiCS, distr=distr,nu=nu,simplify = FALSE)),ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calcs_ui,y=y, x=x, z=z, beta1=beta1, D1=D1,
               sigmae=sigmae,depStruct="CS", phi=phiCS,distr=distr,nu=nu,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiCS,dd,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phiCS",names_dd)
  else names(theta)<- c(colnames(x),"sigma2","phiCS",names_dd,
                        paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae, phi=phiCS,
                                 dsqrt=dd,D=D1),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixCS(y,x,z,ind,beta1,sigmae,phiCS,D1,lambda=rep(0,q1),
                             distr = distr,nu = nu,diagD=diagD,skewind=rep(0,q1)),
                             silent = T)
    #desvios<-try(InfmatrixCS(y,x,z,ind,beta1,sigmae,phiCS,D1,lambda,distr = distr,nu = nu),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}

DAAREM.DEC<- function(formFixed,formRandom,data,groupVar,timeVar,
                      distr,beta1,sigmae,parDEC,D1,nu,lb,lu,luDEC,diagD,
                      precisao,informa,max.iter,showiter,showerroriter,
                      algorithm, #"EM" or "DAAREM"
                      parallelphi,parallelnu,ncores,
                      control.daarem=list()){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  #varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]

  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if (!is.null(parDEC)) {
    if (length(parDEC)!=2) stop ("initial value from phi should have length 2 or NULL")
    if (parDEC[1]<=0||parDEC[1]>=1) stop("invalid initial value from phi1")
    if (parDEC[2]<=0) stop("invalid initial value from phi2")
    if (parDEC[2]>= luDEC) stop("initial value from phi2 must be smaller than luDEC")
  }

  if (is.null(parDEC)) {
    #cat("calculating initial values for DEC... \n")
    thetat<- seq(0.1,2,by=.1)
    phit <- seq(0.1,.9,by=.05)
    vect <-merge(phit,thetat,all=T)
    logveroDECv<-function(phitheta){logveroDECs(y, x, z,time,ind, beta1=beta1, sigmae=sigmae,
                                               phiDEC=phitheta[1],thetaDEC=phitheta[2],
                                               D1=D1,distr=distr, nu=nu)}
    logverovec <- apply(vect,1,logveroDECv)
    parDEC <- as.numeric(vect[which.max(logverovec),])
  }
  #phiDEC=parDEC[1]
  #thetaDEC=parDEC[2]
  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],parDEC,nu)
  ##
  llji <- logveroDECs(y, x, z, time,ind, beta1, sigmae,parDEC[1],parDEC[2], D1, distr, nu)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,2*parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","logveroDECs","matrix.sqrt",
                          "dmvnorm","ljtDECs","ljsDECs","ljcnDECs"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovDEC","logveroDECs","matrix.sqrt",
                          "dmvnorm","ljtDECs","ljsDECs","ljcnDECs"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.DEC,objfn = objfn.DEC,
                    y=y,x=x,z=z,time=time,ind=ind,distr=distr,lb=lb,lu=lu,luDEC=luDEC,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.DEC,objfn = objfn.DEC,
                    y=y,x=x,z=z,time=time,ind=ind,distr=distr,lb=lb,lu=lu,luDEC=luDEC,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,diagD=diagD)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  D1 <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  phiDEC<-EMout$par[(p+2+q2)]
  thetaDEC<-EMout$par[(p+3+q2)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-EMout$par[-(1:(p+3+q2))]
  #
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjDECs,y=y, x=x, z=z,time=time, beta1=beta1, D1=D1,
                             sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC, distr=distr,nu=nu,simplify = FALSE)),
               ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calcs_ui,y=y,time=time, x=x, z=z, beta1=beta1, D1=D1,
               sigmae=sigmae,phi=c(phiDEC,thetaDEC),depStruct="DEC", distr=distr,nu=nu,simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiDEC,thetaDEC,dd,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phi1DEC","phi2DEC",
                                   names_dd)
  else names(theta)<- c(colnames(x),"sigma2","phi1DEC","phi2DEC",
                        names_dd,paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 phi=c(phiDEC,thetaDEC),dsqrt=dd,D=D1),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixDEC(y,x,z,time,ind,beta1,sigmae,phiDEC,thetaDEC,D1,
                              lambda=rep(0,q1),distr = distr,nu = nu,diagD=diagD,
                              skewind=rep(0,q1)),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}

DAAREM.CAR1 <- function(formFixed,formRandom,data,groupVar,timeVar,
                        distr,beta1,sigmae,phiCAR1,D1,nu,lb,lu,diagD,
                        precisao,informa,max.iter,showiter,showerroriter,
                        algorithm, #"EM" or "DAAREM"
                        parallelphi,parallelnu,ncores,
                        control.daarem=list()){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  #varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind
  if (is.null(timeVar)) {
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  } else time <- data[,timeVar]

  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x);  q1<-ncol(z); q2 <- q1*(q1+1)/2
  #
  if (!is.null(phiCAR1) && length(phiCAR1)!=1) stop("initial value from phi must have length 1 or be NULL")
  if (!is.null(phiCAR1)) if (phiCAR1>=1 || phiCAR1<=0) stop ("0<initialValue$phi<1 needed")

  if (is.null(phiCAR1)) {
    lmeCAR = try(lme(formFixed,random=~1|ind,data=data,correlation=corCAR1(form = ~time)),silent=T)
    if (class(lmeCAR)=="try-error") phiDEC =abs(as.numeric(pacf(y-x%*%beta1,lag.max=1,plot=F)$acf))
    else {
      phiDEC = capture.output(lmeCAR$modelStruct$corStruct)[3]
      phiDEC = as.numeric(strsplit(phiDEC, " ")[[1]])
    }
  } else phiDEC <- phiCAR1

  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phiDEC,nu)
  ##
  llji <- logveroCAR1s(y, x, z, time,ind, beta1, sigmae,phiDEC, D1, distr, nu)
  if (is.nan(llji)||is.infinite(abs(llji))) stop("NaN/infinity initial likelihood, please change initial parameter values")

  if (parallelnu||parallelphi) {
    ncores <- min(ncores,1+2*max(length(nu)*parallelnu,parallelphi))
    cl <- makeCluster(ncores) # set the number of processor cores
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    if (parallelphi && !parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","traceM"),
                    envir=environment())
    } else if (!parallelphi && parallelnu) {
      clusterExport(cl, c("n_distinct","CovDEC","logveroCAR1s","matrix.sqrt",
                          "dmvnorm","ljtCAR1s","ljsCAR1s","ljcnCAR1s"),
                    envir=environment())
    } else {
      clusterExport(cl, c("n_distinct","traceM","CovDEC","logveroCAR1s","matrix.sqrt",
                          "dmvnorm","ljtCAR1s","ljsCAR1s","ljcnCAR1s"),
                    envir=environment())
    }
  }

  if (algorithm=="EM") {
    EMout <- fpiter(par=teta,fixptfn = fixpt.CAR1,objfn = objfn.CAR1,
                    y=y,x=x,z=z,time=time,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=list(tol=precisao,maxiter=max.iter),showiter = showiter,
                    showerroriter = showerroriter,parallelphi=parallelphi,
                    parallelnu=parallelnu,diagD=diagD)
  } else {
    control.daarem$tol = precisao
    control.daarem$maxiter = max.iter
    EMout <- daarem(par=teta,fixptfn = fixpt.CAR1,objfn = objfn.CAR1,
                    y=y,x=x,z=z,time=time,ind=ind,distr=distr,lb=lb,lu=lu,
                    control=control.daarem,showiter = showiter,showerroriter = showerroriter,
                    parallelphi=parallelphi,parallelnu=parallelnu,diagD=diagD)
  }

  if (parallelnu||parallelphi) stopCluster(cl)
  if (!EMout$convergence) message("maximum number of iterations reachead \n")
  #
  beta1<-matrix(EMout$par[1:p],ncol=1)
  sigmae<-as.numeric(EMout$par[p+1])
  D1 <- Dmatrix(EMout$par[(p+2):(p+1+q2)])
  phiDEC<-EMout$par[(p+2+q2)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-EMout$par[-(1:(p+2+q2))]
  #
  if (diagD) D1<-diag(diag(D1))
  sD1 <- solve(D1)
  bi <- matrix(unlist(tapply(1:N,ind,calcbi_emjDECs,y=y, x=x, z=z,time=time, beta1=beta1, D1=D1,
                             sigmae=sigmae,phiDEC=phiDEC,thetaDEC=1, distr=distr,nu=nu,simplify = FALSE)),
               ncol=q1,byrow = T)
  ui <- tapply(1:N,ind,calcs_ui,y=y,time=time, x=x, z=z, beta1=beta1, D1=D1,
               sigmae=sigmae,phi=c(phiDEC,1),depStruct="DEC", distr=distr,nu=nu,
               simplify = TRUE)
  #
  if (diagD) {
    dd<-diag(matrix.sqrt(D1))
    names_dd <- paste0("Dsqrt",1:q1,1:q1)
  } else{
    dd<-try(matrix.sqrt(D1)[upper.tri(D1, diag = T)],silent = T)
    if (class(dd)[1]=='try-error') stop("Numerical error, try using algorithm = 'EM' in control")
    names_dd <- matrix(paste0("Dsqrt",rep(1:q1,q1),rep(1:q1,each=q1)),ncol=q1)[upper.tri(D1, diag = T)]
  }
  theta <- c(beta1,sigmae,phiDEC,dd,nu)

  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phiCAR1",names_dd)
  else names(theta)<- c(colnames(x),"sigma2","phiCAR1",names_dd,
                        paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = EMout$fpevals,
                  estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                 phi=phiDEC,dsqrt=dd,D=D1),
                  uhat=ui,loglik.track=EMout$objfn.track) ###

  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixCAR1(y,x,z,time,ind,beta1,sigmae,phiDEC,D1,lambda=rep(0,q1),
                               distr = distr,nu = nu,diagD=diagD,skewind=rep(0,q1)),
                 silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  obj.out$loglik <-as.numeric(EMout$value.objfn)

  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=EMout$criterio
  obj.out
}
