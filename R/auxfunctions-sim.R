################################################################
#Log-likelihood - AR(p)
################################################################
ljnormalARs <-function(j,y,x,z,time,beta1,D1,sigmae,phiAR){
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovARp(phiAR,t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log(dmvnorm(y1,med,Psi))
}
#
ljtARs <-function(j,nu,y,x,z,time,beta1,D1,sigmae,phiAR){
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovARp(phiAR,t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  #dmvt(y1,delta = med,sigma = Psi,df=nu,log=F)
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  log(dtj)
}
#
ljsARs <-function(j,nu,y,x,z,time,beta1,D1,sigmae,phiAR){
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovARp(phiAR,t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))
  resp <- integrate(Vectorize(f2),0,1)$value
  log(nu*resp)
}
#uhat <- pgamma(1,njj/2+nu+1,dj/2)/pgamma(1,njj/2+nu,dj/2)*(njj+2*nu)/dj
#fteste <- function(u) u*dgamma(u,njj/2+nu,dj/2)/pgamma(1,njj/2+nu,dj/2)
#integrate(fteste,0,1)
#
ljcnARs <-function(j,nu,y,x,z,time,beta1,D1,sigmae,phiAR){
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovARp(phiAR,t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
        (1-nu[1])*dmvnorm(y1,med,Psi))
}

logveroARs = function(y,x,z,time,ind,beta1,sigmae,phiAR,D1,distr,nu){ #ind = indicadora de individuo
  N <-length(ind)

  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalARs,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtARs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsARs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnARs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiAR=phiAR))
  lv
}
logveroARpis = function(y,x,z,time,ind,beta1,sigmae,piAR,D1,distr,nu){ #ind = indicadora de individuo

  phiAR <- estphit(piAR)
  N <-length(ind)

  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalARs,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtARs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsARs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnARs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiAR=phiAR))
  lv
}

##############################################################################
# EM - AR(p)
##############################################################################
calcbi_emjARs <- function(jseq,y,x,z,time,beta1,D1,sigmae,piAR,distr,nu) {
  y1=y[jseq]
  t1=time[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1
  nj = length(y1)
  Sigma = sigmae*CovARp(phi = estphit(piAR),t1)
  Psi<-(z1)%*%D1%*%t(z1)+Sigma
  #
  bi<-D1%*%t(z1)%*%solve(Psi)%*%(y1-med)
  bi
}
emjARs = function(jseq, y, x, z,time, beta1, D1,sigmae,piAR,distr,nu) {
  y1=y[jseq]
  t1=time[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1
  nj = length(y1)
  Sigma = sigmae*CovARp(phi = estphit(piAR),t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  #
  if  (distr=="sn"){
    uj<-1
  }
  if (distr=="st"){
    uj<-(nj+nu)/(dj+nu)
  }

  if (distr=="ss"){
    uj<-pgamma(1,nj/2+nu+1,dj/2)/pgamma(1,nj/2+nu,dj/2)*(nj+2*nu)/dj
  }

  if (distr=="scn"){
    fy<-as.numeric((nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))
    uj<-as.numeric((nu[1]*nu[2]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))/fy
  }

  sSigma = solve(Sigma)
  sRi = sSigma*sigmae
  Tbj<-solve(solve(D1)+t(z1)%*%sSigma%*%z1)
  r<-Tbj%*%t(z1)%*%sSigma%*%(y1-x1%*%beta1)
  ub<-uj*r
  ub2j<-Tbj+uj*r%*%t(r)
  #
  sum1<-uj*t(x1)%*%sRi%*%x1 #denom beta
  sum2<-(t(x1)%*%sRi%*%(uj*y1-z1%*%ub)) #num beta
  sum3<-uj*t(y1-x1%*%beta1)%*%sRi%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%sRi%*%z1%*%ub-
    t(ub)%*%t(z1)%*%sRi%*%(y1-x1%*%beta1)+traceM(sRi%*%z1%*%ub2j%*%t(z1)) #soma do sig2
  sum4<-ub2j #soma do D1
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,uj=uj,ubj=ub,ub2j=ub2j)
  #if (calcbi) obj.out$bi=bi
  return(obj.out)
}

EM.AR<- function(formFixed,formRandom,data,groupVar,pAR,timeVar,
                 distr,beta1,sigmae,phiAR,D1,nu,lb,lu,
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
  #
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

  criterio<-10
  count<-0
  llji = logveroARpis(y, x, z, time,ind, beta1, sigmae,piAR, D1, distr, nu)
  if (is.nan(llji)|is.infinite(abs(llji))) stop("NaN/infinity initial likelihood")

  while((criterio > precisao)&(count<max.iter)){
    #print(nu)

    count <- count + 1
    res_emj = revert_list(tapply(1:N,ind,emjARs,y=y, x=x, z=z,time=time, beta1=beta1, D1=D1,
                                 sigmae=sigmae,piAR=piAR, distr=distr,nu=nu))
    sum1 = Reduce("+",res_emj$sum1)
    sum2 = Reduce("+",res_emj$sum2)
    sum3 = sum(unlist(res_emj$sum3))
    sum4 = Reduce("+",res_emj$sum4)
    
    #if (calcbi) bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
    beta1<-solve(sum1)%*%sum2
    sigmae<-as.numeric(sum3)/N
    D1<-sum4/m
    #
    piAR<- optim(piAR,lcAR,gr = NULL,method = "L-BFGS-B", lower =rep(-.9999,pAR),
                 upper = rep(.9999,pAR),control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                 y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
    #
    logvero1<-function(nu){logveroARpis(y = y,x = x, z = z,time = time,ind =  ind,beta1 =  beta1,sigmae = sigmae,piAR = piAR,D1 =  D1, distr, nu)}

    if (distr=="sn"){ nu<-NULL} else
    {
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par
    }
    #
    param <- teta
    teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],piAR,nu)
    criterio2 <- as.numeric(sqrt((teta-param)%*%(teta-param)))
    llj1<-llji
    llji <- logveroARpis(y, x, z, time,ind, beta1, sigmae,piAR, D1, distr, nu)
    criterio <- abs((llji-llj1)/llj1)
    if (showiter&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (count==max.iter) message("\n maximum number of iterations reachead")
  }

  cat("\n")
  bi <- t(bind_cols(tapply(1:N,ind,calcbi_emjARs,y=y, x=x, z=z,time=time, beta1=beta1, D1=D1,
                           sigmae=sigmae,piAR=piAR, distr=distr,nu=nu,simplify = FALSE)))
  dd<-matrix.sqrt(D1)[upper.tri(D1, diag = T)]
  phiAR=estphit(piAR)
  theta = c(beta1,sigmae,phiAR,dd,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",paste0("phiAR",1:length(piAR)),paste0("Dsqrt",1:length(dd)))
  else names(theta)<- c(colnames(x),"sigma2",paste0("phiAR",1:length(piAR)),paste0("Dsqrt",1:length(dd)),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                                           phi=phiAR,dsqrt=dd),
                  uhat=unlist(res_emj$uj))
  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixARs(y,x,z,time,ind,beta1,sigmae,phiAR,D1,distr = distr,nu = nu),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
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
####
predictf.AR<- function(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr,pAR,theta){
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
  if (distr=="sn") {nu<- NULL}
  if (distr=="st") {nu<- theta[p+pAR+q2+2]}
  if (distr=="ss") {nu<- theta[p+pAR+q2+2]}
  if (distr=="scn") {nu<- theta[(p+pAR+q2+2):(p+pAR+q2+3)]}
  if ((p+pAR+1+q2+length(nu))!=length(theta)) stop("theta misspecified")
  D1sqrt <- Dmatrix(dd)
  D1 <- D1sqrt%*%D1sqrt
  #
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
    medFit <- x1%*%beta1
    medPred <- x1Pred%*%beta1
    #
    y1=dataFitj[,all.vars(formFixed)[1]]
    SigmaPlus = sigmae*CovARp(phi = phiAR,c(dataFitj$time,dataPredj$time))
    PsiPlus<-(zPlus1)%*%(D1)%*%t(zPlus1)+SigmaPlus
    dj<-as.numeric(t(y1-medFit)%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit))
    LambdaPlus <- solve(solve(D1)+ t(zPlus1)%*%solve(SigmaPlus)%*%zPlus1)
    sPsiPlus <- solve(PsiPlus)
    Psi22.1 <- PsiPlus[seqPred,seqPred]- PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]
    ypredj <- medPred + PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit)
    ypred[dataPred$ind==indj] <- ypredj
    xpred[dataPred$ind==indj,] <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
    timepred[dataPred$ind==indj] <- dataPredj$time
  }
  colnames(xpred) = colnames(xPlus1)
  if (all(xpred[,1]==1)) xpred=xpred[,-1]
  data.frame(groupVar=dataPred$ind,time=timepred,xpred,ypred)
}

################################################################
#Log-likelihood - independent
################################################################
ljnormals <-function(j,y,x,z,beta1,D1,sigmae){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log(dmvnorm(y1,med,Psi))
}

ljts <-function(j,nu,y,x,z,beta1,D1,sigmae){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  log(dtj)
}

ljss <-function(j,nu,y,x,z,beta1,D1,sigmae){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))
  resp <- integrate(Vectorize(f2),0,1)$value
  log(nu*resp)
}

ljcns <-function(j,nu,y,x,z,beta1,D1,sigmae){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log((nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
         (1-nu[1])*dmvnorm(y1,med,Psi)))
}
logveros = function(y,x,z,ind,beta1,sigmae,D1,distr,nu){ #ind = indicadora de individuo
  m<-n_distinct(ind)
  N<-length(ind)
  p<-dim(x)[2]
  q1<-dim(z)[2]

  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormals,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljts,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljss,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcns,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  lv
}

# ##############################################################################
# # EM - independent
# ##############################################################################
calcbi_emjs = function(jseq, y, x, z, beta1,D1, sigmae,distr,nu) {
  y1=y[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1
  nj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(nj)
  #
  bi<-D1%*%t(z1)%*%solve(Psi)%*%(y1-med)
  bi
 }
emjs = function(jseq, y, x, z, beta1,D1, sigmae,distr,nu) {
  y1=y[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1
  nj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(nj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  #
  if  (distr=="sn"){
    uj<-1
  }
  if (distr=="st"){
    uj<-(nj+nu)/(dj+nu)
  }

  if (distr=="ss"){
    uj<-pgamma(1,nj/2+nu+1,dj/2)/pgamma(1,nj/2+nu,dj/2)*(nj+2*nu)/dj
  }

  if (distr=="scn"){
    fy<-as.numeric((nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))
    uj<-as.numeric((nu[1]*nu[2]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))/fy
  }

  Tbj<-solve(solve(D1)+t(z1)%*%z1/sigmae)
  r<-Tbj%*%t(z1)%*%(y1-x1%*%beta1)/sigmae
  ub<-uj*r
  ub2j<-Tbj+uj*r%*%t(r)
  #
  sum1<-uj*t(x1)%*%x1 #denom beta
  sum2<-(t(x1)%*%(uj*y1-z1%*%ub)) #num beta
  sum3<-uj*t(y1-x1%*%beta1)%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%z1%*%ub-
    t(ub)%*%t(z1)%*%(y1-x1%*%beta1)+traceM(ub2j%*%t(z1)%*%z1) #soma do sig2
  sum4<-ub2j #soma do Gamma
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,uj=uj)
  #if (calcbi) obj.out$bi=bi
  return(obj.out)
}
#
EM.sim<- function(formFixed,formRandom,data,groupVar,distr,beta1,sigmae,D1,nu,lb,lu,
                  precisao,informa,max.iter,showiter,showerroriter){
  ti = Sys.time()
  x <- model.matrix(formFixed,data=data)
  varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind

  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x)
  q1<-ncol(z)
  #
  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],nu)

  criterio<-10
  count<-0
  llji = logveros(y, x, z, ind, beta1, sigmae, D1, distr, nu)
  if (is.nan(llji)|is.infinite(abs(llji))) stop("NaN/infinity initial likelihood")

  while((criterio > precisao)&(count<max.iter)){

    count <- count + 1
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
    #
    logvero1<-function(nu){logveros(y, x, z, ind, beta1, sigmae, D1, distr, nu)}

    if (distr=="sn"){ nu<-NULL} else
    {
      nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par
    }
    param <- teta
    teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],nu)
    criterio2 <- as.numeric(sqrt((teta-param)%*%(teta-param)))
    llj1 <- llji
    llji <- logveros(y, x, z, ind, beta1, sigmae, D1, distr, nu)
    criterio <- abs((llji-llj1)/llj1)
    if (showiter&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (count==max.iter) message("\n maximum number of iterations reachead")
  }
  cat("\n")
  bi <- t(bind_cols(tapply(1:N,ind,calcbi_emjs,y=y, x=x, z=z, beta1=beta1, D1=D1,
                           sigmae=sigmae, distr=distr,nu=nu,simplify = FALSE)))
  dd<-matrix.sqrt(D1)[upper.tri(D1, diag = T)]
  theta = c(beta1,sigmae,dd,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2",paste0("Dsqrt",1:length(dd)))
  else names(theta)<- c(colnames(x),"sigma2",paste0("Dsqrt",1:length(dd)),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                                           dsqrt=dd),
                  uhat=unlist(res_emj$uj))
  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi
  
  if (informa) {
    desvios<-try(Infmatrixs(y,x,z,ind,beta1,sigmae,D1,distr = distr,nu = nu),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
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
#
predictf.sim<- function(formFixed,formRandom,dataFit,dataPred,groupVar,distr,theta){
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
  if (distr=="sn") {nu<- NULL}
  if (distr=="st") {nu<- theta[p+q2+2]}
  if (distr=="ss") {nu<- theta[p+q2+2]}
  if (distr=="scn") {nu<- theta[(p+q2+2):(p+q2+3)]}
  if ((p+1+q2+length(nu))!=length(theta)) stop("theta misspecified")
  D1sqrt <- Dmatrix(dd)
  D1 <- D1sqrt%*%D1sqrt
  #
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
    medFit <- x1%*%beta1
    medPred <- x1Pred%*%beta1
    #
    y1=dataFitj[,all.vars(formFixed)[1]]
    SigmaPlus = sigmae*diag(njFit+njPred)
    PsiPlus<-(zPlus1)%*%(D1)%*%t(zPlus1)+SigmaPlus
    dj<-as.numeric(t(y1-medFit)%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit))
    LambdaPlus <- solve(solve(D1)+ t(zPlus1)%*%solve(SigmaPlus)%*%zPlus1)
    sPsiPlus <- solve(PsiPlus)
    Psi22.1 <- PsiPlus[seqPred,seqPred]- PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]
    ypredj <- medPred + PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit)
    ypred[dataPred$ind==indj] <- ypredj
    xpred[dataPred$ind==indj,] <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
  }
  colnames(xpred) = colnames(xPlus1)
  if (all(xpred[,1]==1)) xpred=xpred[,-1]
  data.frame(groupVar=dataPred$ind,xpred,ypred)
}
#
# ################################################################
# #Log-likelihood - CS
# ################################################################
ljnormalCSs <-function(j,y,x,z,beta1,D1,sigmae,phiCS){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovCS(phiCS,njj)#CovARp(phiAR,t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log(dmvnorm(y1,med,Psi))
}
#
ljtCSs <-function(j,nu,y,x,z,beta1,D1,sigmae,phiCS){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovCS(phiCS,njj)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  log(dtj)
}
#
ljsCSs <-function(j,nu,y,x,z,beta1,D1,sigmae,phiCS){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovCS(phiCS,njj)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))
  resp <- integrate(Vectorize(f2),0,1)$value
  log(nu*resp)
}
#
ljcnCSs <-function(j,nu,y,x,z,beta1,D1,sigmae,phiCS){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovCS(phiCS,njj)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log((nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
         (1-nu[1])*dmvnorm(y1,med,Psi)))
}

logveroCSs = function(y,x,z,ind,beta1,sigmae,phiCS,D1,distr,nu){ #ind = indicadora de individuo

  N <-length(ind)

  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalCSs,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae,phiCS=phiCS))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtCSs,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae,phiCS=phiCS))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsCSs,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae,phiCS=phiCS))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnCSs,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae,phiCS=phiCS))
  lv
}

# ##############################################################################
# # EM - CS
# ##############################################################################
calcbi_emjCSs = function(jseq, y, x, z, beta1, D1, sigmae,phiCS,distr,nu) {
  y1=y[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1
  nj = length(y1)
  Sigma = sigmae*CovCS(phiCS,nj)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  #
  bi<-D1%*%t(z1)%*%solve(Psi)%*%(y1-med)
  bi
}
emjCSs = function(jseq, y, x, z, beta1, D1, sigmae,phiCS,distr,nu) {
  y1=y[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1
  nj = length(y1)
  Sigma = sigmae*CovCS(phiCS,nj)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  #
  if  (distr=="sn"){
    uj<-1
  }
  if (distr=="st"){
    uj<-(nj+nu)/(dj+nu)
  }

  if (distr=="ss"){
    uj<-pgamma(1,nj/2+nu+1,dj/2)/pgamma(1,nj/2+nu,dj/2)*(nj+2*nu)/dj
  }

  if (distr=="scn"){
    fy<-as.numeric((nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))
    uj<-as.numeric((nu[1]*nu[2]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))/fy
  }

  sSigma = solve(Sigma)
  sRi = sSigma*sigmae
  Tbj<-solve(solve(D1)+t(z1)%*%sSigma%*%z1)
  r<-Tbj%*%t(z1)%*%sSigma%*%(y1-x1%*%beta1)
  ub<-uj*r
  ub2j<-Tbj+uj*r%*%t(r)
  #
  sum1<-uj*t(x1)%*%sRi%*%x1 #denom beta
  sum2<-(t(x1)%*%sRi%*%(uj*y1-z1%*%ub)) #num beta
  sum3<-uj*t(y1-x1%*%beta1)%*%sRi%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%sRi%*%z1%*%ub-
    t(ub)%*%t(z1)%*%sRi%*%(y1-x1%*%beta1)+traceM(sRi%*%z1%*%ub2j%*%t(z1)) #soma do sig2
  sum4<-ub2j #soma do D1
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,uj=uj,ubj=ub,ub2j=ub2j)
  #if (calcbi) obj.out$bi=bi
  return(obj.out)
}
#
EM.CS<- function(formFixed,formRandom,data,groupVar,
                 distr,beta1,sigmae,phiCS,D1,nu,lb,lu,
                 precisao,informa,max.iter,showiter,showerroriter){
  ti <- Sys.time()
  x <- model.matrix(formFixed,data=data)
  varsx <- all.vars(formFixed)[-1]
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <- ind

  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x)
  q1<-ncol(z)
  #
  if (!is.null(phiCS) & length(phiCS!=1)) stop ("initial value from phi must have length 1 or be NULL")
  if (!is.null(phiCS)) if (phiCS<=0 | phiCS>=1) stop("0<initialValue$phi<1 needed")
  #
  if (is.null(phiCS)) {
    phiCS = abs(as.numeric(pacf(y-x%*%beta1,lag.max=1,plot=F)$acf))
  }

  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phiCS,nu)

  criterio<-10
  count<-0
  llji = logveroCSs(y, x, z,ind, beta1, sigmae,phiCS, D1, distr, nu)
  if (is.nan(llji)|is.infinite(abs(llji))) stop("NaN/infinity initial likelihood")

  while((criterio > precisao)&(count<max.iter)){
    #print(nu)

    count <- count + 1
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
    #
    phiCS <- optim(phiCS,lcCS,gr = NULL,method = "L-BFGS-B", lower =0,
                   upper = .9999,control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
    #
    logvero1<-function(nu){logveroCSs(y, x, z, ind, beta1, sigmae,phiCS, D1, distr, nu)}

    if (distr=="sn"){ nu<-NULL} else
    {nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par}
    #
    param <- teta
    teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phiCS,nu)
    criterio2 <- as.numeric(sqrt((teta-param)%*%(teta-param)))
    llj1<-llji
    llji <- logveroCSs(y, x, z, ind, beta1, sigmae,phiCS, D1, distr, nu)
    criterio <- abs((llji-llj1)/llj1)
    if (showiter&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (count==max.iter) message("\n maximum number of iterations reachead")
  }

  cat("\n")
  bi <- t(bind_cols(tapply(1:N,ind,calcbi_emjCSs,y=y, x=x, z=z, beta1=beta1, D1=D1,
                           sigmae=sigmae,phiCS=phiCS, distr=distr,nu=nu,simplify = FALSE)))
  
  dd<-matrix.sqrt(D1)[upper.tri(D1, diag = T)]
  theta = c(beta1,sigmae,phiCS,dd,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phiCS",paste0("Dsqrt",1:length(dd)))
  else names(theta)<- c(colnames(x),"sigma2","phiCS",paste0("Dsqrt",1:length(dd)),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                                           phi=phiCS,dsqrt=dd),
                  uhat=unlist(res_emj$uj))
  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi
  

  if (informa) {
    desvios<-try(InfmatrixCS(y,x,z,ind,beta1,sigmae,phiCS,D1,distr = distr,nu = nu),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
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
#
predictf.CS<- function(formFixed,formRandom,dataFit,dataPred,groupVar,distr,theta){
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
  if (distr=="sn") {nu<- NULL}
  if (distr=="st") {nu<- theta[p+q2+3]}
  if (distr=="ss") {nu<- theta[p+q2+3]}
  if (distr=="scn") {nu<- theta[(p+q2+3):(p+q2+4)]}
  if ((p+2+q2+length(nu))!=length(theta)) stop("theta misspecified")
  D1sqrt <- Dmatrix(dd)
  D1 <- D1sqrt%*%D1sqrt
  #
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
    medFit <- x1%*%beta1
    medPred <- x1Pred%*%beta1
    #
    y1=dataFitj[,all.vars(formFixed)[1]]
    SigmaPlus = sigmae*CovCS(phiCS,njFit+njPred)
    PsiPlus<-(zPlus1)%*%(D1)%*%t(zPlus1)+SigmaPlus
    dj<-as.numeric(t(y1-medFit)%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit))
    LambdaPlus <- solve(solve(D1)+ t(zPlus1)%*%solve(SigmaPlus)%*%zPlus1)
    sPsiPlus <- solve(PsiPlus)
    Psi22.1 <- PsiPlus[seqPred,seqPred]- PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]
    ypredj <- medPred + PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit)
    ypred[dataPred$ind==indj] <- ypredj
    xpred[dataPred$ind==indj,] <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
  }
  colnames(xpred) = colnames(xPlus1)
  if (all(xpred[,1]==1)) xpred=xpred[,-1]
  data.frame(groupVar=dataPred$ind,xpred,ypred)
}
#
# ################################################################
# #Log-likelihood - DEC
# ################################################################
ljnormalDECs <-function(j,y,x,z,time,beta1,D1,sigmae,phiDEC,thetaDEC){
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,thetaDEC,t1)#CovARp(phiAR,t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log(dmvnorm(y1,med,Psi))
}
#
ljtDECs <-function(j,nu,y,x,z,time,beta1,D1,sigmae,phiDEC,thetaDEC){
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,thetaDEC,t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  log(dtj)
}
# #
ljsDECs <-function(j,nu,y,x,z,time,beta1,D1,sigmae,phiDEC,thetaDEC){
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,thetaDEC,t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))
  resp <- integrate(Vectorize(f2),0,1)$value
  log(nu*resp)
}
# #
ljcnDECs <-function(j,nu,y,x,z,time,beta1,D1,sigmae,phiDEC,thetaDEC){
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,thetaDEC,t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
        (1-nu[1])*dmvnorm(y1,med,Psi))
}
#
logveroDECs = function(y,x,z,time,ind,beta1,sigmae,phiDEC,thetaDEC,D1,distr,nu){ #ind = indicadora de individuo
  N <-length(ind)

  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalDECs,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtDECs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsDECs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnDECs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  lv
}
# ##############################################################################
# # EM - DEC
# ##############################################################################
calcbi_emjDECs = function(jseq, y, x, z,time, beta1, D1, sigmae,phiDEC,thetaDEC,distr,nu) {
  y1=y[jseq]
  t1=time[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1
  nj = length(y1)
  Sigma = sigmae*CovDEC(phiDEC,thetaDEC,t1)#CovARp(phi = estphit(piAR),t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  #
  bi<-D1%*%t(z1)%*%solve(Psi)%*%(y1-med)
  bi
}
emjDECs = function(jseq, y, x, z,time, beta1, D1, sigmae,phiDEC,thetaDEC,distr,nu) {
  y1=y[jseq]
  t1=time[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1
  nj = length(y1)
  Sigma = sigmae*CovDEC(phiDEC,thetaDEC,t1)#CovARp(phi = estphit(piAR),t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  #
  if  (distr=="sn"){
    uj<-1
  }
  if (distr=="st"){
    uj<-(nj+nu)/(dj+nu)
  }

  if (distr=="ss"){
    uj<-pgamma(1,nj/2+nu+1,dj/2)/pgamma(1,nj/2+nu,dj/2)*(nj+2*nu)/dj
  }

  if (distr=="scn"){
    fy<-as.numeric((nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))
    uj<-as.numeric((nu[1]*nu[2]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))/fy
  }

  sSigma = solve(Sigma)
  sRi = sSigma*sigmae
  Tbj<-solve(solve(D1)+t(z1)%*%sSigma%*%z1)
  r<-Tbj%*%t(z1)%*%sSigma%*%(y1-x1%*%beta1)
  ub<-uj*r
  ub2j<-Tbj+uj*r%*%t(r)
  #
  sum1<-uj*t(x1)%*%sRi%*%x1 #denom beta
  sum2<-(t(x1)%*%sRi%*%(uj*y1-z1%*%ub)) #num beta
  sum3<-uj*t(y1-x1%*%beta1)%*%sRi%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%sRi%*%z1%*%ub-
    t(ub)%*%t(z1)%*%sRi%*%(y1-x1%*%beta1)+traceM(sRi%*%z1%*%ub2j%*%t(z1)) #soma do sig2
  sum4<-ub2j #soma do D1
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,uj=uj,ubj=ub,ub2j=ub2j)
  #if (calcbi) obj.out$bi=bi
  return(obj.out)
}

#
EM.DEC<- function(formFixed,formRandom,data,groupVar,timeVar,
                  beta1,sigmae,D1,distr,nu,parDEC,lb,lu,luDEC,
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
  #
  if (is.null(parDEC)) {
    #cat("calculating initial values for DEC... \n")
    thetat<- seq(0.1,2,by=.1)
    phit <- seq(0.1,.9,by=.05)
    vect <-merge(phit,thetat,all=T)
    logveroDECvs<-function(phitheta){logveroDECs(y, x, z,time,ind, beta1=beta1, sigmae=sigmae,
                                                 phiDEC=phitheta[1],thetaDEC=phitheta[2],
                                                 D1=D1, distr=distr, nu=nu)}
    logverovec <- apply(vect,1,logveroDECvs)
    parDEC <- as.numeric(vect[which.max(logverovec),])
  }
  phiDEC=parDEC[1]
  thetaDEC=parDEC[2]

  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phiDEC,thetaDEC,nu)

  criterio<-10
  count<-0
  llji = logveroDECs(y, x, z, time,ind, beta1, sigmae,phiDEC,thetaDEC, D1, distr, nu)
  if (is.nan(llji)|is.infinite(abs(llji))) stop("NaN/infinity initial likelihood")

  while((criterio > precisao)&(count<max.iter)){
    #print(nu)

    count <- count + 1
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
    #
    parDEC<- optim(c(phiDEC,thetaDEC),lcDEC,gr = NULL,method = "L-BFGS-B", lower =rep(0.0001,2),
                   upper = c(.9999,luDEC),control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
    phiDEC<-parDEC[1]; thetaDEC<-parDEC[2]
    #
    logvero1<-function(nu){logveroDECs(y, x, z,time, ind, beta1, sigmae,phiDEC,thetaDEC, D1, distr, nu)}

    if (distr=="sn"){ nu<-NULL} else
    {nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par}
    #
    param <- teta
    teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phiDEC,thetaDEC,nu)
    criterio2 <- as.numeric(sqrt((teta-param)%*%(teta-param)))
    llj1<-llji
    llji <- logveroDECs(y, x, z, time,ind, beta1, sigmae,phiDEC,thetaDEC, D1, distr, nu)
    criterio <- abs((llji-llj1)/llj1)
    if (showiter&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (count==max.iter) message("\n maximum number of iterations reachead")
  }

  cat("\n")
  bi <- t(bind_cols(tapply(1:N,ind,calcbi_emjDECs,y=y, x=x, z=z,time=time, beta1=beta1, D1=D1,
                           sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC, distr=distr,nu=nu,simplify = FALSE)))
  dd<-matrix.sqrt(D1)[upper.tri(D1, diag = T)]
  theta = c(beta1,sigmae,phiDEC,thetaDEC,dd,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phiDEC","thetaDEC",paste0("Dsqrt",1:length(dd)))
  else names(theta)<- c(colnames(x),"sigma2","phiDEC","thetaDEC",paste0("Dsqrt",1:length(dd)),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                            phi=c(phiDEC,thetaDEC),dsqrt=dd),
                  uhat=unlist(res_emj$uj))
  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixDECs(y,x,z,time,ind,beta1,sigmae,phiDEC,thetaDEC,D1,distr = distr,nu = nu),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
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

predictf.DEC<- function(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr,theta){
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
  if (distr=="sn") {nu<- NULL}
  if (distr=="st") {nu<- theta[p+q2+4]}
  if (distr=="ss") {nu<- theta[p+q2+4]}
  if (distr=="scn") {nu<- theta[(p+q2+4):(p+q2+5)]}
  if ((p+3+q2+length(nu))!=length(theta)) stop("theta misspecified")
  D1sqrt <- Dmatrix(dd)
  D1 <- D1sqrt%*%D1sqrt
  #
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
    medFit <- x1%*%beta1
    medPred <- x1Pred%*%beta1
    #
    y1=dataFitj[,all.vars(formFixed)[1]]
    SigmaPlus = sigmae*CovDEC(phiDEC,thetaDEC,c(dataFitj$time,dataPredj$time))
    PsiPlus<-(zPlus1)%*%(D1)%*%t(zPlus1)+SigmaPlus
    dj<-as.numeric(t(y1-medFit)%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit))
    LambdaPlus <- solve(solve(D1)+ t(zPlus1)%*%solve(SigmaPlus)%*%zPlus1)
    sPsiPlus <- solve(PsiPlus)
    Psi22.1 <- PsiPlus[seqPred,seqPred]- PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]
    ypredj <- medPred + PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit)
    ypred[dataPred$ind==indj] <- ypredj
    xpred[dataPred$ind==indj,] <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
    timepred[dataPred$ind==indj] <- dataPredj$time
  }
  colnames(xpred) = colnames(xPlus1)
  if (all(xpred[,1]==1)) xpred=xpred[,-1]
  data.frame(groupVar=dataPred$ind,time=timepred,xpred,ypred)
}
#
# ################################################################
# #Log-likelihood - CAR(1)
# ################################################################
ljnormalCAR1s <-function(j,y,x,z,time,beta1,D1,sigmae,phiDEC){
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,1,t1)#CovARp(phiAR,t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log(dmvnorm(y1,med,Psi))
}
# #
ljtCAR1s <-function(j,nu,y,x,z,time,beta1,D1,sigmae,phiDEC){
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,1,t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  #log(2*dmvt(y1,delta = med, sigma = Psi, df = nu,log=F)*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))#veroST1(Psi,Ajj,dj,nu,pp=njj))
  log(dtj)
}
#
ljsCAR1s <-function(j,nu,y,x,z,time,beta1,D1,sigmae,phiDEC){
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,1,t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))
  resp <- integrate(Vectorize(f2),0,1)$value
  log(nu*resp)
}
#
ljcnCAR1s <-function(j,nu,y,x,z,time,beta1,D1,sigmae,phiDEC){
  y1=y[j]
  t1= time[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Sigma=sigmae*CovDEC(phiDEC,1,t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
        (1-nu[1])*dmvnorm(y1,med,Psi))
}

logveroCAR1s = function(y,x,z,time,ind,beta1,sigmae,phiDEC,D1,distr,nu){ #ind = indicadora de individuo
  N <-length(ind)

  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalCAR1s,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtCAR1s,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsCAR1s,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnCAR1s,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC))
  lv
}
##############################################################################
# EM - CAR(1)
##############################################################################
emjCAR1s = function(jseq, y, x, z,time, beta1, D1, sigmae,phiDEC,distr,nu) {
  y1=y[jseq]
  t1=time[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1
  nj = length(y1)
  Sigma = sigmae*CovDEC(phiDEC,1,t1)#CovARp(phi = estphit(piAR),t1)
  Psi<-(z1)%*%(D1)%*%t(z1)+Sigma
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  #
  if  (distr=="sn"){
    uj<-1
  }
  if (distr=="st"){
    uj<-(nj+nu)/(dj+nu)
  }

  if (distr=="ss"){
    uj<-pgamma(1,nj/2+nu+1,dj/2)/pgamma(1,nj/2+nu,dj/2)*(nj+2*nu)/dj
  }

  if (distr=="scn"){
    fy<-as.numeric((nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))
    uj<-as.numeric((nu[1]*nu[2]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))/fy
  }

  sSigma = solve(Sigma)
  sRi = sSigma*sigmae
  Tbj<-solve(solve(D1)+t(z1)%*%sSigma%*%z1)
  r<-Tbj%*%t(z1)%*%sSigma%*%(y1-x1%*%beta1)
  ub<-uj*r
  ub2j<-Tbj+uj*r%*%t(r)
  #
  sum1<-uj*t(x1)%*%sRi%*%x1 #denom beta
  sum2<-(t(x1)%*%sRi%*%(uj*y1-z1%*%ub)) #num beta
  sum3<-uj*t(y1-x1%*%beta1)%*%sRi%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%sRi%*%z1%*%ub-
    t(ub)%*%t(z1)%*%sRi%*%(y1-x1%*%beta1)+traceM(sRi%*%z1%*%ub2j%*%t(z1)) #soma do sig2
  sum4<-ub2j #soma do D1
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,uj=uj,ubj=ub,ub2j=ub2j)
  #if (calcbi) obj.out$bi=bi
  return(obj.out)
}
#
EM.CAR1<- function(formFixed,formRandom,data,groupVar,timeVar,
                   distr,beta1,sigmae,phiCAR1,D1,nu,lb,lu,
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
  #
  if (is.null(phiCAR1)) {
    lmeCAR = try(lme(formFixed,random=~1|ind,data=data,correlation=corCAR1(form = ~time)),silent=T)
    if (class(lmeCAR)=="try-error") phiDEC =abs(as.numeric(pacf(y-x%*%beta1,lag.max=1,plot=F)$acf))
    else {
      phiDEC = capture.output(lmeCAR$modelStruct$corStruct)[3]
      phiDEC = as.numeric(strsplit(phiDEC, " ")[[1]])
    }
  } else phiDEC <- phiCAR1

  teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phiDEC,nu)

  criterio<-10
  count<-0
  llji = logveroCAR1s(y, x, z, time,ind, beta1, sigmae,phiDEC, D1, distr, nu)
  if (is.nan(llji)|is.infinite(abs(llji))) stop("NaN/infinity initial likelihood")

  while((criterio > precisao)&(count<max.iter)){
    #print(nu)

    count <- count + 1
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
    #
    phiDEC<- optim(phiDEC,lcCAR1,gr = NULL,method = "L-BFGS-B", lower =0.0001,
                   upper = .9999,control = list(fnscale=-1),beta1=beta1,sigmae=sigmae,
                   y=y,x=x,z=z,time=time,ind=ind,u=res_emj$uj,ub=res_emj$ubj,ub2=res_emj$ub2j)$par
    #
    logvero1<-function(nu){logveroCAR1s(y, x, z,time, ind, beta1, sigmae,phiDEC, D1, distr, nu)}

    if (distr=="sn"){ nu<-NULL} else
    {nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par}
    #
    param <- teta
    teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],phiDEC,nu)
    criterio2 <- as.numeric(sqrt((teta-param)%*%(teta-param)))
    llj1<-llji
    llji <- logveroCAR1s(y, x, z, time,ind, beta1, sigmae,phiDEC, D1, distr, nu)
    criterio <- abs((llji-llj1)/llj1)
    if (showiter&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio,"\r") #  criterium ",criterio," or ",criterio2,"\r")
    if (count==max.iter) message("\n maximum number of iterations reachead")
  }

  cat("\n")
  bi <- t(bind_cols(tapply(1:N,ind,calcbi_emjDECs,y=y, x=x, z=z,time=time, beta1=beta1, D1=D1,
                           sigmae=sigmae,phiDEC=phiDEC,thetaDEC=1, distr=distr,nu=nu,simplify = FALSE)))
  dd<-matrix.sqrt(D1)[upper.tri(D1, diag = T)]
  theta = c(beta1,sigmae,phiDEC,dd,nu)
  if (is.null(colnames(x))) colnames(x) <- paste0("beta",1:p-1)
  if (distr=="sn") names(theta)<-c(colnames(x),"sigma2","phiCAR1",paste0("Dsqrt",1:length(dd)))
  else names(theta)<- c(colnames(x),"sigma2","phiCAR1",paste0("Dsqrt",1:length(dd)),paste0("nu",1:length(nu)))

  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                                           phi=phiDEC,dsqrt=dd),
                  uhat=unlist(res_emj$uj))
  if (distr != "sn") obj.out$estimates$nu = nu
  colnames(bi) <- colnames(z)
  obj.out$random.effects<- bi

  if (informa) {
    desvios<-try(InfmatrixCAR1s(y,x,z,time,ind,beta1,sigmae,phiDEC,D1,distr = distr,nu = nu),silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
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
#
predictf.CAR1<- function(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr,theta){
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
  if (distr=="sn") {nu<- NULL}
  if (distr=="st") {nu<- theta[p+q2+3]}
  if (distr=="ss") {nu<- theta[p+q2+3]}
  if (distr=="scn") {nu<- theta[(p+q2+3):(p+q2+4)]}
  if ((p+2+q2+length(nu))!=length(theta)) stop("theta misspecified")
  D1sqrt <- Dmatrix(dd)
  D1 <- D1sqrt%*%D1sqrt
  #
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
    medFit <- x1%*%beta1
    medPred <- x1Pred%*%beta1
    #
    y1=dataFitj[,all.vars(formFixed)[1]]
    SigmaPlus = sigmae*CovDEC(phiDEC,1,c(dataFitj$time,dataPredj$time))
    PsiPlus<-(zPlus1)%*%(D1)%*%t(zPlus1)+SigmaPlus
    dj<-as.numeric(t(y1-medFit)%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit))
    LambdaPlus <- solve(solve(D1)+ t(zPlus1)%*%solve(SigmaPlus)%*%zPlus1)
    sPsiPlus <- solve(PsiPlus)
    Psi22.1 <- PsiPlus[seqPred,seqPred]- PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]
    ypredj <- medPred + PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit)
    ypred[dataPred$ind==indj] <- ypredj
    xpred[dataPred$ind==indj,] <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
    timepred[dataPred$ind==indj] <- dataPredj$time
  }
  colnames(xpred) = colnames(xPlus1)
  if (all(xpred[,1]==1)) xpred=xpred[,-1]
  data.frame(groupVar=dataPred$ind,time=timepred,xpred,ypred)
}
####
#inf mat (using first derivative and codes from the asymmetric case)
#Information matrix for SMN-LMM and SMN-LMM-AR(p) with E(bi)=0
scoreis <- function(jseq,y,x,z,beta1,sigmae,D1,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  lambda=rep(0,nrow(D1))
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


Infmatrixs <- function(y,x,z,ind,beta1,sigmae,D1,distr,nu){
  N <-length(y)
  score_list=tapply(1:N,ind,scoreis,y=y, x=x, z=z, beta1=beta1, sigmae=sigmae,D1=D1,distr=distr,nu=nu)
  mi_list = lapply(score_list,function(tt) {xm = matrix(tt,ncol=1);xm%*%t(xm)})
  infmat <- Reduce("+",mi_list)
  npar <- nrow(infmat);nd<-nrow(D1)
  infmat<- infmat[1:(npar-nd),1:(npar-nd)]
  if (abs(det(infmat))<1e-5) infmat= infmat+1e-20*diag(nrow(infmat))
  sqrt(diag(solve(infmat)))
}

#########################################################################
#ar(p)
#########################################################################
scoreARis <- function(jseq,y,x,z,time,beta1,sigmae,phiAR,D1,distr,nu) {
  lambda=rep(0,nrow(D1))
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


InfmatrixARs <- function(y,x,z,time,ind,beta1,sigmae,phiAR,D1,distr,nu){
  N <-length(y)
  score_list=tapply(1:N,ind,scoreARis,y=y, x=x, z=z,time=time, beta1=beta1, sigmae=sigmae,phiAR=phiAR,D1=D1,distr=distr,nu=nu)
  mi_list = lapply(score_list,function(tt) {xm = matrix(tt,ncol=1);xm%*%t(xm)})
  infmat <- Reduce("+",mi_list)
  npar <- nrow(infmat);nd<-nrow(D1)
  infmat<- infmat[1:(npar-nd),1:(npar-nd)]
  if (abs(det(infmat))<1e-5) infmat= infmat+1e-10*diag(nrow(infmat))
  sqrt(diag(solve(infmat)))
}

#########################################################################
#CS
#########################################################################
scoreCSis <- function(jseq,y,x,z,beta1,sigmae,phiCS,D1,distr,nu) {
  lambda=rep(0,nrow(D1))
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


InfmatrixCSs <- function(y,x,z,ind,beta1,sigmae,phiCS,D1,distr,nu){
  N <-length(y)
  score_list=tapply(1:N,ind,scoreCSis,y=y, x=x, z=z, beta1=beta1, sigmae=sigmae,phiCS=phiCS,D1=D1,distr=distr,nu=nu)
  mi_list = lapply(score_list,function(tt) {xm = matrix(tt,ncol=1);xm%*%t(xm)})
  infmat <- Reduce("+",mi_list)
  npar <- nrow(infmat);nd<-nrow(D1)
  infmat<- infmat[1:(npar-nd),1:(npar-nd)]
  if (abs(det(infmat))<2e-5) infmat= infmat+1e-10*diag(nrow(infmat))
  sqrt(diag(solve(infmat)))
}

#########################################################################
#DEC
#########################################################################
scoreDECis <- function(jseq,y,x,z,time,beta1,sigmae,phiDEC,thetaDEC,D1,distr,nu) {
  lambda = rep(0,nrow(D1))
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


InfmatrixDECs <- function(y,x,z,time,ind,beta1,sigmae,phiDEC,thetaDEC,D1,distr,nu){
  N <-length(y)
  score_list=tapply(1:N,ind,scoreDECis,y=y, x=x, z=z,time=time, beta1=beta1, sigmae=sigmae,
                    phiDEC=phiDEC,thetaDEC=thetaDEC,D1=D1,distr=distr,nu=nu)
  mi_list = lapply(score_list,function(tt) {xm = matrix(tt,ncol=1);xm%*%t(xm)})
  infmat <- Reduce("+",mi_list)
  npar <- nrow(infmat);nd<-nrow(D1)
  infmat<- infmat[1:(npar-nd),1:(npar-nd)]
  if (abs(det(infmat))<1e-5) infmat= infmat+1e-10*diag(nrow(infmat))
  sqrt(diag(solve(infmat)))
}

#########################################################################
#CAR(1)
#########################################################################
scoreCAR1is <- function(jseq,y,x,z,time,beta1,sigmae,phiDEC,D1,distr,nu) {
  lambda=rep(0,nrow(D1))
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

InfmatrixCAR1s <- function(y,x,z,time,ind,beta1,sigmae,phiDEC,D1,distr,nu){
  N <-length(y)
  score_list=tapply(1:N,ind,scoreCAR1is,y=y, x=x, z=z,time=time, beta1=beta1, sigmae=sigmae,
                    phiDEC=phiDEC,D1=D1,distr=distr,nu=nu)
  mi_list = lapply(score_list,function(tt) {xm = matrix(tt,ncol=1);xm%*%t(xm)})
  infmat <- Reduce("+",mi_list)
  npar <- nrow(infmat);nd<-nrow(D1)
  infmat<- infmat[1:(npar-nd),1:(npar-nd)]
  if (abs(det(infmat))<1e-5) infmat= infmat+1e-10*diag(nrow(infmat))
  sqrt(diag(solve(infmat)))
}


