residuals.SMSN<- function(object,level="conditional",type="response",...){
  if (!(level %in% c("marginal","conditional"))) stop("Accepted levels: marginal, conditional")
  if (!(type %in% c("response","modified","normalized"))) stop("Accepted types: response, normalized or modified")
  data <- object$data
  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<-object$groupVar
  timeVar <- object$timeVar
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  if (!is.null(timeVar)) {
    time<-data[,timeVar]
  } else{
    time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  }
  p<-ncol(x)
  q1<-ncol(z)
  N<- nrow(data)
  ind_levels <- levels(ind)
  distr <- object$distr
  #depStruct <- object$depStruct
  #
  res <- numeric(N)
  if (type=="response") {
    if (level=="marginal") {
      lab <- "marginal raw residuals"
      for (i in seq_along(ind_levels)) {
        seqi <- ind==ind_levels[i]
        xfiti <- matrix(x[seqi,],ncol=p)
        zfiti <- matrix(z[seqi,],ncol=q1)
        res[seqi]<- y[seqi]- xfiti%*%object$estimates$beta
      }
    } else{
      lab <- "conditional raw residuals"
      for (i in seq_along(ind_levels)) {
        seqi <- ind==ind_levels[i]
        xfiti <- matrix(x[seqi,],ncol=p)
        zfiti <- matrix(z[seqi,],ncol=q1)
        res[seqi]<- y[seqi]- (xfiti%*%object$estimates$beta + zfiti%*%object$random.effects[i,])
      }
    }
  }
  if (type=="modified"){
    if (level=="marginal") {
      lab <- "marginal modified residuals"
      Dest <- Dmatrix(object$estimates$dsqrt)%*%Dmatrix(object$estimates$dsqrt)
      for (i in seq_along(ind_levels)) {
        seqi <- ind==ind_levels[i]
        xfiti <- matrix(x[seqi,],ncol=p)
        zfiti <- matrix(z[seqi,],ncol=q1)
        timei <- time[seqi]
        Sigmaest <- errorVar(timei,object)
        vary <- Sigmaest+ (zfiti)%*%Dest%*%t(zfiti)
        sigFitinv <- matrix.sqrt(solve(vary))
        res[seqi]<- sigFitinv%*%(y[seqi]- xfiti%*%object$estimates$beta)
      }
    } else{
      lab <- "conditional modified residuals"
      for (i in seq_along(ind_levels)) {
        seqi <- ind==ind_levels[i]
        xfiti <- matrix(x[seqi,],ncol=p)
        zfiti <- matrix(z[seqi,],ncol=q1)
        timei <- time[seqi]
        Sigmaest <- errorVar(timei,object)
        sigeFitinv <- matrix.sqrt(solve(Sigmaest))
        res[seqi]<- sigeFitinv%*%(y[seqi]- (xfiti%*%object$estimates$beta + zfiti%*%object$random.effects[i,]))
      }
    }
  }
  if (type=="normalized"){
    if (distr=="st"|distr=="t") if (object$estimates$nu<=2) stop("normalized residual not defined for nu<=2")
    if (distr=="ssl"|distr=="sl") if (object$estimates$nu<=1) stop("normalized residual not defined for nu<=1")
    if (distr=="sn"|distr=="norm") {c.=-sqrt(2/pi);k2=1}
    if (distr=="st"|distr=="t") {c.=-sqrt(object$estimates$nu/pi)*
      gamma((object$estimates$nu-1)/2)/gamma(object$estimates$nu/2);
          k2=object$estimates$nu/(object$estimates$nu-2)}
    if (distr=="ssl"|distr=="sl") {c.=-sqrt(2/pi)*object$estimates$nu/(object$estimates$nu-.5);
          k2=object$estimates$nu/(object$estimates$nu-1)}
    if (distr=="scn"|distr=="cn") {c.=-sqrt(2/pi)*(1+object$estimates$nu[1]*
                                         (object$estimates$nu[2]^(-.5)-1));
          k2=1+object$estimates$nu[1]*(object$estimates$nu[2]^(-1)-1)}
    if (level=="marginal") {
      if (inherits(object,"SMN")) {
        object$estimates$lambda <- rep(0,q1);c.=0
      }
      Dest <- Dmatrix(object$estimates$dsqrt)%*%Dmatrix(object$estimates$dsqrt)
      delta = object$estimates$lambda/as.numeric(
        sqrt(1+t(object$estimates$lambda)%*%(object$estimates$lambda)))
      Delta = Dmatrix(object$estimates$dsqrt)%*%delta
      lab <- "marginal standardized residuals"
      for (i in seq_along(ind_levels)) {
        seqi <- ind==ind_levels[i]
        xfiti <- matrix(x[seqi,],ncol=p)
        zfiti <- matrix(z[seqi,],ncol=q1)
        timei <- time[seqi]
        Sigmaest <- errorVar(timei,object)
        vary <- Sigmaest+ (zfiti)%*%Dest%*%t(zfiti)
        sigFitinv <- matrix.sqrt(solve(k2*vary - c.^2*zfiti%*%Delta%*%t(zfiti%*%Delta)))
        res[seqi]<- sigFitinv%*%(y[seqi]- xfiti%*%object$estimates$beta)
      }
    } else{
      lab <- "conditional standardized residuals"
      for (i in seq_along(ind_levels)) {
        seqi <- ind==ind_levels[i]
        xfiti <- matrix(x[seqi,],ncol=p)
        zfiti <- matrix(z[seqi,],ncol=q1)
        timei <- time[seqi]
        Sigmaest <- errorVar(timei,object)
        sigeFitinv <- matrix.sqrt(solve(k2*Sigmaest))
        res[seqi]<- sigeFitinv%*%(y[seqi]- (xfiti%*%object$estimates$beta + zfiti%*%object$random.effects[i,]))
      }
    }
  }
  attr(res, "label") <- lab
  res
}

residuals.SMN<- residuals.SMSN

acfresid <- function(object,maxLag,resLevel="marginal",resType="normalized",
                     calcCI=FALSE,IClevel,MCiter,seed){
  if (!missing(seed)) set.seed(seed)
  if(!inherits(object,c("SMSN","SMN"))) stop("object must inherit from class SMSN or SMN")
  #
  #object$data <- object$data[order(object$data[,object$groupVar]),]
  res <- residuals(object,level=resLevel,type = resType)
  data <- object$data
  data$res <- res
  data <- data[order(data[,object$groupVar]),]
  res<-data$res
  groupVar<-object$groupVar
  timeVar <- object$timeVar
  ind <-data[,groupVar]
  if (!is.null(timeVar)) {
    time<-data[,timeVar]
  } else{
    time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  }
  if (any(!is.wholenumber(time))) {
    time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
    warning("time must contain positive integer numbers, timeVar will be ignored and the order will be used instead")
  }
  if(missing(maxLag)) {
    maxLag <- min(c(maxL <- max(tapply(time,ind, function(x) diff(range(x)))),
                    as.integer(10 * log10(maxL + 1))))
  } else{
    if (maxLag> (maxL <- max(tapply(time,ind, function(x) diff(range(x)))))) maxLag <-
        min(c(maxL,as.integer(10 * log10(maxL + 1))))
  }
  if (calcCI) {
    if ((resType =="response")|(resLevel!="marginal"))
      stop("ICs are obtained through an empirical method that is only
      appropriate for marginal modified/normalized residuals, please use
      (resType='normalized' or resType='modified') and resLevel='marginal',
      or calcCI=FALSE")
    if (!missing(IClevel)) {
      if (IClevel>=1|IClevel<=0) stop("0<IClevel<1 needed")
    } else IClevel <- .95
    if (!missing(MCiter)) {
      if (!is.wholenumber(MCiter)) stop("MCiter must be an integer positive number")
      if (MCiter<=0) stop("MCiter must be an integer positive number")
    } else MCiter<- 300
  }
  N<- nrow(data)
  ind_levels <- levels(ind)
  #
  corind <- nind<-matrix(ncol=maxLag+1,nrow=length(ind_levels))
  for (i in seq_along(ind_levels)) {
    seqi <- ind==ind_levels[i]
    resim <-res[seqi]
    timei <- time[seqi]
    resi <- rep(NA,diff(range(timei))+1)
    resi[timei] <- resim
    corind[i,1] <- sum(resi^2,na.rm=T); nind[i,1]<- sum(!is.na(resi))
    kmaxi <- min(maxLag,diff(range(timei)))#nind[i,1]-1)
    for (k in seq_len(kmaxi)) {
      nind[i,k+1]<- sum(!is.na(diff(resi,k)))
      corind[i,k+1] <- sum(resi[-seq_len(k)]*
                             resi[-seq.int(length(resi),by=-1,length.out = k)],na.rm=T)
    }
  } #
  nvec <- apply(nind, 2, sum,na.rm=T)
  cormean <- apply(corind, 2, sum,na.rm=T)/nvec
  corvec<-cormean/cormean[1]
  out <- data.frame(lag=0:maxLag,ACF=corvec,n.used=nvec)
  ##############################
  #MC IC
  if (calcCI) {
    x <- model.matrix(object$formula$formFixed,data=data)
    z<-model.matrix(object$formula$formRandom,data=data)
    p<-ncol(x)
    q1<-ncol(z)
    distr <- object$distr
    #
    acfMC <- matrix(nrow=MCiter,ncol=maxLag)
    #Dfit <- Dmatrix(object$estimates$dsqrt)
    Dfit2 <- object$estimates$D #Dfit%*%Dfit
    Dfit <- matrix.sqrt(Dfit2)
    #sigma2<-object$estimates$sigma2
    sigma2<-as.numeric(errorVar(1,object))
    lambda <- object$estimates$lambda
    if (is.null(lambda)) lambda<- rep(0,nrow(Dfit))
    beta<- object$estimates$beta
    nu <- object$estimates$nu
    #
    if (resType=="modified") {#modified residual
      for (sample in 1:MCiter) {
        dadosi = tapply(1:N,ind,gerar_ind_smsnACF,x=x,z=z,sigma2=sigma2,Dsqrti=Dfit,
                        beta1=beta,lambda=lambda,distr=distr,nu=nu,ind=ind,time=time) %>% bind_rows()
        #dadosi$ind <- ind
        #calculating residuals
        resMC<- numeric(N)
        for (i in seq_along(ind_levels)) {
          seqi <- dadosi$ind==ind_levels[i]#;print(sum(seqi))
          xfiti <- matrix(x[seqi,],ncol=p)
          zfiti <- matrix(z[seqi,],ncol=q1)
          Sigmaest <- sigma2*diag(sum(seqi))#errorVar(timei,object)
          vary <- Sigmaest+ (zfiti)%*%Dfit2%*%t(zfiti)
          sigFitinv <- matrix.sqrt(solve(vary))
          resMC[seqi]<- sigFitinv%*%(dadosi$y[seqi]- xfiti%*%beta)
        }
        corindi <- nindi<-matrix(ncol=maxLag+1,nrow=length(ind_levels))
        for (i in seq_along(ind_levels)) {
          seqi <- dadosi$ind==ind_levels[i]
          resim <-resMC[seqi]
          timei <- dadosi$time[seqi]
          resi <- rep(NA,diff(range(timei))+1)
          resi[timei] <- resim
          corindi[i,1] <- sum(resi^2,na.rm=T); nindi[i,1]<- sum(!is.na(resi))
          kmaxi <- min(maxLag,diff(range(timei)))#nind[i,1]-1)
          for (k in seq_len(kmaxi)) {
            nindi[i,k+1]<- sum(!is.na(diff(resi,k)))
            corindi[i,k+1] <- sum(resi[-seq_len(k)]*
                                    resi[-seq.int(length(resi),by=-1,length.out = k)],na.rm=T)
          }
        } #
        #ormeani <-apply(corindi/nindi,2,sum,na.rm=T)
        cormeani <- apply(corindi, 2, sum,na.rm=T)/apply(nindi, 2, sum,na.rm=T)
        acfMC[sample,] <-cormeani[-1]/cormeani[1]
      }
      MCacfIC<-apply(acfMC,2,function(x) quantile(x,probs = c((1-IClevel)/2,(1+IClevel)/2)))
      out <- cbind(out,IC = rbind(rep(NA,2),t(MCacfIC)))
    } else { #normalized residual
      if (distr=="st"|distr=="t") if (nu<=2) stop("normalized residual not defined for nu<=2")
      if (distr=="ssl"|distr=="sl") if (nu<=1) stop("normalized residual not defined for nu<=1")
      if (distr=="sn"|distr=="norm") {c.=-sqrt(2/pi);k2=1}
      if (distr=="st"|distr=="t") {c.=-sqrt(nu/pi)*
        gamma((nu-1)/2)/gamma(nu/2);
      k2=nu/(nu-2)}
      if (distr=="ssl"|distr=="sl") {c.=-sqrt(2/pi)*nu/(nu-.5);
      k2=nu/(nu-1)}
      if (distr=="scn"|distr=="cn") {c.=-sqrt(2/pi)*(1+nu[1]*
                                                       (nu[2]^(-.5)-1));
      k2=1+nu[1]*(nu[2]^(-1)-1)}
      if (inherits(object,"SMN")) {
        c.=0
      }
      delta = lambda/as.numeric(sqrt(1+t(lambda)%*%(lambda)))
      Delta = Dfit%*%delta
      for (sample in 1:MCiter) {
        dadosi = tapply(1:N,ind,gerar_ind_smsnACF,x=x,z=z,sigma2=sigma2,Dsqrti=Dfit,
                        beta1=beta,lambda=lambda,distr=distr,nu=nu,ind=ind,time=time) %>% bind_rows()
        #calculating residuals
        resMC<- numeric(N)
        for (i in seq_along(ind_levels)) {
          seqi <- dadosi$ind==ind_levels[i]#;print(sum(seqi))
          xfiti <- matrix(x[seqi,],ncol=p)
          zfiti <- matrix(z[seqi,],ncol=q1)
          Sigmaest <- sigma2*diag(sum(seqi))#errorVar(timei,object)
          vary <- Sigmaest+ (zfiti)%*%Dfit2%*%t(zfiti)
          sigFitinv <- matrix.sqrt(solve(k2*vary - c.^2*zfiti%*%Delta%*%t(zfiti%*%Delta)))
          #sigFitinv <- matrix.sqrt(solve(vary))
          resMC[seqi]<- sigFitinv%*%(dadosi$y[seqi]- xfiti%*%beta)
        }
        corindi <- nindi<-matrix(ncol=maxLag+1,nrow=length(ind_levels))
        for (i in seq_along(ind_levels)) {
          seqi <- dadosi$ind==ind_levels[i]
          resim <-resMC[seqi]
          timei <- dadosi$time[seqi]
          resi <- rep(NA,diff(range(timei))+1)
          resi[timei] <- resim
          corindi[i,1] <- sum(resi^2,na.rm=T); nindi[i,1]<- sum(!is.na(resi))
          kmaxi <- min(maxLag,diff(range(timei)))#nind[i,1]-1)
          for (k in seq_len(kmaxi)) {
            nindi[i,k+1]<- sum(!is.na(diff(resi,k)))
            corindi[i,k+1] <- sum(resi[-seq_len(k)]*
                                    resi[-seq.int(length(resi),by=-1,length.out = k)],na.rm=T)
          }
        } #
        cormeani <- apply(corindi, 2, sum,na.rm=T)/apply(nindi, 2, sum,na.rm=T)
        acfMC[sample,] <-cormeani[-1]/cormeani[1]
      }
      MCacfIC<-apply(acfMC,2,function(x) quantile(x,probs = c((1-IClevel)/2,(1+IClevel)/2)))
      out <- cbind(out,IC = rbind(rep(NA,2),t(MCacfIC)))
    }
  }
  class(out) <- c("acfresid","data.frame")
  out
}

mahalDist<- function(object,decomposed=FALSE,dataPlus=NULL){
  if(!inherits(object,c("SMSN","SMN"))) stop("object must inherit from class SMSN or SMN")
  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<-object$groupVar
  timeVar <- object$timeVar
  if (!is.null(dataPlus)) {
    data <- dataPlus
    if (!is.data.frame(data)) stop("data must be a data.frame")
    vars_used<-unique(c(all.vars(formFixed),all.vars(formRandom),groupVar,timeVar))
    vars_miss <- which(!(vars_used %in% names(data)))
    if (length(vars_miss)>0) stop(paste(vars_used[vars_miss],"not found in dataPlus"))
  } else data <- object$data
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  if (!is.null(timeVar)) {
    time<-data[,timeVar]
  } else{
    time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  }
  p<-ncol(x)
  q1<-ncol(z)
  N<- nrow(data)
  ind_levels <- levels(ind)
  #
  if (inherits(object,"SMN")) {
    object$estimates$lambda <- rep(0,q1); c.<-1
  }
  #
  distr <- object$distr
  if (distr=="sn") {c.=-sqrt(2/pi)}
  if (distr=="st") {c.=-sqrt(object$estimates$nu/pi)*
                  gamma((object$estimates$nu-1)/2)/gamma(object$estimates$nu/2)}
  if (distr=="ssl") {c.=-sqrt(2/pi)*object$estimates$nu/(object$estimates$nu-.5)}
  if (distr=="scn") {c.=-sqrt(2/pi)*(1+object$estimates$nu[1]*
                                       (object$estimates$nu[2]^(-.5)-1))}
  #
  mahaldist <- distbi<- distei<- numeric(length(ind_levels))
  Dest <- object$estimates$D
  delta = object$estimates$lambda/as.numeric(
    sqrt(1+t(object$estimates$lambda)%*%(object$estimates$lambda)))
  Delta = matrix.sqrt(Dest)%*%delta
  for (i in seq_along(ind_levels)) {
    seqi <- ind==ind_levels[i]
    #ibi <- which(row.names(object$random.effects)==ind_levels[i])
    xfiti <- x[seqi,]
    zfiti <- z[seqi,]
    timei <- time[seqi]
    Sigmaest <- errorVar(timei,object)
    Psiy <- Sigmaest+ (zfiti)%*%Dest%*%t(zfiti)
    #
    ytil <- y[seqi]-xfiti%*%object$estimates$beta-c.*zfiti%*%Delta
    if (decomposed) {
      mub <- Dest%*%t(zfiti)%*%solve(Psiy)%*%ytil+c.*Delta
      ei <- y[seqi]-xfiti%*%object$estimates$beta-zfiti%*%mub
      distei[i]<-t(ei)%*%solve(Sigmaest)%*%ei
      distbi[i]<-t(mub-c.*Delta)%*%solve(Dest)%*%(mub-c.*Delta)
    } else{
      mahaldist[i] <- t(ytil)%*%solve(Psiy)%*%ytil
    }
  }
  if (decomposed) {
    out <- data.frame(md.error=distei,md.b=distbi,md=distei+distbi)
    row.names(out) <- row.names(object$random.effects)
    class(out) <- c("mahalDist","data.frame")
  } else {
    out<- mahaldist
    names(out) <- row.names(object$random.effects)
    class(out) <- c("mahalDist","numeric")
  }
  attr(out,'call')<-match.call()
  out
}
###
#IC "monte carlo" para acf
gerar_ind_smsnACF = function(jvec,x,z,sigma2,Dsqrti,beta1,lambda,distr,nu,ind,time) {
  if (distr=="sn"|distr=="norm") {ui=1; c.=-sqrt(2/pi)}
  if (distr=="st"|distr=="t") {ui=rgamma(1,nu/2,nu/2); c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ssl"|distr=="sl") {ui=rbeta(1,nu,1); c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn"|distr=="cn") {ui=ifelse(runif(1)<nu[1],nu[2],1);
                    c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jvec,  ],ncol=p)
  z1=matrix(z[jvec,  ],ncol=q1)
  nj = nrow(x1)
  Sig = sigma2*diag(nj)
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%(lambda)))
  Delta = Dsqrti%*%delta
  Gammab = Dsqrti%*%Dsqrti - Delta%*%t(Delta)
  Beta = matrix(beta1,ncol=1)
  ti = c.+abs(rnorm(1,0,ui^-.5))
  bi = t(rmvnorm(1,Delta*ti,sigma=ui^(-1)*Gammab))
  Yi = t(rmvnorm(1,x1%*%Beta+z1%*%bi,sigma=ui^(-1)*Sig))
  return(data.frame(y=Yi,ind=ind[jvec],time=time[jvec]))
}

plot.acfresid <-  function(x,...) {
  if (ncol(x)<=3) {
    ggplot(data = x,aes_string(x = "lag",y="ACF"))+
      theme_minimal()+geom_hline(aes(yintercept = 0))+
      geom_segment(mapping = aes(xend = lag, yend = 0))+ ylim(c(-1,1))
  } else {
    lagmax <- max(x$lag)
    datIC <- rbind(x[,c(1,4,5)],
                   c(lagmax+1,min(x$`IC.2.5%`,na.rm=T),max(x$`IC.97.5%`,na.rm=T)))
    datIC$lag <- datIC$lag-.5
    names(datIC) <- c("lag","inf","sup")
    ggplot(x,aes_string(x = "lag",y="ACF")) +
    #ggplot(x,aes(x = lag,y=ACF)) +
    geom_hline(aes(yintercept = 0)) +
      coord_cartesian(xlim=c(0,max(x$lag)))+
      geom_segment(mapping = aes_string(xend = "lag", yend = 0)) +
      geom_step(aes_string(x="lag",y="inf"),data=datIC,color=4,linetype="dashed",na.rm=TRUE)+
      geom_step(aes_string(x="lag",y="sup"),data=datIC,color=4,linetype="dashed",na.rm=TRUE)+
      theme_minimal()+ ylim(c(-1,1))
      }
}

#plot obj
plot.SMSN <- function(x,type="response",level="conditional",alpha=.3,...) {
  resid <- residuals(x,type=type,level=level)
  qplot(fitted(x),resid,alpha=I(alpha),...=...) +theme_minimal() +
    geom_hline(yintercept = 0,linetype="dashed") + ylab(attr(resid,"label")) +
    xlab("fitted values")
}
plot.SMN <- plot.SMSN

#plot mahalDist
plot.mahalDist <- function(x,fitobject,type,level=.99,nlabels=3,...){
  if(missing(fitobject)) fitobject<-eval(str2lang(as.character(attr(x,'call')[2])))
  #if(missing(fitobject)) stop("please provide fitobject used to compute the Mahalanobis distance")
  if(missing(type)) {
    type<-"total"
  } else if (!(type %in% c("total","error","b"))) stop("type must be one of the following: total,error,b")
  if ((!is.data.frame(x))&(type!="total")) stop(paste("for type",type,"x must have decomposed=TRUE"))
  if (!is.data.frame(x)) x <- data.frame(md=x)
  if (level>=1|level<=0) stop("0<level<1 needed")
  #
  if(!inherits(fitobject,c("SMSN","SMN"))) stop("fitobject must inherit from class SMSN or SMN")
  data <- fitobject$data
  formFixed <- fitobject$formula$formFixed
  formRandom <- fitobject$formula$formRandom
  groupVar<-fitobject$groupVar
  timeVar <- fitobject$timeVar
  ind <-data[,groupVar]
  if (!is.null(timeVar)) {
    time<-data[,timeVar]
  } else{
    time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  }
  distr <- fitobject$distr
  nu <- fitobject$estimates$nu
  #
  x$nj <- tapply(time,ind,length)
  x$index <- seq_along(x$nj)
  x$ind <- levels(ind)
  #
  if (type!="total") {
    if (n_distinct(x$nj) ==1){
      if (type=="error") {
        plotout<-ggplot(x,aes_string("index","md.error")) +
          geom_point(shape=1) + ylab("error distance")+
          geom_text_repel(aes_string(label="ind"),data=subset(x,rank(x$md.error)>length(x$nj)-nlabels),
                          nudge_x=1.5,nudge_y = .5,size=3)
      } else {
        plotout<-ggplot(x,aes_string("index","md.b")) +
          geom_point(shape=1) + ylab("R.E. distance")+
          geom_text_repel(aes_string(label="ind"),data=subset(x,rank(x$md.b)>length(x$nj)-nlabels),
                          nudge_x=1.5,nudge_y = 0,size=3)
      }
    } else {
      njvec <- sort(unique(x$nj))
      if (type=="error") {
        plotout<-ggplot(x,aes_string("nj","md.error")) +
          geom_point(position = position_jitter(width = .3,height = 0),
                     shape=1) + ylab("error distance")+xlab("number of observations")+
          geom_text(aes_string(label="ind"),data=subset(x,rank(x$md.error)>length(x$nj)-nlabels),
                    nudge_x=0,nudge_y = .5,size=3)+
          scale_x_continuous(breaks=njvec)
      } else {
        plotout<-ggplot(x,aes_string("nj","md.b")) +
          geom_point(position = position_jitter(width = .3,height = 0),
                     shape=1) + ylab("R.E. distance")+xlab("number of observations")+
          geom_text(aes_string(label="ind"),data=subset(x,rank(x$md.b)>length(x$nj)-nlabels),
                    nudge_x=0,nudge_y = .5,size=3)+
          scale_x_continuous(breaks=njvec)
      }
    }
  } else {
    if (n_distinct(x$nj) ==1) {
      nj1 <- x$nj[1]
      if (distr=="sn"|distr=="norm") mdquantile <- qchisq(level,nj1)
      if (distr=="st"|distr=="t") mdquantile <- nj1*qf(level,nj1,nu)
      if (distr=="ssl"|distr=="sl") {
        mdquantile <- uniroot(function(x.) pMDsl(x.,nu,nj1) -level,
                              c(1,qchisq(level,nj1)^2))$root
        #y.<-seq(1,qchisq(level,nj1)^2,length.out =1000)
        #py.<-pchisq(y.,nj1)-2^nu*gamma(nj1/2+nu)/(y.^nu)/gamma(nj1/2)*
        #                      pchisq(y.,nj1+2*nu)
        #mdquantile <- min(y.[py.>level])
      }
      if (distr=="scn"|distr=="cn") {
        mdquantile <- uniroot(function(x.) pMDcn(x.,nu,nj1) -level,
                              c(1,qchisq(level,nj1)^2))$root
      }

      plotout<-ggplot(x,aes_string("index","md")) +
        geom_point(shape=1) + ylab("Mahalanobis distance")+
        geom_text_repel(aes_string(label="ind"),data=subset(x,rank(x$md)>length(x$nj)-nlabels),
                        nudge_x=1.5,nudge_y = .5,size=3) +
        geom_hline(yintercept = mdquantile,col=4,linetype="dashed")
      attr(plotout,"info") <- data.frame(nj= nj1,quantile=c(mdquantile))
    } else {
      njvec <- sort(unique(x$nj))
      if (distr=="sn"|distr=="norm") mdquantile <- qchisq(level,njvec)
      if (distr=="st"|distr=="t") mdquantile <- njvec*qf(level,njvec,nu)
      if (distr=="ssl"|distr=="sl") {
        mdquantile <- unlist(lapply(as.list(njvec),function(nj1) uniroot(function(x) pMDsl(x,nu,nj1) -level,
                                                                         c(1,qchisq(level,nj1)^2))$root))
      }
      if (distr=="scn"|distr=="cn") { #nu<-c(.25,.3)
        mdquantile <- unlist(lapply(as.list(njvec),function(nj1) uniroot(function(x) pMDcn(x,nu,nj1) -level,
                                                                         c(1,qchisq(level,nj1)^2))$root))
      }
      datline <- data.frame(nj= c(njvec,max(njvec)+1),quantile=c(mdquantile,max(mdquantile)))
      datline$nj2<-datline$nj-.5
      plotout<-ggplot(x,aes_string("nj","md")) +
        geom_point(position = position_jitter(width = .25,height = 0),
                   shape=1) + ylab("Mahalanobis distance")+ xlab("number of observations")+
        geom_text(aes_string(label="ind"),data=subset(x,rank(x$md)>length(x$nj)-nlabels),
                  nudge_x=0,nudge_y = 0,size=3) +
        geom_step(aes_string(x="nj2",y="quantile"),data=datline,color=4,linetype="dashed")+
        scale_x_continuous(breaks=njvec)
      attr(plotout,"info") <- data.frame(nj= c(njvec),quantile=c(mdquantile))
    }
  }
  plotout + theme_minimal()
}
pMDsl <- function(x,nu,nj) {
    pchisq(x,nj)- 2^nu*gamma(nj/2+nu)/(x^nu)/gamma(nj/2)*
                  pchisq(x,nj+2*nu)
  }
pMDcn <- function(x,nu,nj) {
  nu[1]*pchisq(x*nu[2],nj) + (1-nu[1])*pchisq(x,nj)
}

# Healy's plot
gerar_smsn_healy = function(jvec,x,z,sigma2,Dsqrti,beta1,lambda,distr,nu,
                            ind,time,depStruct,phi) {
  if (distr=="sn"|distr=="norm") {ui=1; c.=-sqrt(2/pi)}
  if (distr=="st"|distr=="t") {ui=rgamma(1,nu/2,nu/2); c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ssl"|distr=="sl") {ui=rbeta(1,nu,1); c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn"|distr=="cn") {ui=ifelse(runif(1)<nu[1],nu[2],1);
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
  return(data.frame(y=Yi,ind=ind[jvec],time=time[jvec]))
}

healy.plot <- function(object,dataPlus=NULL,dotsize=0.4,calcCI = FALSE,
                       IClevel, MCiter, seed,...) {
  if (!missing(seed)) set.seed(seed)
  if(!inherits(object,c("SMSN","SMN"))) stop("object must inherit from class SMSN or SMN")
  if (!is.null(dataPlus)) {
    data <- dataPlus
  } else data <- object$data
  groupVar<-object$groupVar
  ind <-data[,groupVar]
  data <- data[order(data[,object$groupVar]),]
  njvec <- tapply(ind,ind,length)
  nj1 <- as.numeric(njvec[1])
  if (!all(njvec==nj1)) stop("for this plot all subjects must have the same number of observations, you can predict the missing values and enter the full dataset using the dataPlus argument")
  #
  distr <- object$distr
  nu <- object$estimates$nu
  depStruct <- object$depStruct
  #
  if (calcCI) {
    if (!missing(IClevel)) {
      if (IClevel>=1|IClevel<=0) stop("0<IClevel<1 needed")
    } else IClevel <- .95
    if (!missing(MCiter)) {
      if (!is.wholenumber(MCiter)) stop("MCiter must be an integer positive number")
      if (MCiter<=0) stop("MCiter must be an integer positive number")
    } else MCiter<- 300
  }
  #
  mahaldist<- sort(mahalDist(object,dataPlus = data))
  #
  m <- length(mahaldist)
  empprob <- (1:m)/(m)
  #
  if (distr=="sn"|distr=="norm") theoprob <- pchisq(mahaldist,nj1)
  if (distr=="st"|distr=="t") theoprob <- pf(mahaldist/nj1,nj1,nu)
  if (distr=="ssl"|distr=="sl") {
    theoprob <- pchisq(mahaldist,nj1) - 2^nu*gamma(nj1/2+nu)/(mahaldist^nu)/gamma(nj1/2)*
      pchisq(mahaldist,nj1+2*nu)
  }
  if (distr=="scn"|distr=="cn") {
    theoprob <- nu[1]*pchisq(nu[2]*mahaldist,nj1) + (1-nu[1])*pchisq(mahaldist,nj1)
  }
  ####
  #print(cor(empprob,theoprob))
  distrp <- toupper(distr)
  if (distrp=='NORM') distrp<-"N"
  depStructp <- depStruct
  if (depStruct=="ARp") depStructp <- paste0("AR(",length(object$estimates$phi),")")
  p1<-qplot(empprob,theoprob,size=I(dotsize),...=...) + geom_abline(intercept=0,slope=1) +
    ylab(NULL) + xlab(NULL) + theme_minimal() + ggtitle(paste0(depStructp,'-',distrp,'-LMM')) +
    theme(plot.title = element_text( face="italic", size=10))
  #
  if (calcCI){
    timeVar <- object$timeVar
    if (!is.null(timeVar)) {
      time<-data[,timeVar]
    } else{
      time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
    }
    x <- model.matrix(object$formula$formFixed,data=data)
    z<-model.matrix(object$formula$formRandom,data=data)
    p<-ncol(x)
    q1<-ncol(z)
    #
    phi <- object$estimates$phi
    Dfit2 <- object$estimates$D #Dfit%*%Dfit
    Dfit <- matrix.sqrt(Dfit2)
    sigma2<- object$estimates$sigma2
    lambda <- object$estimates$lambda
    if (is.null(lambda)) lambda<- rep(0,nrow(Dfit))
    beta<- object$estimates$beta
    N<- nrow(data)
    ind_levels <- levels(ind)
    vars <- unique(c(all.vars(object$formula$formFixed)[-1],
                     all.vars(object$formula$formRandom),
                     groupVar,timeVar))
    mahalSim <- matrix(ncol=length(ind_levels),nrow=MCiter)
    for (sample in 1:MCiter){
      dadosi <- tapply(1:N,ind,gerar_smsn_healy,x=x,z=z,sigma2=sigma2,Dsqrti=Dfit,
                       beta1=beta,lambda=lambda,distr=distr,nu=nu,ind=ind,time=time,
                       depStruct=depStruct,phi=phi) %>% bind_rows()
      names(dadosi)[1] <- all.vars(object$formula$formFixed)[1]
      dadosi[,vars] <- data[,vars]
      mahaldisti<-sort(mahalDist(object,dataPlus = dadosi))
      #
      if (distr=="sn"||distr=="norm") theoprobi <- pchisq(mahaldisti,nj1)
      if (distr=="st"||distr=="t") theoprobi <- pf(mahaldisti/nj1,nj1,nu)
      if (distr=="ssl"||distr=="sl") {
        theoprobi <- pchisq(mahaldisti,nj1) - 2^nu*gamma(nj1/2+nu)/(mahaldisti^nu)/gamma(nj1/2)*
          pchisq(mahaldisti,nj1+2*nu)
      }
      if (distr=="scn"||distr=="cn") {
        theoprobi <- nu[1]*pchisq(nu[2]*mahaldisti,nj1) + (1-nu[1])*pchisq(mahaldisti,nj1)
      }
      mahalSim[sample,]<-theoprobi
    }
    MCmahalIC<-apply(mahalSim,2,function(x) quantile(x,probs = c((1-IClevel)/2,(1+IClevel)/2)))
    MCmahalIC <-t(MCmahalIC) %>% as.data.frame()
    colnames(MCmahalIC) <- c('inf','sup')
    p1 <- p1+geom_line(aes(y=.data$inf),data = MCmahalIC,linetype=2,color=4)+
      geom_line(aes(y=.data$sup),data = MCmahalIC,linetype=2,color=4)
  }
  p1
}
