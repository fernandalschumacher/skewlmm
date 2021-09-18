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
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
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
      Dest <- object$estimates$D
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
      Dest <-object$estimates$D
      delta = object$estimates$lambda/as.numeric(
        sqrt(1+t(object$estimates$lambda)%*%(object$estimates$lambda)))
      Delta = matrix.sqrt(object$estimates$D)%*%delta
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
                     calcCI=FALSE,levelCI,MCiter,seed){
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
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  }
  if (any(!is.wholenumber(time))||any(time<=0)) {
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
    warning("for computing ACF, time must contain positive integer numbers \ntimeVar was ignored and the order was be used instead")
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
      stop("CIs are obtained through an empirical method that is only
      appropriate for marginal modified/normalized residuals, please use
      (resType='normalized' or resType='modified') and resLevel='marginal',
      or calcCI=FALSE")
    if (!missing(levelCI)) {
      if (levelCI>=1|levelCI<=0) stop("0<levelCI<1 needed")
    } else levelCI <- .95
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
      MCacfIC<-apply(acfMC,2,function(x) quantile(x,probs = c((1-levelCI)/2,(1+levelCI)/2)))
      row.names(MCacfIC) <- NULL
      out <- cbind(out,CI = rbind(rep(NA,2),t(MCacfIC)))
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
      MCacfIC<-apply(acfMC,2,function(x) quantile(x,probs = c((1-levelCI)/2,(1+levelCI)/2)))
      row.names(MCacfIC)<-NULL
      out <- cbind(out,CI = rbind(rep(NA,2),t(MCacfIC)))
    }
  }
  class(out) <- c("acfresid","data.frame")
  attr(out,'distr') <- object$distr
  attr(out,'depStruct') <- object$depStruct
  if (object$depStruct=="ARp") attr(out,'pAR') <- length(object$estimates$phi)
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
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
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

#plot ACF
# plot.ACF <- function(ACFobj,...) {
#   plot(ACFobj$lag,ACFobj$ACF,type="h",...,xlab="lag",ylab="ACF");abline(h=0)
#   if (ncol(ACFobj)>3) {
#     lagPlus <- c(ACFobj$lag,max(ACFobj$lag)+1)
#     ICinf <- c(ACFobj[,4],0)
#     ICsup <- c(ACFobj[,5],0)
#     lines(lagPlus-.5,ICinf,lty=2,col=4,type="s")
#     lines(lagPlus-.5,ICsup,lty=2,col=4,type="s")
#   }
# }
plot.acfresid <-  function(x,...) {
  distrp <- toupper(attr(x,'distr'))
  if (distrp=='NORM') distrp<-"N"
  depStructp <- attr(x,'depStruct')
  if (depStructp=="ARp") depStructp <- paste0("AR(",attr(x,'pAR'),")")
  if (ncol(x)<=3) {
    ggplot(data = x,aes_string(x = "lag",y="ACF"))+
      theme_minimal()+geom_hline(aes(yintercept = 0))+
      geom_segment(mapping = aes(xend = lag, yend = 0))+ ylim(c(-1,1))+
      ggtitle(paste0(depStructp,'-',distrp,'-LMM')) +
      theme(plot.title = element_text( face="italic", size=10))
  } else {
    lagmax <- max(x$lag)
    datIC <- rbind(x[,c(1,4,5)],
                   c(lagmax+1,min(x$`CI.1`,na.rm=T),max(x$`CI.2`,na.rm=T)))
    datIC$lag <- datIC$lag-.5
    names(datIC) <- c("lag","inf","sup")
    ggplot(x,aes_string(x = "lag",y="ACF")) +
    #ggplot(x,aes(x = lag,y=ACF)) +
    geom_hline(aes(yintercept = 0)) +
      coord_cartesian(xlim=c(0,max(x$lag)))+
      geom_segment(mapping = aes_string(xend = "lag", yend = 0)) +
      geom_step(aes_string(x="lag",y="inf"),data=datIC,color=4,linetype="dashed",na.rm=TRUE)+
      geom_step(aes_string(x="lag",y="sup"),data=datIC,color=4,linetype="dashed",na.rm=TRUE)+
      theme_minimal()+ ylim(c(-1,1))+
      ggtitle(paste0(depStructp,'-',distrp,'-LMM')) +
      theme(plot.title = element_text( face="italic", size=10))
      }
}

#plot obj
plot.SMSN <- function(x,type="response",level="conditional",useweight=TRUE,alpha=.3,...) {
  resid <- residuals(x,type=type,level=level)
  distrp <- toupper(x$distr)
  if (distrp=='NORM') distrp<-"N"
  depStructp <- x$depStruct
  if (depStructp=="ARp") depStructp <- paste0("AR(",length(x$estimates$phi),")")
  if (useweight){
    peso <-data.frame(weight=x$uhat)
    peso$ind <- row.names(peso)
    peso <- left_join(x$data,peso,by='ind')
    peso$fitted <- fitted(x)
    peso$resid <- resid
    ggplot(peso, aes_string(x="fitted",y="resid",color="weight"))+geom_point()+theme_minimal() +
      geom_hline(yintercept = 0,linetype="dashed") + ylab(attr(resid,"label")) +
      xlab("fitted values") +
      scale_color_continuous(high = "#132B43", low = "#56B1F7") +
      ggtitle(paste0(depStructp,'-',distrp,'-LMM')) +
      theme(plot.title = element_text( face="italic", size=10))
  } else{
    qplot(fitted(x),resid,alpha=I(alpha),...=...) +theme_minimal() +
      geom_hline(yintercept = 0,linetype="dashed") + ylab(attr(resid,"label")) +
      xlab("fitted values")+ ggtitle(paste0(depStructp,'-',distrp,'-LMM')) +
      theme(plot.title = element_text( face="italic", size=10))
  }
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
    time <- numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  }
  distr <- fitobject$distr
  nu <- fitobject$estimates$nu
  #
  x$nj <- tapply(time,ind,length)
  x$index <- seq_along(x$nj)
  x$ind <- levels(ind)
  #
  distrp <- toupper(fitobject$distr)
  if (distrp=='NORM') distrp<-"N"
  depStructp <- fitobject$depStruct
  if (depStructp=="ARp") depStructp <- paste0("AR(",length(fitobject$estimates$phi),")")
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
      if (distr=="sn"||distr=="norm") mdquantile <- qchisq(level,nj1)
      if (distr=="st"||distr=="t") mdquantile <- nj1*qf(level,nj1,nu)
      if (distr=="ssl"||distr=="sl") {
        mdquantile <- try(uniroot(function(x.) pMDsl(x.,nu,nj1) -level,
                              c(1,qchisq(level,nj1)^2))$root,silent = T)
        if (class(mdquantile)=='try-error') {
          mdquantile <- try(uniroot(function(x.) pMDsl(x.,nu,nj1) -level,
                                    c(1,qchisq(level,nj1)^5))$root,silent = T)
          if (class(mdquantile)=='try-error') mdquantile<-NULL
        }
        #y.<-seq(1,qchisq(level,nj1)^2,length.out =1000)
        #py.<-pchisq(y.,nj1)-2^nu*gamma(nj1/2+nu)/(y.^nu)/gamma(nj1/2)*
        #                      pchisq(y.,nj1+2*nu)
        #mdquantile <- min(y.[py.>level])
      }
      if (distr=="scn"||distr=="cn") {
        mdquantile <- try(uniroot(function(x.) pMDcn(x.,nu,nj1) -level,
                              c(1,qchisq(level,nj1)^2))$root,silent = T)
        if (class(mdquantile)=='try-error') {
          mdquantile <- try(uniroot(function(x.) pMDcn(x.,nu,nj1) -level,
                                    c(1,qchisq(level,nj1)^5))$root,silent = T)
          if (class(mdquantile)=='try-error') mdquantile<-NULL
        }
      }

      plotout<-ggplot(x,aes_string("index","md")) +
        geom_point(shape=1) + ylab("Mahalanobis distance")+
        geom_text_repel(aes_string(label="ind"),data=subset(x,rank(x$md)>length(x$nj)-nlabels),
                  nudge_x=1.5,nudge_y = .5,size=3) +
        geom_hline(yintercept = mdquantile,col=4,linetype="dashed")
      attr(plotout,"info") <- data.frame(nj= nj1,quantile=c(mdquantile))
    } else {
      njvec <- sort(unique(x$nj))
      if (distr=="sn"||distr=="norm") mdquantile <- qchisq(level,njvec)
      if (distr=="st"||distr=="t") mdquantile <- njvec*qf(level,njvec,nu)
      if (distr=="ssl"||distr=="sl") {
        mdquantile <- try(unlist(lapply(as.list(njvec),function(nj1) uniroot(function(x) pMDsl(x,nu,nj1) -level,
                                    c(1,qchisq(level,nj1)^2))$root)),silent = T)
        if (class(mdquantile)=='try-error') {
          mdquantile <- try(unlist(lapply(as.list(njvec),function(nj1) uniroot(function(x) pMDsl(x,nu,nj1) -level,
                                  c(1,qchisq(level,nj1)^5))$root)),silent = T)
          if (class(mdquantile)=='try-error') mdquantile<-NULL
        }
      }
      if (distr=="scn"||distr=="cn") { #nu<-c(.25,.3)
        mdquantile <- try(unlist(lapply(as.list(njvec),function(nj1) uniroot(function(x) pMDcn(x,nu,nj1) -level,
                                       c(1,qchisq(level,nj1)^2))$root)),silent = T)
        if (class(mdquantile)=='try-error') {
          mdquantile <- try(unlist(lapply(as.list(njvec),function(nj1) uniroot(function(x) pMDcn(x,nu,nj1) -level,
                                   c(1,qchisq(level,nj1)^5))$root)),silent = T)
          if (class(mdquantile)=='try-error') mdquantile<-NULL
        }
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
  plotout + theme_minimal()+ ggtitle(paste0(depStructp,'-',distrp,'-LMM')) +
    theme(plot.title = element_text( face="italic", size=10))
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
  return(data.frame(y=Yi,ind=ind[jvec],time=time[jvec]))
}

healy.plot <- function(object,dataPlus=NULL,dotsize=0.4,calcCI = FALSE,
                       levelCI, MCiter, seed,...) {
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
    if (!missing(levelCI)) {
      if (levelCI>=1|levelCI<=0) stop("0<levelCI<1 needed")
    } else levelCI <- .95
    if (!missing(MCiter)) {
      if (!is.wholenumber(MCiter)) stop("MCiter must be an integer positive number")
      if (MCiter<=0) stop("MCiter must be an integer positive number")
    } else MCiter<- 300
  }
  #
  mahaldist<- sort(mahalDist(object,dataPlus = data))
  #
  m <- length(mahaldist)
  empprob <- (1:m)/(m)#seq(0,1,len=m)#
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
  p1<-qplot(empprob,theoprob,size=I(dotsize),...=...) + #geom_abline(intercept=0,slope=1) +
    ylab(NULL) + xlab(NULL) + theme_minimal() + ggtitle(paste0(depStructp,'-',distrp,'-LMM')) +
    theme(plot.title = element_text( face="italic", size=10))
  #
  if (calcCI){
    timeVar <- object$timeVar
    if (!is.null(timeVar)) {
      time<-data[,timeVar]
    } else{
      time <- numeric(length = length(ind))
      for (indi in levels(ind)) time[ind==indi] <- seq_len(sum(ind==indi))
      #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
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
    #MCmahalIC<-apply(mahalSim,2,function(x) quantile(x,probs = c((1-levelCI)/2,(1+levelCI)/2)))
    MCmahalIC<-apply(mahalSim,2,function(x) quantile(x,probs = c((1-levelCI)/2,.5,(1+levelCI)/2)))
    MCmahalIC <-t(MCmahalIC) %>% as.data.frame()
    colnames(MCmahalIC) <- c('inf','med','sup')
    p1 <- p1+geom_line(aes(y=.data$inf),data = MCmahalIC,linetype=2,color=4)+
      geom_line(aes(y=.data$sup),data = MCmahalIC,linetype=2,color=4)+
      geom_line(aes(y=.data$med),data = MCmahalIC,linetype=2,color=1)
  } else {
    p1 <- p1+geom_abline(intercept=0,slope=1)
  }
  p1
}

## bootstrap

gerar_smn_healy <- function(jvec,x,z,sigma2,Dsqrti,beta1,distr,nu,
                            ind,time,depStruct,phi) {
  if (distr=="sn"||distr=="norm") {ui=1}
  if (distr=="st"||distr=="t") {ui=rgamma(1,nu/2,nu/2)}
  if (distr=="ssl"||distr=="sl") {ui=rbeta(1,nu,1)}
  if (distr=="scn"||distr=="cn") {ui=ifelse(runif(1)<nu[1],nu[2],1)}
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jvec,  ],ncol=p)
  z1=matrix(z[jvec,  ],ncol=q1)
  nj = nrow(x1)
  Sig = errorVar(time[jvec],sigma2=sigma2,depStruct=depStruct,phi=phi)#sigma2*diag(nj)
  Gammab = Dsqrti%*%Dsqrti
  Beta = matrix(beta1,ncol=1)
  bi = t(rmvnorm(1,sigma=ui^(-1)*Gammab))
  Yi = t(rmvnorm(1,x1%*%Beta+z1%*%bi,sigma=ui^(-1)*Sig))
  return(data.frame(y=Yi,ind=ind[jvec],time=time[jvec]))
}

gen_fit <- function(fit1, ...) {
  x <- model.matrix(fit1$formula$formFixed,data=fit1$data)
  #varsx <- all.vars(fit1$formula$formFixed)[-1]
  y <-fit1$data[,all.vars(fit1$formula$formFixed)[1]]
  z<-model.matrix(fit1$formula$formRandom,data=fit1$data)
  ind <-fit1$data[,fit1$groupVar]
  #fit1$data$ind <- ind
  time <- fit1$data$time

  m<-n_distinct(ind)
  N<-length(ind)
  p<-ncol(x)
  q1<-ncol(z)
  distr <- fit1$distr
  if (distr=='ssl') distr='ss'
  if (is.null(fit1$estimates$lambda)) {
    #fit1$estimates$lambda=rep(0,q1)
    if (distr=="norm") distrs="sn"
    if (distr=="t") distrs="st"
    if (distr=="sl") distrs="ss"
    if (distr=="cn") distrs="scn"
    sym <- TRUE
  } else sym <- FALSE

  vars <- unique(c(all.vars(fit1$formula$formFixed)[-1],
                   all.vars(fit1$formula$formRandom),
                   #fit1$groupVar,fit1$timeVar,
                   "time","ind"))
  #
  if (sym) {
    dadosi = tapply(seq_len(N),ind,gerar_smn_healy,x=x,z=z,sigma2=fit1$estimates$sigma2,
                    Dsqrti=Dmatrix(fit1$estimates$dsqrt),beta1=fit1$estimates$beta,
                    distr=distr,nu=fit1$estimates$nu,ind=ind,time=time,
                    depStruct=fit1$depStruct,phi=fit1$estimates$phi) %>% bind_rows()
    names(dadosi)[1] <- all.vars(fit1$formula$formFixed)[1]
    #dadosi <- dadosi[order(dadosi$ind)]
    dadosi <- left_join(dadosi,fit1$data[,vars],by=c('ind','time'))
    #dadosi[,vars] <- fit1$data[,vars]
    #
    # if (!use_as_start) {
    #   lmefit = try(lme(fit1$formula$formFixed,
    #                    random=formula(paste('~',as.character(fit1$formula$formRandom)[length(fit1$formula$formRandom)],
    #                                         '|',"ind")),data=dadosi),silent=T)
    #   if (class(lmefit)=="try-error") {
    #     lmefit = try(lme(fit1$formula$formFixed,random=~1|ind,data=dadosi),silent=TRUE)
    #     if (class(lmefit)=="try-error") {
    #       stop("error in calculating initial values")
    #     } else {
    #       #lambdainit <- rep(1,q1)*sign(as.numeric(skewness(random.effects(lmefit))))
    #       D1init <- diag(q1)*as.numeric(var(random.effects(lmefit)))
    #     }
    #   } else {
    #     #lambdainit <- sign(as.numeric(skewness(random.effects(lmefit))))*1
    #     D1init <- (var(random.effects(lmefit)))
    #   }
    #   beta1 <- as.numeric(lmefit$coefficients$fixed)
    #   sigmae <- as.numeric(lmefit$sigma^2)
    #   D1 <- D1init
    #   #lambda <- lambdainit*fit1$skewind
    # } else{
      beta1 <- fit1$estimates$beta
      sigmae <- fit1$estimates$sigma2
      D1 <- fit1$estimates$D
    #}
    #
    if (fit1$depStruct=="UNC") {
      obj.out <- try(DAAREM.UNC(formFixed = fit1$formula$formFixed,
                                formRandom = fit1$formula$formRandom,
                                data = dadosi, groupVar = 'ind',
                                distr = distrs, beta1 = beta1,
                                sigmae = sigmae, D1 = D1,
                                nu = fit1$estimates$nu, lb=fit1$control$lb, lu=fit1$control$lu,
                                diagD=fit1$diagD, precisao=fit1$control$tol, informa = FALSE,
                                max.iter=fit1$control$max.iter, showiter = FALSE,
                                showerroriter = FALSE,
                                algorithm=fit1$control$algorithm,
                                control.daarem = fit1$control$control.daarem,
                                parallelnu=FALSE, ncores=NULL),silent = TRUE)
    }  else if (fit1$depStruct=="ARp") {
      obj.out <- try(DAAREM.AR(formFixed = fit1$formula$formFixed,
                               formRandom = fit1$formula$formRandom,
                               data = dadosi, groupVar = 'ind',
                               pAR=length(fit1$estimates$phi),timeVar = "time",
                               distr = distrs, beta1 = beta1,
                               sigmae = sigmae, phiAR = fit1$estimates$phi,
                               D1 = D1,  nu = fit1$estimates$nu,
                               lb=fit1$control$lb, lu=fit1$control$lu,
                               diagD=fit1$diagD,
                               precisao=fit1$control$tol, informa = FALSE,
                               max.iter=fit1$control$max.iter, showiter = FALSE,
                               showerroriter = FALSE,
                               algorithm=fit1$control$algorithm,
                               control.daarem = fit1$control$control.daarem,
                               parallelphi = FALSE,
                               parallelnu=FALSE, ncores=NULL),silent = TRUE)
    } else if (fit1$depStruct=="CS") {
      obj.out <-try(DAAREM.CS(formFixed = fit1$formula$formFixed,
                              formRandom = fit1$formula$formRandom,
                              data = dadosi, groupVar = 'ind',
                              distr = distrs, beta1 = beta1,
                              sigmae = sigmae, phiCS = fit1$estimates$phi,
                              D1 = D1, nu = fit1$estimates$nu,
                              lb=fit1$control$lb, lu=fit1$control$lu,
                              diagD=fit1$diagD,
                              precisao=fit1$control$tol, informa = FALSE,
                              max.iter=fit1$control$max.iter, showiter = FALSE,
                              showerroriter = FALSE,
                              algorithm=fit1$control$algorithm, parallelphi = FALSE,
                              control.daarem = fit1$control$control.daarem,
                              parallelnu=FALSE, ncores=NULL),silent = TRUE)
    } else if (fit1$depStruct=="DEC") {
      obj.out <-try(DAAREM.DEC(formFixed = fit1$formula$formFixed,
                               formRandom = fit1$formula$formRandom,
                               data = dadosi, groupVar = 'ind',
                               timeVar = "time", distr = distrs,
                               beta1 = beta1,sigmae = sigmae,
                               parDEC = fit1$estimates$phi,D1 = D1,
                               nu = fit1$estimates$nu, lb=fit1$control$lb, lu=fit1$control$lu,
                               luDEC=fit1$control$luDEC, diagD=fit1$diagD,
                               precisao=fit1$control$tol, informa = FALSE,
                               max.iter=fit1$control$max.iter, showiter = FALSE,
                               showerroriter = FALSE,
                               algorithm=fit1$control$algorithm,
                               control.daarem = fit1$control$control.daarem,
                               parallelphi = FALSE,
                               parallelnu=FALSE, ncores=NULL),silent = TRUE)
    } else if (fit1$depStruct=="CAR1") {
      obj.out <-try(DAAREM.CAR1(formFixed = fit1$formula$formFixed,
                                formRandom = fit1$formula$formRandom,
                                data = dadosi, groupVar = 'ind',
                                timeVar = "time", distr = distrs,
                                beta1 = beta1,sigmae = sigmae,
                                phiCAR1 = fit1$estimates$phi,D1 = D1,
                                nu = fit1$estimates$nu,lb=fit1$control$lb, lu=fit1$control$lu,
                                diagD=fit1$diagD, precisao=fit1$control$tol, informa = FALSE,
                                max.iter=fit1$control$max.iter, showiter = FALSE,
                                showerroriter = FALSE,
                                algorithm=fit1$control$algorithm,
                                control.daarem = fit1$control$control.daarem,
                                parallelphi = FALSE,
                                parallelnu=FALSE, ncores=NULL),silent = TRUE)
    }
  } else{
    dadosi = tapply(seq_len(N),ind,gerar_smsn_healy,x=x,z=z,sigma2=fit1$estimates$sigma2,
                    Dsqrti=Dmatrix(fit1$estimates$dsqrt),beta1=fit1$estimates$beta,
                    lambda=as.matrix(fit1$estimates$lambda),distr=distr,
                    nu=fit1$estimates$nu,ind=ind,time=time,depStruct=fit1$depStruct,
                    phi=fit1$estimates$phi) %>% bind_rows()
    names(dadosi)[1] <- all.vars(fit1$formula$formFixed)[1]
    dadosi[,vars] <- fit1$data[,vars]
    #
    # if (!use_as_start) {
    #   lmefit = try(lme(fit1$formula$formFixed,
    #                    random=formula(paste('~',as.character(fit1$formula$formFixed)[length(fit1$formula$formFixed)],
    #                                         '|',"ind")),data=dadosi),silent=T)
    #   if (class(lmefit)=="try-error") {
    #     lmefit = try(lme(fit1$formula$formFixed,random=~1|ind,data=dadosi),silent=TRUE)
    #     if (class(lmefit)=="try-error") {
    #       stop("error in calculating initial values")
    #     } else {
    #       lambdainit <- rep(1,q1)*sign(as.numeric(skewness(random.effects(lmefit))))
    #       D1init <- diag(q1)*as.numeric(var(random.effects(lmefit)))
    #     }
    #   } else {
    #     lambdainit <- sign(as.numeric(skewness(random.effects(lmefit))))*1
    #     D1init <- (var(random.effects(lmefit)))
    #   }
    #   beta1 <- as.numeric(lmefit$coefficients$fixed)
    #   sigmae <- as.numeric(lmefit$sigma^2)
    #   D1 <- D1init
    #   lambda <- lambdainit*fit1$skewind
    # } else{
      beta1 <- fit1$estimates$beta
      sigmae <- fit1$estimates$sigma2
      D1 <- fit1$estimates$D
      lambda <- fit1$estimates$lambda
    #}
    #
    if (fit1$depStruct=="UNC") {
      obj.out <- try(DAAREM.SkewUNC(formFixed = fit1$formula$formFixed,
                                    formRandom = fit1$formula$formRandom,
                                    data = dadosi, groupVar = 'ind',
                                    distr = distr, beta1 = beta1,
                                    sigmae = sigmae, D1 = D1,lambda = lambda,
                                    nu = fit1$estimates$nu,
                                    lb=fit1$control$lb, lu=fit1$control$lu,
                                    skewind = fit1$skewind, diagD=fit1$diagD,
                                    precisao=fit1$control$tol, informa = FALSE,
                                    max.iter=fit1$control$max.iter, showiter = FALSE,
                                    showerroriter = FALSE,
                                    algorithm=fit1$control$algorithm,
                                    control.daarem = fit1$control$control.daarem,
                                    parallelnu=FALSE, ncores=NULL),silent = TRUE)
    } else if (fit1$depStruct=="ARp") {
      obj.out <- try(DAAREM.SkewAR(formFixed = fit1$formula$formFixed,
                                   formRandom = fit1$formula$formRandom,
                                   data = dadosi, groupVar = 'ind',
                                   pAR=length(fit1$estimates$phi),timeVar = "time",
                                   distr = distr, beta1 = beta1,
                                   sigmae = sigmae, phiAR = fit1$estimates$phi,
                                   D1 = D1,lambda = lambda,
                                   nu = fit1$estimates$nu,
                                   lb=fit1$control$lb, lu=fit1$control$lu,
                                   skewind = fit1$skewind, diagD=fit1$diagD,
                                   precisao=fit1$control$tol, informa = FALSE,
                                   max.iter=fit1$control$max.iter, showiter = FALSE,
                                   showerroriter = FALSE,
                                   algorithm=fit1$control$algorithm,
                                   control.daarem = fit1$control$control.daarem,
                                   parallelphi = FALSE,
                                   parallelnu=FALSE, ncores=NULL),silent = TRUE)
    } else if (fit1$depStruct=="CS") {
      obj.out <-try(DAAREM.SkewCS(formFixed = fit1$formula$formFixed,
                                  formRandom = fit1$formula$formRandom,
                                  data = dadosi, groupVar = 'ind',
                                  distr = distr, beta1 = beta1,
                                  sigmae = sigmae, phiCS = fit1$estimates$phi,
                                  D1 = D1, lambda = lambda,
                                  nu = fit1$estimates$nu,
                                  lb=fit1$control$lb, lu=fit1$control$lu,
                                  skewind = fit1$skewind, diagD=fit1$diagD,
                                  precisao=fit1$control$tol, informa = FALSE,
                                  max.iter=fit1$control$max.iter, showiter = FALSE,
                                  showerroriter = FALSE,
                                  algorithm=fit1$control$algorithm,
                                  control.daarem = fit1$control$control.daarem,
                                  parallelphi = FALSE,
                                  parallelnu=FALSE, ncores=NULL),silent = TRUE)
    } else if (fit1$depStruct=="DEC") {
      obj.out <-try(DAAREM.SkewDEC(formFixed = fit1$formula$formFixed,
                                   formRandom = fit1$formula$formRandom,
                                   data = dadosi, groupVar = 'ind',
                                   timeVar = "time", distr = distr,
                                   beta1 = beta1,sigmae = sigmae,
                                   parDEC = fit1$estimates$phi,D1 = D1,
                                   lambda = lambda, nu = fit1$estimates$nu,
                                   lb=fit1$control$lb, lu=fit1$control$lu,
                                   luDEC=fit1$control$luDEC,
                                   skewind = fit1$skewind, diagD=fit1$diagD,
                                   precisao=fit1$control$tol, informa = FALSE,
                                   max.iter=fit1$control$max.iter, showiter = FALSE,
                                   showerroriter = FALSE,
                                   algorithm=fit1$control$algorithm,
                                   control.daarem = fit1$control$control.daarem,
                                   parallelphi = FALSE,
                                   parallelnu=FALSE, ncores=NULL),silent = TRUE)
    } else if (fit1$depStruct=="CAR1") {
      obj.out <-try(DAAREM.SkewCAR1(formFixed = fit1$formula$formFixed,
                                    formRandom = fit1$formula$formRandom,
                                    data = dadosi, groupVar = 'ind',
                                    timeVar = "time", distr = distr,
                                    beta1 = beta1,sigmae = sigmae,
                                    phiCAR1 = fit1$estimates$phi,D1 = D1,
                                    lambda = lambda, nu = fit1$estimates$nu,
                                    lb=fit1$control$lb, lu=fit1$control$lu,
                                    skewind = fit1$skewind, diagD=fit1$diagD,
                                    precisao=fit1$control$tol, informa = FALSE,
                                    max.iter=fit1$control$max.iter, showiter = FALSE,
                                    showerroriter = FALSE,
                                    algorithm=fit1$control$algorithm,
                                    control.daarem = fit1$control$control.daarem,
                                    parallelphi = FALSE,
                                    parallelnu=FALSE, ncores=NULL),silent = TRUE)
    }
  }
  if (class(obj.out)[1]!='try-error') {
    return(obj.out$theta)
  } else return(NULL)
}

boot_par <- function(object, B=100, seed = 123) {#method
  if (!(inherits(object, "SMN")||inherits(object, "SMSN"))) stop("object must inherit from class SMSN or SMN")
  if (!is.null(object$timeVar)) {
    object$data$time<-object$data[,object$timeVar]
  } else{
    object$data$time <- numeric(length = length(object$data$ind))
    for (indi in levels(object$data$ind)) {
      object$data$time[object$data$ind==indi] <- seq_len(sum(object$data$ind==indi))
    }
    #time<- flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
  }
  #object$data$ind<-object$data[,object$groupVar]
  plan(multisession)
  a1<-suppressMessages(future_map(.x = seq_len(B),.f = gen_fit,
                  fit1=object,.options = furrr_options(seed = seed)))
  plan(sequential)
  a1<-bind_rows(a1)
  class(a1) <- c("lmmBoot", class(a1))
  return(a1)
}

boot_ci <- function(object, conf = 0.95){
  if (!inherits(object, "lmmBoot")) stop("object must inherit from class lmmBoot")
  if (nrow(object)==0) stop("boot_par did not return a valid sample")
  objout<-apply(object,2,quantile,probs=c((1-conf)/2,1-(1-conf)/2))
  attr(objout,'nsamples') <- nrow(object)
  return(objout)
}
