#main functions from lmmsmsn package
smsn.lmm <- function(data,formFixed,groupVar,formRandom=~1,depStruct = "CI", timeVar=NULL,
                     distr="sn",pAR=1,luDEC=10,
                     tol=1e-6,max.iter=200,calc.se=T,calc.bi=T,lb=NULL,lu=NULL,
                     initialValues =list(beta=NULL,sigma2=NULL,D=NULL,lambda=NULL,phi=NULL,nu=NULL)) {
                     #beta1,sigmae,D1,lambda,pAR=length(phiAR),phiAR=NULL,
                     #phi1DEC=NULL,phi2DEC=NULL,phiCAR=NULL,phiCS=NULL,
  if (!is.list(initialValues)) stop("initialValues must be a list")
  if (any(!(names(initialValues) %in% c("beta","sigma2","D","lambda","phi","nu")))) stop("initialValues must be a list with
                                                                                         named elements beta, sigma2, D, lambda, phi and/or nu")
  #
  if (!is.character(groupVar)) stop("groupVar must be a character containing the name of the grouping variable in data")
  if (!is.null(timeVar)&!is.character(timeVar)) stop("timeVar must be a character containing the name of the time variable in data")
  if (length(formFixed)!=3) stop("formFixed must be a two-sided linear formula object")
  if (sum(!(c(all.vars(formFixed),all.vars(formRandom),groupVar,timeVar) %in% names(data)))>0) stop("Variable not found in data")
  #
  if (!is.data.frame(data)) stop("data must be a data.frame")
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <-data[,groupVar]
  p<-ncol(x)
  q1<-ncol(z)
  if ((sum(is.na(x))+sum(is.na(z))+sum(is.na(y))+sum(is.na(ind)))>0) stop ("NAs not allowed")
  if (!is.null(timeVar) & sum(is.na(data[,timeVar]))) stop ("NAs not allowed")
  #
  if (!(distr %in% c("sn","st","ss","scn"))) stop("Invalid distribution")
  if ((!is.null(lb))&distr!="sn") if((distr=="st"&(lb<=1))|(distr=="ss"&(lb<=.5))) stop("Invalid lb")
  if (is.null(lb)&distr!="sn") lb = ifelse(distr=="scn",rep(.01,2),ifelse(distr=="st",1.01,.51))
  if (is.null(lu)&distr!="sn") lu = ifelse(distr=="scn",rep(.99,2),ifelse(distr=="st",100,50))
  #
  if (depStruct=="ARp" & !is.null(timeVar) & ((sum(!is.wholenumber(data[,timeVar]))>0)|(sum(data[,timeVar]<=0)>0))) stop("timeVar must contain positive integer numbers when using ARp dependency")
  if (depStruct=="ARp" & !is.null(timeVar)) if (min(data[,timeVar])!=1) warning("consider using a transformation such that timeVar starts at 1")
  if (!(depStruct %in% c("CI","ARp","CS","DEC","CAR1"))) stop("accepted depStruct: CI, ARp, CS, DEC or CAR1")
  #
  if (sum(unlist(lapply(initialValues,is.null)))>0) {
    lmefit = try(lme(formFixed,random=~1|ind,data=data),silent=T)
    if (class(lmefit)=="try-error") stop("error in calculating initial values")
  }
  if (!is.null(initialValues$beta)) {
    beta1 <- initialValues$beta
  } else beta1 <- as.numeric(lmefit$coefficients$fixed)
  if (!is.null(initialValues$sigma2)) {
    sigmae <- initialValues$sigma2
  } else sigmae <- as.numeric(lmefit$sigma^2)
  if (!is.null(initialValues$D)) {
    D1 <- initialValues$D
  } else D1 <- diag(q1)
  if (!is.null(initialValues$lambda)) {
    lambda <- initialValues$lambda
  } else lambda <- rep(3,q1)*sign(as.numeric(skewness(random.effects(lmefit))))
  #
  if (length(D1)==1 & !is.matrix(D1)) D1=as.matrix(D1)
  #
  if (length(beta1)!=p) stop ("wrong dimension of beta")
  if (length(lambda)!=q1) stop ("wrong dimension of lambda")
  if (!is.matrix(D1)) stop("D must be a matrix")
  if ((ncol(D1)!=q1)|(nrow(D1)!=q1)) stop ("wrong dimension of D")
  if (length(sigmae)!=1) stop ("wrong dimension of sigma2")
  if (sigmae<=0) stop("sigma2 must be positive")
  #
  if (depStruct=="ARp") phiAR<- initialValues$phi
  if (depStruct=="CS") phiCS<- initialValues$phi
  if (depStruct=="DEC") parDEC<- initialValues$phi
  if (depStruct=="CAR1") phiCAR1<- initialValues$phi
  #
  nu = initialValues$nu
  #
  if (distr=="st"&is.null(nu)) nu=10
  if (distr=="ss"&is.null(nu)) nu=5
  if (distr=="scn"&is.null(nu)) nu=c(.05,.8)
  #
  if (distr=="st"&length(nu)!=1) stop ("wrong dimension of nu")
  if (distr=="ss"&length(nu)!=1) stop ("wrong dimension of nu")
  if (distr=="scn"&length(nu)!=2) stop ("wrong dimension of nu")
  ###
  if (depStruct=="CI") obj.out <- EM.Skew(formFixed,formRandom,data,groupVar,distr,beta1,sigmae,D1,
                                    lambda,nu,lb,lu,precisao=tol,informa=calc.se,calcbi=calc.bi,max.iter=max.iter)
  if (depStruct=="ARp") obj.out <- EM.SkewAR(formFixed,formRandom,data,groupVar,pAR,timeVar,
                                      distr,beta1,sigmae,phiAR,D1,lambda,nu,lb,lu,
                                      precisao=tol,informa=calc.se,calcbi=calc.bi,max.iter=max.iter)
  if (depStruct=="CS") obj.out <-EM.SkewCS(formFixed,formRandom,data,groupVar,
                                            distr,beta1,sigmae,phiCS,D1,lambda,nu,lb,lu,
                                         precisao=tol,informa=calc.se,calcbi=calc.bi,max.iter=max.iter)
  if (depStruct=="DEC") obj.out <-EM.SkewDEC(formFixed,formRandom,data,groupVar,timeVar,
                                          beta1,sigmae,D1,lambda,distr,nu,parDEC,lb,lu,luDEC,
                                          precisao=tol,informa=calc.se,calcbi=calc.bi,max.iter=max.iter)
  if (depStruct=="CAR1") obj.out <-EM.SkewCAR1(formFixed,formRandom,data,groupVar,timeVar,
                                        distr,beta1,sigmae,phiCAR1,D1,lambda,nu,lb,lu,
                                        precisao=tol,informa=calc.se,calcbi=calc.bi,max.iter=max.iter)
  obj.out$call <- match.call()

  npar<-length(obj.out$theta);N<-nrow(data)
  obj.out$criteria$AIC <- 2*npar-2*obj.out$loglik
  obj.out$criteria$BIC <- log(N)*npar - 2*obj.out$loglik
  obj.out$data = data
  obj.out$formula$formFixed=formFixed
  obj.out$formula$formRandom=formRandom
  obj.out$depStruct = depStruct
  obj.out$distr=distr
  obj.out$N = nrow(data)
  obj.out$n = n_distinct(ind)
  obj.out$groupVar = groupVar
  obj.out$timeVar = timeVar
  #
  if (calc.bi) {
    fitted <- numeric(N)
    ind_levels <- levels(ind)
    for (i in seq_along(ind_levels)) {
      seqi <- ind==ind_levels[i]
      xfiti <- matrix(x[seqi,],ncol=p)
      zfiti <- matrix(z[seqi,],ncol=q1)
      fitted[seqi]<- xfiti%*%obj.out$estimates$beta + zfiti%*%obj.out$random.effects[i,]
    }
    obj.out$fitted <- fitted
  }

  class(obj.out)<- c("SMSN","list")
  obj.out
}

print.SMSN <- function(x,...){
  cat("Linear mixed models with distribution", x$distr, "and dependency structure",x$depStruct,"\n")
  cat("Call:\n")
  print(x$call)
  cat("\nFixed:")
  print(x$formula$formFixed)
  #print(x$theta)
  cat("Random:")
  print(x$formula$formRandom)
  cat("  Estimated variance (D):\n")
  D1 = Dmatrix(x$estimates$dsqrt)%*%Dmatrix(x$estimates$dsqrt)
  colnames(D1)=row.names(D1)= colnames(model.matrix(x$formula$formRandom,data=x$data))
  print(D1)
  cat("\nEstimated parameters:\n")
  if (!is.null(x$std.error)) {
    tab = round(rbind(x$theta,x$std.error),4)
    colnames(tab) = names(x$theta)
    rownames(tab) = c("","s.e.")
  }
  else {
    tab = round(rbind(x$theta),4)
    colnames(tab) = x$theta
    rownames(tab) = c("")
  }
  print(tab)
  cat('\n')
  cat('Model selection criteria:\n')
  critFin <- c(x$loglik, x$criteria$AIC, x$criteria$BIC)
  critFin <- round(t(as.matrix(critFin)),digits=3)
  dimnames(critFin) <- list(c(""),c("logLik", "AIC", "BIC"))
  print(critFin)
  cat('\n')
  cat('Number of observations:',x$N,'\n')
  cat('Number of groups:',x$n,'\n')
}

summary.SMSN <- function(object,confint.level=.95,...){
  cat("Linear mixed models with distribution", object$distr, "and dependency structure",object$depStruct,"\n")
  cat("Call:\n")
  print(object$call)
  cat("\nDistribution", object$distr)
  if (object$distr!="sn") cat(" with nu =", object$estimates$nu,"\n")
  cat("\nRandom effects: ")
  print(object$formula$formRandom)
  cat("  Estimated variance (D):\n")
  D1 = Dmatrix(object$estimates$dsqrt)%*%Dmatrix(object$estimates$dsqrt)
  colnames(D1)=row.names(D1)= colnames(model.matrix(object$formula$formRandom,data=object$data))
  print(D1)
  cat("\nFixed effects: ")
  print(object$formula$formFixed)
  cat("with approximate confidence intervals\n")
  if (!is.null(object$std.error)) {
    p<-length(object$estimates$beta)
    qIC <- qnorm(.5+confint.level/2)
    ICtab <- cbind(object$estimates$beta-qIC*object$std.error[1:p],
                  object$estimates$beta+qIC*object$std.error[1:p])
    tab = (cbind(object$estimates$beta,object$std.error[1:p],
                      ICtab))
    rownames(tab) = names(object$theta[1:p])
    colnames(tab) = c("Value","Std.error",paste0("IC ",confint.level*100,"% lower"),
                      paste0("IC ",confint.level*100,"% upper"))
  }
  else {
    tab = (rbind(object$theta))
    colnames(tab) = names(object$theta[1:p])
    rownames(tab) = c("Value")
  }
  print(tab)
  cat("\nDependency structure:", object$depStruct)
  cat("\n  Estimate(s):\n")
  covParam <- c(object$estimates$sigma2, object$estimates$phi)
  if (object$depStruct=="CI") names(covParam) <- "sigma2"
  else names(covParam) <- c("sigma2",paste0("phi",1:(length(covParam)-1)))
  print(covParam)
  cat("\nSkewness parameter estimate:", object$estimates$lambda)
  cat('\n')
  cat('\nModel selection criteria:\n')
  criteria <- c(object$loglik, object$criteria$AIC, object$criteria$BIC)
  criteria <- round(t(as.matrix(criteria)),digits=3)
  dimnames(criteria) <- list(c(""),c("logLik", "AIC", "BIC"))
  print(criteria)
  cat('\n')
  cat('Number of observations:',object$N,'\n')
  cat('Number of groups:',object$n,'\n')
  invisible(list(varRandom=D1,varFixed=covParam,tableFixed=tab,criteria=criteria))
}

fitted.SMSN <- function(object,...) object$fitted
ranef.SMSN <- function(object,...) object$random.effects

predict.SMSN <- function(object,newData,...){
  dataFit <- object$data
  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<-object$groupVar
  timeVar <- object$timeVar
  dataPred<- newData
  if (sum(!(c(all.vars(formFixed),all.vars(formRandom),groupVar,timeVar) %in% names(newData)))>0) stop("Variable not found in newData")
  depStruct <- object$depStruct
  if (depStruct=="CI") obj.out <- predictf.skew(formFixed,formRandom,dataFit,dataPred,groupVar,distr=object$distr,theta=object$theta)
  if (depStruct=="ARp") obj.out <- predictf.skewAR(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr=object$distr,
                                                  pAR=length(object$estimates$phi),theta=object$theta)
  if (depStruct=="CS") obj.out <-predictf.skewCS(formFixed,formRandom,dataFit,dataPred,groupVar,distr=object$distr,theta=object$theta)
  if (depStruct=="DEC") obj.out <-predictf.skewDEC(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr=object$distr,theta=object$theta)
  if (depStruct=="CAR1") obj.out <-predictf.skewCAR1(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr=object$distr,theta=object$theta)
  obj.out
}

errorVar<- function(times,object=NULL,sigma2=NULL,depStruct=NULL,phi=NULL) {
  if (is.null(object)&is.null(depStruct)) stop("object or depStruct must be provided")
  if (is.null(object)&is.null(sigma2)) stop("object or sigma2 must be provided")
  if (is.null(depStruct)) depStruct<-object$depStruct
  if (depStruct!="CI" & is.null(object)&is.null(phi)) stop("object or phi must be provided")
  if (!(depStruct %in% c("CI","ARp","CS","DEC","CAR1"))) stop("accepted depStruct: CI, ARp, CS, DEC or CAR1")
  if (is.null(sigma2)) sigma2<-object$estimates$sigma2
  if (is.null(phi)&depStruct!="CI") phi<-object$estimates$phi
  if (depStruct=="ARp" & (any(!is.wholenumber(times))|any(times<=0))) stop("times must contain positive integer numbers when using ARp dependency")
  if (depStruct=="ARp" & any(tphitopi(phi)< -1|tphitopi(phi)>1)) stop("AR(p) non stationary, choose other phi")
  #
  if (depStruct=="CI") var.out<- sigma2*diag(length(times))
  if (depStruct=="ARp") var.out<- sigma2*CovARp(phi,times)
  if (depStruct=="CS") var.out<- sigma2*CovCS(phi,length(times))
  if (depStruct=="DEC") var.out<- sigma2*CovDEC(phi[1],phi[2],times)
  if (depStruct=="CAR1") var.out<- sigma2*CovDEC(phi,1,times)
  var.out
}

rsmsn.lmm <- function(time1,x1,z1,sigma2,D1,beta,lambda,depStruct="CI",phi=NULL,distr="sn",nu=NULL) {
  if (length(D1)==1 & !is.matrix(D1)) D1=as.matrix(D1)
  q1 = nrow(D1)
  p = length(beta)
  if (ncol(as.matrix(x1))!=p) stop("incompatible dimension of x1/beta")
  if (ncol(as.matrix(z1))!=q1) stop ("incompatible dimension of z1/D1")
  if (length(lambda)!=q1) stop ("incompatible dimension of lambda/D1")
  if (!is.matrix(D1)) stop("D must be a matrix")
  if ((ncol(D1)!=q1)|(nrow(D1)!=q1)) stop ("wrong dimension of D")
  if (length(sigma2)!=1) stop ("wrong dimension of sigma2")
  if (sigma2<=0) stop("sigma2 must be positive")
  Sig <- errorVar(time1,depStruct = depStruct,sigma2=sigma2,phi=phi)
  #
  if (!(distr %in% c("sn","st","ss","scn"))) stop("Invalid distribution")
  if (distr=="sn") {ui=1; c.=-sqrt(2/pi)}
  if (distr=="st") {ui=rgamma(1,nu/2,nu/2); c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {ui=rbeta(1,nu,1); c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {ui=ifelse(runif(1)<nu[1],nu[2],1);
                      c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%(lambda)))
  Delta = matrix.sqrt(D1)%*%delta
  Gammab = D1 - Delta%*%t(Delta)
  Xi = matrix(x1,ncol=p)
  Zi = matrix(z1,ncol=q1)
  Beta = matrix(beta,ncol=1)
  ti = c.+abs(rnorm(1,0,ui^-.5))
  bi = t(rmvnorm(1,Delta*ti,sigma=ui^(-1)*Gammab))
  Yi = t(rmvnorm(1,Xi%*%Beta+Zi%*%bi,sigma=ui^(-1)*Sig))
  if (all(Xi[,1]==1)) Xi = Xi[,-1]
  if (all(Zi[,1]==1)) Zi = Zi[,-1]
  return(data.frame(time=time1,y=Yi,x=Xi,z=Zi))
}

