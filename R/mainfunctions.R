#main functions from skewlmm package - SMSN-LMM
smsn.lmm <- function(data,formFixed,groupVar,formRandom=~1,depStruct = "UNC", timeVar=NULL,
                     distr="sn",covRandom = 'pdSymm',skewind,pAR=1,
                     control = lmmControl()
                     ) {
  if (!is(formFixed,"formula")) stop("formFixed must be a formula")
  if (!is(formRandom,"formula")) stop("formRandom must be a formula")
  if (!inherits(control,"lmmControl")) stop("control must be a list generated with lmmControl()")
  #
  if (!is.character(groupVar)) stop("groupVar must be a character containing the name of the grouping variable in data")
  if (!is.null(timeVar)&&!is.character(timeVar)) stop("timeVar must be a character containing the name of the time variable in data")
  if (length(formFixed)!=3) stop("formFixed must be a two-sided linear formula object")
  if (!is.data.frame(data)) stop("data must be a data.frame")
  if (length(class(data))>1) data <- as.data.frame(data)
  vars_used<-unique(c(all.vars(formFixed),all.vars(formRandom),groupVar,timeVar))
  vars_miss <- which(!(vars_used %in% names(data)))
  if (length(vars_miss)>0) stop(paste(vars_used[vars_miss],"not found in data"))
  data = data[,vars_used]
  #
  #data <- data[order(data[,groupVar]),]
  if (!is.factor(data[,groupVar])) data[,groupVar]<-haven::as_factor(data[,groupVar])
  data$ind <-data[,groupVar]
  ind <-data[,groupVar]
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  m<-nlevels(ind)#n_distinct(ind)
  if (m<=1) stop(paste(groupVar,"must have more than 1 level"))
  if (all(table(ind)==1)) stop(paste(groupVar,"must have more than 1 observation by level"))
  p<-ncol(x)
  q1<-ncol(z)
  if ((sum(is.na(x))+sum(is.na(z))+sum(is.na(y))+sum(is.na(ind)))>0) stop ("NAs not allowed")
  if (!is.null(timeVar) && sum(is.na(data[,timeVar]))) stop ("NAs not allowed")
  #
  if (distr=="ssl") distr<-"ss"
  distr <- match.arg(distr,c("sn","st","ss","scn"))
  #if (!(distr %in% c("sn","st","ss","scn"))) stop("Accepted distributions: sn, st, ssl, scn")
  if ((!is.null(control$lb))&&distr!="sn") if((distr=="st"&&(control$lb<=1))||(distr=="ss"&&(control$lb<=.5))) stop("Invalid lb")
  if (is.null(control$lb)&&distr!="sn") control$lb = ifelse(distr=="scn",rep(.01,2),ifelse(distr=="st",1.01,.51))
  if (is.null(control$lu)&&distr!="sn") control$lu = ifelse(distr=="scn",rep(.99,2),ifelse(distr=="st",100,50))
  #
  if (depStruct=="ARp" && !is.null(timeVar) &&
      ((sum(!is.wholenumber(data[,timeVar]))>0)||(sum(data[,timeVar]<=0)>0))) stop("timeVar must contain positive integer numbers when using ARp dependence")
  if (depStruct=="ARp" && !is.null(timeVar)) if (min(data[,timeVar])!=1) warning("consider using a transformation such that timeVar starts at 1")
  if (depStruct=="CI") depStruct = "UNC"
  depStruct <- match.arg(depStruct,c("UNC","ARp","CS","DEC","CAR1"))
  #if (!(depStruct %in% c("UNC","ARp","CS","DEC","CAR1"))) stop("accepted depStruct: UNC, ARp, CS, DEC or CAR1")
  #
  covRandom <- match.arg(covRandom, c('pdSymm','pdDiag'))
  #if (!(covRandom %in% c('pdSymm','pdDiag'))) stop("accepted covRandom: pdSymm or pdDiag")
  diagD <- covRandom=='pdDiag'
  if (q1==1) diagD=FALSE
  #
  if (missing(skewind)) skewind<-rep(1,q1)
  if (!all(skewind %in% c(0,1))) stop('skewind must be a vector containing 0s and/or 1s of length equal to the number of random effects')
  if (length(skewind)!=q1) stop('skewind must be a vector containing 0s and/or 1s of length equal to the number of random effects')
  if (all(skewind==0)) stop('for the symmetrical model please use the function smn.lmm')
  #
  if (is.null(control$parallelphi)) control$parallelphi <- ifelse(m>30,TRUE,FALSE)
  if (depStruct=="UNC") control$parallelphi <- FALSE
  if (is.null(control$parallelnu)) {
    if (distr=='st'||distr=='sn') control$parallelnu <- FALSE
    else if (distr=='ss') control$parallelnu <- ifelse(m>30,TRUE,FALSE)
    else control$parallelnu <- ifelse(m>50,TRUE,FALSE)
  }
  if (distr=='sn') control$parallelnu <- FALSE
  if (is.null(control$ncores)) if (control$parallelnu||control$parallelphi) control$ncores <- max(parallel::detectCores() - 1, 1, na.rm = TRUE)
  #
  if (is.null(control$initialValues$beta)||is.null(control$initialValues$sigma2)||
      is.null(control$initialValues$lambda)||is.null(control$initialValues$D)) {
    lmefit = try(lme(formFixed,random=formula(paste('~',as.character(formRandom)[length(formRandom)],
                                                    '|',"ind")),data=data),silent=T)
    if (is(lmefit,"try-error")) {
      lmefit = try(lme(formFixed,random=~1|ind,data=data),silent=TRUE)
      if (is(lmefit,"try-error")) {
        stop("error in calculating initial values")
      } else {
        lambdainit <- rep(1,q1)*sign(as.numeric(skewness(random.effects(lmefit))))
        D1init <- diag(q1)*as.numeric(var(random.effects(lmefit)))
      }
    } else {
      lambdainit <- sign(as.numeric(skewness(random.effects(lmefit))))*1
      D1init <- (var(random.effects(lmefit)))
    }
  }
  if (!is.null(control$initialValues$beta)) {
    beta1 <- control$initialValues$beta
  } else beta1 <- as.numeric(lmefit$coefficients$fixed)
  if (!is.null(control$initialValues$sigma2)) {
    sigmae <- control$initialValues$sigma2
  } else sigmae <- as.numeric(lmefit$sigma^2)
  if (!is.null(control$initialValues$D)) {
    D1 <- control$initialValues$D
  } else D1 <- D1init
  if (!is.null(control$initialValues$lambda)) {
    lambda <- control$initialValues$lambda
  } else lambda <- lambdainit
  lambda <- lambda*skewind
  #
  if (length(D1)==1 && !is.matrix(D1)) D1=as.matrix(D1)
  #
  if (length(beta1)!=p) stop ("wrong dimension of beta")
  if (length(lambda)!=q1) stop ("wrong dimension of lambda")
  if (!is.matrix(D1)) stop("D must be a matrix")
  if ((ncol(D1)!=q1)||(nrow(D1)!=q1)) stop ("wrong dimension of D")
  if (length(sigmae)!=1) stop ("wrong dimension of sigma2")
  if (sigmae<=0) stop("sigma2 must be positive")
  #
  if (depStruct=="ARp") phiAR<- control$initialValues$phi
  if (depStruct=="CS") phiCS<- control$initialValues$phi
  if (depStruct=="DEC") parDEC<- control$initialValues$phi
  if (depStruct=="CAR1") phiCAR1<- control$initialValues$phi
  #
  nu = control$initialValues$nu
  #
  if (distr=="st"&&is.null(nu)) nu=10
  if (distr=="ss"&&is.null(nu)) nu=5
  if (distr=="scn"&&is.null(nu)) nu=c(.05,.8)
  #
  if (distr=="st"&&length(nu)!=1) stop ("wrong dimension of nu")
  if (distr=="ss"&&length(nu)!=1) stop ("wrong dimension of nu")
  if (distr=="scn"&&length(nu)!=2) stop ("wrong dimension of nu")
  ###
  if (depStruct=="UNC") obj.out <- DAAREM.SkewUNC(formFixed, formRandom, data, groupVar,
                                                  distr, beta1, sigmae, D1, lambda, nu,
                                                  lb=control$lb, lu=control$lu, skewind = skewind, diagD=diagD,
                                                  precisao=control$tol, informa=control$calc.se,
                                                  max.iter=control$max.iter, showiter=!control$quiet,
                                                  showerroriter = (!control$quiet)&&control$showCriterium,
                                                  algorithm=control$algorithm,
                                                  control.daarem = control$control.daarem,
                                                  parallelnu=control$parallelnu, ncores=control$ncores)
  if (depStruct=="ARp") obj.out <- DAAREM.SkewAR(formFixed, formRandom, data, groupVar, pAR,
                                                 timeVar, distr, beta1, sigmae,phiAR, D1,
                                                 lambda, nu, lb=control$lb, lu=control$lu, skewind = skewind, diagD=diagD,
                                                 precisao=control$tol, informa=control$calc.se,
                                                 max.iter=control$max.iter, showiter=!control$quiet,
                                                 showerroriter = (!control$quiet)&&control$showCriterium,
                                                 algorithm=control$algorithm,
                                                 control.daarem = control$control.daarem,
                                                 parallelphi = control$parallelphi,
                                                 parallelnu = control$parallelnu, ncores=control$ncores)
  if (depStruct=="CS") obj.out <-DAAREM.SkewCS(formFixed, formRandom, data, groupVar,
                                               distr, beta1, sigmae,phiCS, D1, lambda, nu,
                                               lb=control$lb, lu=control$lu, skewind = skewind, diagD=diagD,
                                               precisao=control$tol, informa=control$calc.se,
                                               max.iter=control$max.iter, showiter=!control$quiet,
                                               showerroriter = (!control$quiet)&&control$showCriterium,
                                               algorithm=control$algorithm,
                                               control.daarem = control$control.daarem,
                                               parallelphi = control$parallelphi,
                                               parallelnu = control$parallelnu, ncores=control$ncores)
  if (depStruct=="DEC") obj.out <-DAAREM.SkewDEC(formFixed, formRandom, data, groupVar,
                                                 timeVar, distr, beta1, sigmae,parDEC, D1,
                                                 lambda, nu, lb=control$lb, lu=control$lu,luDEC=control$luDEC,
                                                 skewind = skewind, diagD=diagD,
                                                 precisao=control$tol, informa=control$calc.se,
                                                 max.iter=control$max.iter, showiter=!control$quiet,
                                                 showerroriter = (!control$quiet)&&control$showCriterium,
                                                 algorithm=control$algorithm,
                                                 control.daarem = control$control.daarem,
                                                 parallelphi = control$parallelphi,
                                                 parallelnu = control$parallelnu, ncores=control$ncores)
  if (depStruct=="CAR1") obj.out <-DAAREM.SkewCAR1(formFixed, formRandom, data, groupVar,
                                                   timeVar, distr, beta1, sigmae,phiCAR1, D1,
                                                   lambda, nu, lb=control$lb, lu=control$lu, skewind = skewind,
                                                   diagD=diagD,precisao=control$tol, informa=control$calc.se,
                                                   max.iter=control$max.iter, showiter=!control$quiet,
                                                   showerroriter = (!control$quiet)&&control$showCriterium,
                                                   algorithm=control$algorithm,
                                                   control.daarem = control$control.daarem,
                                                   parallelphi = control$parallelphi,
                                                   parallelnu = control$parallelnu, ncores=control$ncores)
  obj.out$call <- match.call()

  npar<-length(obj.out$theta);N<-nrow(data)
  obj.out$criteria$AIC <- 2*npar-2*obj.out$loglik
  obj.out$criteria$BIC <- log(N)*npar - 2*obj.out$loglik
  obj.out$data <- data
  obj.out$formula$formFixed <- formFixed
  obj.out$formula$formRandom <- formRandom
  obj.out$formula$groupVar <- groupVar
  obj.out$depStruct <- depStruct
  obj.out$covRandom <- covRandom
  if (distr=="ss") distr<-"ssl"
  obj.out$distr <- distr
  obj.out$N <- N
  obj.out$n <- m#n_distinct(ind)
  obj.out$groupVar <- groupVar
  obj.out$timeVar <- timeVar
  obj.out$control <- control
  obj.out$diagD <- diagD
  obj.out$skewind <- skewind
  #
  fitted <- numeric(N)
  ind_levels <- levels(ind)
  for (i in seq_along(ind_levels)) {
    seqi <- ind==ind_levels[i]
    xfiti <- matrix(x[seqi,],ncol=p)
    zfiti <- matrix(z[seqi,],ncol=q1)
    fitted[seqi]<- xfiti%*%obj.out$estimates$beta + zfiti%*%obj.out$random.effects[i,]
  }
  obj.out$fitted <- fitted
  names(obj.out$estimates$beta) <- colnames(x)
  colnames(obj.out$estimates$D)<- row.names(obj.out$estimates$D) <- colnames(z)

  class(obj.out)<- c("SMSN","list")
  obj.out
}

print.SMSN <- function(x,...){
  cat("Linear mixed models with distribution", x$distr, "and dependence structure",x$depStruct,"\n")
  cat("Log-likelihood value at convergence:", x$loglik)
  cat("\nDistribution", x$distr)
  if (x$distr!="sn") cat(" with nu =", x$estimates$nu)
  cat("\nFixed: ")
  print(x$formula$formFixed)
  print(x$estimates$beta)
  cat("Random effects:\n")
  cat("  Formula: ")
  cat(as.character(x$formula$formRandom), "by", x$groupVar,"\n")
  cat("  Structure:",ifelse(x$covRandom=='pdSymm','General positive-definite',
                             'Diagonal'),'\n')
  cat("  Estimated variance (D):\n")
  D1 = x$estimates$D
  colnames(D1)=row.names(D1)= colnames(model.matrix(x$formula$formRandom,data=x$data))
  print(D1)
  cat("  Skewness parameter:", x$estimates$lambda,"\n")
  cat("Error dependence structure:", x$depStruct)
  cat("\n  Estimate(s):\n")
  covParam <- c(x$estimates$sigma2, x$estimates$phi)
  if (x$depStruct=="UNC") names(covParam) <- "sigma2"
  else names(covParam) <- c("sigma2",paste0("phi",1:(length(covParam)-1)))
  print(covParam)
  cat('Number of observations:',x$N,'\n')
  cat('Number of groups:',x$n,'\n')
}

summary.SMSN <- function(object, confint.level=.95, ...){
  D1 = object$estimates$D
  colnames(D1)=row.names(D1)= colnames(model.matrix(object$formula$formRandom,data=object$data))
  p<-length(object$estimates$beta)
  if (!is.null(object$std.error)) {
    qIC <- qnorm(.5+confint.level/2)
    ICtab <- cbind(object$estimates$beta-qIC*object$std.error[1:p],
                   object$estimates$beta+qIC*object$std.error[1:p])
    tab = (cbind(object$estimates$beta,object$std.error[1:p],
                 ICtab))
    rownames(tab) = names(object$theta[1:p])
    colnames(tab) = c("Value","Std.error",paste0("CI ",confint.level*100,"% lower"),
                      paste0("CI ",confint.level*100,"% upper"))
  }
  else {
    tab = rbind(object$estimates$beta)
    colnames(tab) = names(object$theta[1:p])
    rownames(tab) = c("Value")
  }
  covParam <- c(object$estimates$sigma2, object$estimates$phi)
  if (object$depStruct=="UNC") names(covParam) <- "sigma2"
  else names(covParam) <- c("sigma2",paste0("phi",1:(length(covParam)-1)))
  criteria <- c(object$loglik, object$criteria$AIC, object$criteria$BIC)
  criteria <- round(t(as.matrix(criteria)),digits=3)
  dimnames(criteria) <- list(c(""),c("logLik", "AIC", "BIC"))
  outobj<- list(varRandom=D1,varFixed=covParam,tableFixed=tab,criteria=criteria,
                call = object$call, distr = object$distr, formula = object$formula,
                D = D1, depStruct = object$depStruct,
                estimates = object$estimates, n = object$n, N = object$N,
                covParam = covParam)
  class(outobj) <- c("SMSNsumm","list")
  outobj
}

print.SMSNsumm <- function(x,...){
  cat("Linear mixed models with distribution", x$distr, "and dependence structure",x$depStruct,"\n")
  cat("Call:\n")
  print(x$call)
  cat("\nDistribution", x$distr)
  if (x$distr!="sn") cat(" with nu =", x$estimates$nu,"\n")
  cat("\nRandom effects: \n")
  cat("  Formula: ")
  print(x$formula$formRandom)
  cat("  Structure:",ifelse(x$covRandom=='pdSymm','General positive-definite',
                            'Diagonal'),'\n')
  cat("  Estimated variance (D):\n")
  D1 = x$D
  print(D1)
  cat("\nFixed effects: ")
  print(x$formula$formFixed)
  if (nrow(x$tab)>1) cat("with approximate confidence intervals\n")
  else cat(" (std errors not estimated)\n")
  print(x$tab)
  cat("\nDependence structure:", x$depStruct)
  cat("\n  Estimate(s):\n")
  print(x$covParam)
  cat("\nSkewness parameter estimate:", x$estimates$lambda)
  cat('\n')
  cat('\nModel selection criteria:\n')
  print(x$criteria)
  cat('\n')
  cat('Number of observations:',x$N,'\n')
  cat('Number of groups:',x$n,'\n')
}

fitted.SMSN <- function(object,...) object$fitted
ranef <- function(object) object$random.effects
#
ranef <- function(object, ...) UseMethod("ranef")
nobs <- function(object, ...) UseMethod("nobs")
fixef <- function(object, ...) UseMethod("fixef")
coef <- function(object, ...) UseMethod("coef")
#
ranef.SMN <- ranef.SMSN <- ranef.SMNclmm <- function(object,...) object$random.effects
logLik.SMN <- logLik.SMSN <- logLik.SMNclmm <- function(object,...) object$loglik
fixef.SMN <- fixef.SMSN <- fixef.SMNclmm <- function(object,...) object$estimates$beta
formula.SMN <- formula.SMSN <- formula.SMNclmm <- function(x,...) x$formula
nobs.SMN <- nobs.SMSN <- nobs.SMNclmm <- function(object,...) object$N
sigma.SMN <- sigma.SMSN <- sigma.SMNclmm <- function(object,...) sqrt(object$estimates$sigma2)
#
predict.SMSN <- function(object,newData,...){
  if (missing(newData)||is.null(newData)) return(fitted(object))
    #stop("newData must be a dataset containing the covariates, groupVar and timeVar (when used) from data that should be predicted")
  if (!is.data.frame(newData)) stop("newData must be a data.frame object")
  if (nrow(newData)==0) stop("newData can not be an empty dataset")
  dataFit <- object$data
  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<-object$groupVar
  timeVar <- object$timeVar
  dataPred<- newData
  vars_used<-unique(c(all.vars(formFixed)[-1],all.vars(formRandom),groupVar,timeVar))
  vars_miss <- which(!(vars_used %in% names(newData)))
  if (length(vars_miss)>0) stop(paste(vars_used[vars_miss],"not found in newData"))
  depStruct <- object$depStruct
  if (any(!(dataPred[,groupVar] %in% dataFit[,groupVar]))) stop("subjects for which future values should be predicted must also be at fitting data")
  if (!is.factor(dataFit[,groupVar])) dataFit[,groupVar]<-haven::as_factor(dataFit[,groupVar])
  if (!is.factor(dataPred[,groupVar])) dataPred[,groupVar]<-factor(dataPred[,groupVar],levels=levels(dataFit[,groupVar]))
  #
  if (object$distr=="ssl") object$distr<-"ss"
  dd<-matrix.sqrt(object$estimates$D)[upper.tri(object$estimates$D, diag = T)]
  theta <- c(object$estimates$beta,object$estimates$sigma2,object$estimates$phi,
             dd,object$estimates$lambda,object$estimates$nu)
  #
  if (depStruct=="UNC") obj.out <- predictf.skew(formFixed,formRandom,dataFit,dataPred,groupVar,distr=object$distr,theta=theta)
  if (depStruct=="ARp") obj.out <- predictf.skewAR(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr=object$distr,
                                                  pAR=length(object$estimates$phi),theta=theta)
  if (depStruct=="CS") obj.out <-predictf.skewCS(formFixed,formRandom,dataFit,dataPred,groupVar,distr=object$distr,theta=theta)
  if (depStruct=="DEC") obj.out <-predictf.skewDEC(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr=object$distr,theta=theta)
  if (depStruct=="CAR1") obj.out <-predictf.skewCAR1(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr=object$distr,theta=theta)
  obj.out
}

errorVar<- function(times,object=NULL,sigma2=NULL,depStruct=NULL,phi=NULL) {
  if((!is.null(object))&&(!inherits(object,c("SMSN","SMN")))) stop("object must inherit from class SMSN or SMN")
  if (is.null(object)&&is.null(depStruct)) stop("object or depStruct must be provided")
  if (is.null(object)&&is.null(sigma2)) stop("object or sigma2 must be provided")
  if (is.null(depStruct)) depStruct<-object$depStruct
  if (depStruct=="CI") depStruct = "UNC"
  if (depStruct!="UNC" && is.null(object)&is.null(phi)) stop("object or phi must be provided")
  depStruct <- match.arg(depStruct, c("UNC","ARp","CS","DEC","CAR1"))
  #if (!(depStruct %in% c("UNC","ARp","CS","DEC","CAR1"))) stop("accepted depStruct: UNC, ARp, CS, DEC or CAR1")
  if (is.null(sigma2)) sigma2<-object$estimates$sigma2
  if (is.null(phi)&&depStruct!="UNC") phi<-object$estimates$phi
  if (depStruct=="ARp" && (any(!is.wholenumber(times))|any(times<=0))) stop("times must contain positive integer numbers when using ARp dependence")
  if (depStruct=="ARp" && any(tphitopi(phi)< -1|tphitopi(phi)>1)) stop("AR(p) non stationary, choose other phi")
  #
  if (depStruct=="UNC") var.out<- sigma2*diag(length(times))
  if (depStruct=="ARp") var.out<- sigma2*CovARp(phi,times)
  if (depStruct=="CS") var.out<- sigma2*CovCS(phi,length(times))
  if (depStruct=="DEC") var.out<- sigma2*CovDEC(phi[1],phi[2],times)
  if (depStruct=="CAR1") var.out<- sigma2*CovDEC(phi,1,times)
  var.out
}

rsmsn.lmm <- function(time1,x1,z1,sigma2,D1,beta,lambda,depStruct="UNC",phi=NULL,distr="sn",nu=NULL) {
  if (length(D1)==1 && !is.matrix(D1)) D1=as.matrix(D1)
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
  if (distr=="ssl") distr<-"ss"
  distr <- match.arg(distr, c("sn","st","ss","scn"))
  #if (!(distr %in% c("sn","st","ss","scn"))) stop("Invalid distribution")
  if (distr!="sn"&is.null(nu)) stop("nu must be provided")
  if (distr=="sn") {ui=1; c.=-sqrt(2/pi)}
  if (distr=="st") {ui=rgamma(1,nu/2,nu/2); c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {ui=rbeta(1,nu,1); c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {ui=ifelse(runif(1)<nu[1],nu[2],1);
                      c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  #if (all(lambda==0)) c.=0
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

lr.test <- function(obj1,obj2,level=0.05) {
  if (!all(c(class(obj1)[1],class(obj2)[1])%in%c("SMN","SMSN"))) stop("obj1 and obj2 should be smsn.lmm or smn.lmm objects")
  if (level<=0 || level>=1) stop("0<level<1 needed")
  if (obj1$N!=obj2$N) stop("obj1 and obj2 should refer to the same data set")
  if (obj1$n!=obj2$n) stop("obj1 and obj2 should refer to the same data set")
  npar1 <- length(obj1$theta)
  npar2 <- length(obj2$theta)
  #
  criteria <- rbind(c(obj1$loglik, obj1$criteria$AIC, obj1$criteria$BIC),
                    c(obj2$loglik, obj2$criteria$AIC, obj2$criteria$BIC))
  criteria <- round((as.matrix(criteria)),digits=3)
  dimnames(criteria) <- list(c(deparse(substitute(obj1)),deparse(substitute(obj2))),
                             c("logLik", "AIC", "BIC"))
  if (npar1==npar2) {
    warning("obj1 and obj2 do not contain nested models with different number of parameters")
    return(criteria)
  }
  if (npar1<npar2) {objB <- obj2;objS<-obj1} else {objB <- obj1;objS<-obj2}
  if (objB$depStruct=='DEC') {
    names(objB$theta)[substr(names(objB$theta), 1, 3)=='phi'] = 'phi1'
    names(objB$theta)[substr(names(objB$theta), 1, 5)=='theta'] = 'phi2'
  }
  if (objS$depStruct=='DEC') {
    names(objS$theta)[substr(names(objS$theta), 1, 3)=='phi'] = 'phi1'
    names(objS$theta)[substr(names(objS$theta), 1, 5)=='theta'] = 'phi2'
  }
  if (objB$depStruct =="ARp" || objB$depStruct =="CAR1") names(objB$theta)[substr(names(objB$theta), 1, 3)=='phi'] = paste0('phi',1:length(objB$estimates$phi))
  if (objS$depStruct =="ARp" || objS$depStruct =="CAR1") names(objS$theta)[substr(names(objS$theta), 1, 3)=='phi'] = paste0('phi',1:length(objS$estimates$phi))
  if (!all(names(objS$theta)%in%names(objB$theta))) {
    warning("obj1 and obj2 do not contain nested models")
    return(criteria)
  }
  if ((objB$loglik-objS$loglik)<=0) {
    warning("logLik from model with more parameters is not bigger than the one
    with less parameters. This probably indicates problems on convergence,
    try changing the initial values and/or maximum number of iteration")
    return(criteria)
  }
  lrstat <- 2*(objB$loglik-objS$loglik)
  pval <- pchisq(lrstat,df=abs(npar1-npar2),lower.tail = FALSE)
  out<-list(criteria=criteria,statistic=lrstat,
            p.value=pval,df=abs(npar1-npar2),level=level)
  class(out) <- "lmmLRT"
  out
}

print.lmmLRT <- function(x, ...) {
  cat('\nModel selection criteria:\n')
  print(x$criteria)
  if (!is.null(x$statistic)) {
    cat("\n")
    cat("    Likelihood-ratio Test\n\n")
    cat("chi-square statistics = ",x$statistic,"\n")
    cat("df = ",x$df,"\n")
    cat("p-value = ",x$p.value,"\n")
    if (x$p.value<=x$level) cat("\nThe null hypothesis that both models represent the \ndata equally well is rejected at level ",x$level,'\n')
    else cat("\nThe null hypothesis that both models represent the \ndata equally well is not rejected at level ",x$level,'\n')
  }
}

criteria = function(lobjects) {
  if (!is(lobjects,"list")) stop("lobjects must be a list of SMN, SMSN, or SMNclmm objects")
  if (!all(sapply(lobjects,function(x) class(x)[1] %in% c("SMN","SMSN", "SMNclmm")))) stop("lobjects must be a list of SMN, SMSN, or SMNclmm objects")
  #
  crit.out = t(sapply(lobjects, function(x) c(x$loglik, length(x$theta),
                                             x$criteria$AIC, x$criteria$BIC)))
  colnames(crit.out) = c("logLik", "npar", "AIC", "BIC")
  as.data.frame(crit.out)
}

lmmControl <- function(tol=1e-6,max.iter=300,calc.se=TRUE,
                       lb=NULL,lu=NULL,luDEC=10,
                       initialValues =list(beta=NULL,sigma2=NULL,D=NULL,
                                           lambda=NULL,phi=NULL,nu=NULL),
                       quiet=!interactive(),showCriterium=FALSE,algorithm="DAAREM",
                       parallelphi=NULL, parallelnu=NULL, ncores=NULL,
                       control.daarem=list()) {
  if ((!is.numeric(tol))||(length(tol)>1)) stop("tol must be a small number")
  if (max.iter%%1!=0) {
    max.iter <- ceiling(max.iter)
    warning("using max.iter = ", max.iter)
  }
  if (!is.logical(calc.se)) stop("calc.se must be TRUE or FALSE")
  if ((!is.numeric(luDEC))||luDEC<0) stop("luDEC must be a (not too small) number")
  if (!is.list(initialValues)) stop("initialValues must be a list")
  if (all(c("D","dsqrt") %in% names(initialValues))) initialValues$dsqrt<-NULL
  if (any(!(names(initialValues) %in% c("beta","sigma2","lambda","D","phi","nu")))) warning("initialValues must be a list with named elements beta, sigma2, D, phi and/or nu, elements with other names are ignored")
  if (!is.list(control.daarem)) stop("control.daarem must be a list")
  if (!(algorithm%in%c("DAAREM","EM"))) algorithm = match.arg(algorithm, c("DAAREM","EM"))#stop("algorithm must be either 'EM' or 'DAAREM'")
  if (!is.logical(quiet)) stop("quiet must be TRUE or FALSE")
  if (!is.logical(showCriterium)) stop("showCriterium must be TRUE or FALSE")
  if (!is.null(parallelphi)) if (!is.logical(parallelphi)) stop("parallelphi must be TRUE or FALSE")
  if (!is.null(parallelnu)) if (!is.logical(parallelnu)) stop("parallelnu must be TRUE or FALSE")
  if (!is.null(ncores)) if (ncores%%1!=0) {
    ncores <- ceiling(ncores)
    warning("using ncores = ", ncores)
  }
  out<-list(tol = tol, max.iter = max.iter, calc.se = calc.se,
            lb = lb, lu = lu, luDEC = luDEC, initialValues = initialValues,
            quiet = quiet, showCriterium = showCriterium, algorithm = algorithm,
            parallelphi = parallelphi, parallelnu = parallelnu, ncores = ncores,
            control.daarem = control.daarem)
  class(out) <- c("lmmControl","list")
  return(out)
}

# based on nlme update.lme
update.SMSN <- update.SMN <- function (object, ..., evaluate = TRUE){
  call <- object$call
  if (is.null(call))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if(length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    ## do these individually to allow NULL to remove entries.
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if(evaluate) eval(call, parent.frame())
  else call
}

# adding coef function
coef.SMSN <- coef.SMN <- coef.SMNclmm <- function(object, ...){
  fef <- data.frame(matrix(object$estimates$beta, ncol=length(object$estimates$beta),
                           nrow=object$n, byrow=T), check.names = FALSE)
  names(fef) <- names(object$estimates$beta)
  ref <- data.frame(skewlmm::ranef(object), check.names = FALSE)
  ##
  varsnames <- union(colnames(fef), colnames(ref))
  fef[,setdiff(varsnames, colnames(fef))] <- 0
  ref[,setdiff(varsnames, colnames(ref))] <- 0
  ##
  common_cols <- intersect(colnames(fef), colnames(ref))
  coef <- data.frame(sapply(common_cols, function(col) fef[[col]] + ref[[col]]), check.names = FALSE)
  class(coef) <- c("coef", "data.frame")
  coef
}

# adding confint method
confint.SMSN <- confint.SMN <- function(object, parm, level = 0.95, method = "asymptotic", ...){
  if (is.null(object$std.error)) stop("A numerical error prevented calculation of standard errors. Please consider changing the model, the algorithm, or the initial values")
  if (missing(parm)) {
    parm = "all"
  } else parm <- match.arg(parm, c("beta","all"))
  method <- match.arg(method, c("asymptotic","bootstrap"))
  if (level>=1|level<=0) stop("level must be a number between 0 and 1")
  p <- length(object$estimates$beta)
  if (method == "asymptotic") {
    qIC <- qnorm(.5+level/2)
    if (parm == "beta") {
      ICtab <- cbind(object$estimates$beta-qIC*object$std.error[1:p],
                     object$estimates$beta+qIC*object$std.error[1:p])
      tab = (cbind(object$estimates$beta, ICtab))
      rownames(tab) = names(object$theta[1:p])
      colnames(tab) = c("Estimate",paste0("CI ",level*100,"% lower"),
                        paste0("CI ",level*100,"% upper"))
    }
    else {
      tab <- cbind(object$theta,
                     object$theta-qIC*object$std.error,
                     object$theta+qIC*object$std.error)
      rownames(tab) = names(object$theta)
      colnames(tab) = c("Estimate",paste0("CI ",level*100,"% lower"),
                        paste0("CI ",level*100,"% upper"))
    }
  } else {
    message("Computing bootstrap intervals...")
    boot_sample <- boot_par(object, ...)
    tab <- cbind(object$theta, t(boot_ci(boot_sample)))
    colnames(tab)[1] <- "Estimate"
    if (parm == "beta") {
      tab <- tab[1:p,]
    }
  }
  return(tab)
}
