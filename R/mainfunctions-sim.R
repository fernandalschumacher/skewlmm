#main functions from skewlmm package - SMN-LMM
smn.lmm <- function(data,formFixed,groupVar,formRandom=~1,depStruct = "UNC",
                    timeVar=NULL,distr="norm",covRandom='pdSymm',pAR=1,
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
  #if (is_tibble(data)) data=as.data.frame(data)
  if (length(class(data))>1) data=as.data.frame(data)
  vars_used<-unique(c(all.vars(formFixed),all.vars(formRandom),groupVar,timeVar))
  vars_miss <- which(!(vars_used %in% names(data)))
  if (length(vars_miss)>0) stop(paste(vars_used[vars_miss],"not found in data"))
  data = data[,vars_used]
  #
  #data <- data[order(data[,groupVar]),]
  if (!is.factor(data[,groupVar])) data[,groupVar]<-haven::as_factor(data[,groupVar])
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <-data[,groupVar]
  m<-nlevels(ind)#n_distinct(ind)
  if (m<=1) stop(paste(groupVar,"must have more than 1 level"))
  if (all(table(ind)==1)) stop(paste(groupVar,"must have more than 1 observation by level"))
  p<-ncol(x)
  q1<-ncol(z)
  if ((sum(is.na(x))+sum(is.na(z))+sum(is.na(y))+sum(is.na(ind)))>0) stop ("NAs not allowed")
  if (!is.null(timeVar) && sum(is.na(data[,timeVar]))) stop ("NAs not allowed")
  #
  distr <- match.arg(distr, c("norm","t","sl","cn"))
  #if (!(distr %in% c("norm","t","sl","cn"))) stop("Accepted distributions: norm, t, sl, cn")
  if ((!is.null(control$lb))&&distr!="norm") if((distr=="t"&&(control$lb<=1))||(distr=="sl"&&(control$lb<=.5))) stop("Invalid lb")
  if (is.null(control$lb)&&distr!="norm") control$lb = ifelse(distr=="cn",rep(.01,2),ifelse(distr=="t",1.01,.51))
  if (is.null(control$lu)&&distr!="norm") control$lu = ifelse(distr=="cn",rep(.99,2),ifelse(distr=="t",100,50))
  #
  if (depStruct=="ARp" && !is.null(timeVar) &&
      ((sum(!is.wholenumber(data[,timeVar]))>0)||(sum(data[,timeVar]<=0)>0))) stop("timeVar must contain positive integer numbers when using ARp dependence")
  if (depStruct=="ARp" && !is.null(timeVar)) if (min(data[,timeVar])!=1) warning("consider using a transformation such that timeVar starts at 1")
  if (depStruct=="CI") depStruct = "UNC"
  depStruct <- match.arg(depStruct, c("UNC","ARp","CS","DEC","CAR1"))
  #if (!(depStruct %in% c("UNC","ARp","CS","DEC","CAR1"))) stop("accepted depStruct: UNC, ARp, CS, DEC or CAR1")
  #
  covRandom <- match.arg(covRandom, c('pdSymm','pdDiag'))
  #if (!(covRandom %in% c('pdSymm','pdDiag'))) stop("accepted covRandom: pdSymm or pdDiag")
  diagD <- covRandom=='pdDiag'
  if (q1==1) diagD=FALSE
  #
  if (is.null(control$parallelphi)) control$parallelphi <- ifelse(m>30,TRUE,FALSE)
  if (depStruct=="UNC") control$parallelphi <- FALSE
  if (is.null(control$parallelnu)) {
    if (distr=='t'||distr=='norm') control$parallelnu <- FALSE
    else if (distr=='sl') control$parallelnu <- ifelse(m>30,TRUE,FALSE)
    else control$parallelnu <- ifelse(m>50,TRUE,FALSE)
  }
  if (distr=='norm') control$parallelnu <- FALSE
  if (is.null(control$ncores)) if (control$parallelnu||control$parallelphi) control$ncores <- max(parallel::detectCores() - 1, 1, na.rm = TRUE)
  #
  if (is.null(control$initialValues$beta)||is.null(control$initialValues$sigma2)||
      is.null(control$initialValues$D)) {
    lmefit = try(lme(formFixed,random=formula(paste('~',as.character(formRandom)[length(formRandom)],
                                                    '|',"ind")),data=data),silent=T)
    if (is(lmefit,"try-error")) {
      lmefit = try(lme(formFixed,random=~1|ind,data=data),silent=TRUE)
      if (is(lmefit,"try-error")) {
        stop("error in calculating initial values")
      } else {
        D1init <- diag(q1)*as.numeric(var(random.effects(lmefit)))
      }
    } else {
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
  #
  if (length(D1)==1 && !is.matrix(D1)) D1=as.matrix(D1)
  #
  if (length(beta1)!=p) stop ("wrong dimension of beta")
  if (!is.matrix(D1)) stop("D must be a matrix")
  if ((ncol(D1)!=q1)|(nrow(D1)!=q1)) stop ("wrong dimension of D")
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
  if (distr=="t"&&is.null(nu)) nu=10
  if (distr=="sl"&&is.null(nu)) nu=5
  if (distr=="cn"&&is.null(nu)) nu=c(.05,.8)
  #
  if (distr=="t"&&length(nu)!=1) stop ("wrong dimension of nu")
  if (distr=="sl"&&length(nu)!=1) stop ("wrong dimension of nu")
  if (distr=="cn"&&length(nu)!=2) stop ("wrong dimension of nu")
  #
  if (distr=="norm") distrs="sn"
  if (distr=="t") distrs="st"
  if (distr=="sl") distrs="ss"
  if (distr=="cn") distrs="scn"
  ###
  if (depStruct=="UNC") obj.out <- DAAREM.UNC(formFixed, formRandom, data, groupVar,
                                              distr=distrs, beta1, sigmae, D1, nu, lb=control$lb, lu=control$lu,
                                              diagD = diagD, precisao=control$tol, informa=control$calc.se,
                                              max.iter=control$max.iter,showiter=!control$quiet,
                                              showerroriter = (!control$quiet)&&control$showCriterium,
                                              algorithm=control$algorithm, control.daarem = control$control.daarem,
                                              parallelnu = control$parallelnu, ncores = control$ncores)
  if (depStruct=="ARp") obj.out <- DAAREM.AR(formFixed, formRandom, data, groupVar,pAR,timeVar,
                                             distr=distrs, beta1, sigmae,phiAR, D1, nu, lb=control$lb, lu=control$lu,
                                             diagD = diagD,precisao=control$tol, informa=control$calc.se,
                                             max.iter=control$max.iter,showiter=!control$quiet,
                                             showerroriter = (!control$quiet)&&control$showCriterium,
                                             algorithm=control$algorithm, control.daarem = control$control.daarem,
                                             parallelphi = control$parallelphi, parallelnu = control$parallelnu,
                                             ncores = control$ncores)
  if (depStruct=="CS") obj.out <- DAAREM.CS(formFixed, formRandom, data, groupVar,
                                           distr=distrs, beta1, sigmae,phiCS, D1, nu,
                                           lb=control$lb, lu=control$lu, diagD = diagD,
                                           precisao=control$tol, informa=control$calc.se,
                                           max.iter=control$max.iter,showiter=!control$quiet,
                                           showerroriter = (!control$quiet)&&control$showCriterium,
                                           algorithm=control$algorithm, control.daarem = control$control.daarem,
                                           parallelphi = control$parallelphi, parallelnu = control$parallelnu,
                                           ncores = control$ncores)
  if (depStruct=="DEC") obj.out <- DAAREM.DEC(formFixed, formRandom, data, groupVar,timeVar,
                                             distr=distrs, beta1, sigmae,parDEC, D1, nu,
                                             lb=control$lb, lu=control$lu, luDEC=control$luDEC,
                                             diagD = diagD,precisao=control$tol, informa=control$calc.se,
                                             max.iter=control$max.iter,showiter=!control$quiet,
                                             showerroriter = (!control$quiet)&&control$showCriterium,
                                             algorithm=control$algorithm, control.daarem = control$control.daarem,
                                             parallelphi = control$parallelphi, parallelnu = control$parallelnu,
                                             ncores = control$ncores)
  if (depStruct=="CAR1") obj.out <-DAAREM.CAR1(formFixed, formRandom, data, groupVar,timeVar,
                                              distr=distrs, beta1, sigmae,phiCAR1, D1, nu,
                                              lb=control$lb, lu=control$lu, diagD = diagD,
                                              precisao=control$tol, informa=control$calc.se,
                                              max.iter=control$max.iter,showiter=!control$quiet,
                                              showerroriter = (!control$quiet)&&control$showCriterium,
                                              algorithm=control$algorithm, control.daarem = control$control.daarem,
                                              parallelphi = control$parallelphi, parallelnu = control$parallelnu,
                                              ncores = control$ncores)
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
  obj.out$distr <- distr
  obj.out$N <- N
  obj.out$n <- m#n_distinct(ind)
  obj.out$groupVar <- groupVar
  obj.out$timeVar <- timeVar
  obj.out$control <- control
  obj.out$diagD <- diagD
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

  class(obj.out)<- c("SMN","list")
  obj.out
}

print.SMN <- function(x,...){
  cat("Linear mixed models with distribution", x$distr, "and dependence structure",x$depStruct,"\n")
  cat("Log-likelihood value at convergence:", x$loglik)
  cat("\nDistribution", x$distr)
  if (x$distr!="norm") cat(" with nu =", x$estimates$nu,"\n")
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
  cat("Error dependence structure:", x$depStruct)
  cat("\n  Estimate(s):\n")
  covParam <- c(x$estimates$sigma2, x$estimates$phi)
  if (x$depStruct=="UNC") names(covParam) <- "sigma2"
  else names(covParam) <- c("sigma2",paste0("phi",1:(length(covParam)-1)))
  print(covParam)
  cat('Number of observations:',x$N,'\n')
  cat('Number of groups:',x$n,'\n')
}

summary.SMN <- function(object, confint.level=.95, ...){
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
  class(outobj) <- c("SMNsumm","list")
  outobj
}

print.SMNsumm <- function(x, ...){
  cat("Linear mixed models with distribution", x$distr, "and dependence structure",x$depStruct,"\n")
  cat("Call:\n")
  print(x$call)
  cat("\nDistribution", x$distr)
  if (x$distr!="norm") cat(" with nu =", x$estimates$nu,"\n")
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
  cat('\nModel selection criteria:\n')
  print(x$criteria)
  cat('\n')
  cat('Number of observations:',x$N,'\n')
  cat('Number of groups:',x$n,'\n')
}

fitted.SMN <- function(object,...) object$fitted
#
predict.SMN <- function(object,newData,...){
  if (missing(newData)||is.null(newData)) return(fitted(object))
  if (!is.data.frame(newData)) stop("newData must be a data.frame object")
  if (nrow(newData)==0) stop("newData can not be an empty dataset")
  dataFit <- object$data
  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<-object$groupVar
  timeVar <- object$timeVar
  dataPred<- newData
  #
  if (object$distr=="norm") distrs="sn"
  if (object$distr=="t") distrs="st"
  if (object$distr=="sl") distrs="ss"
  if (object$distr=="cn") distrs="scn"
  #
  vars_used<-unique(c(all.vars(formFixed)[-1],all.vars(formRandom),groupVar,timeVar))
  vars_miss <- which(!(vars_used %in% names(newData)))
  if (length(vars_miss)>0) stop(paste(vars_used[vars_miss],"not found in newData"))
  depStruct <- object$depStruct
  if (depStruct=="CI") depStruct = "UNC"
  if (any(!(dataPred[,groupVar] %in% dataFit[,groupVar]))) stop("subjects for which future values should be predicted must also be at fitting data")
  if (!is.factor(dataFit[,groupVar])) dataFit[,groupVar]<-haven::as_factor(dataFit[,groupVar])
  if (!is.factor(dataPred[,groupVar])) dataPred[,groupVar]<-factor(dataPred[,groupVar],levels=levels(dataFit[,groupVar]))
  dd<-matrix.sqrt(object$estimates$D)[upper.tri(object$estimates$D, diag = T)]
  theta <- c(object$estimates$beta,object$estimates$sigma2,object$estimates$phi,
             dd,object$estimates$nu)
  #
  if (depStruct=="UNC") obj.out <- predictf.sim(formFixed,formRandom,dataFit,dataPred,groupVar,distr=distrs,theta=theta)
  if (depStruct=="ARp") obj.out <- predictf.AR(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr=distrs,
                                                   pAR=length(object$estimates$phi),theta=theta)
  if (depStruct=="CS") obj.out <-predictf.CS(formFixed,formRandom,dataFit,dataPred,groupVar,distr=distrs,theta=theta)
  if (depStruct=="DEC") obj.out <-predictf.DEC(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr=distrs,theta=theta)
  if (depStruct=="CAR1") obj.out <-predictf.CAR1(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr=distrs,theta=theta)
  obj.out
}
