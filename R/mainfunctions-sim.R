#main functions from skewlmm package - SMN-LMM
smn.lmm <- function(data,formFixed,groupVar,formRandom=~1,depStruct = "CI", timeVar=NULL,
                    distr="norm",pAR=1,luDEC=10,
                    tol=1e-6,max.iter=200,calc.se=TRUE,lb=NULL,lu=NULL,
                    initialValues =list(beta=NULL,sigma2=NULL,D=NULL,phi=NULL,nu=NULL),
                    quiet=FALSE,showCriterium=FALSE) {
  if (class(formFixed)!="formula") stop("formFixed must be a formula")
  if (class(formRandom)!="formula") stop("formRandom must be a formula")
  if (!is.list(initialValues)) stop("initialValues must be a list")
  if (all(c("D","dsqrt") %in% names(initialValues))) initialValues$dsqrt<-NULL
  if (any(!(names(initialValues) %in% c("beta","sigma2","D","lambda","phi","nu")))) warning("initialValues must be a list with named elements beta, sigma2, D, lambda, phi and/or nu, elements with other names are ignored")
  #
  if (!is.character(groupVar)) stop("groupVar must be a character containing the name of the grouping variable in data")
  if (!is.null(timeVar)&!is.character(timeVar)) stop("timeVar must be a character containing the name of the time variable in data")
  if (length(formFixed)!=3) stop("formFixed must be a two-sided linear formula object")
  if (!is.data.frame(data)) stop("data must be a data.frame")
  #if (is_tibble(data)) data=as.data.frame(data)
  if (length(class(data))>1) data=as.data.frame(data)
  vars_used<-unique(c(all.vars(formFixed),all.vars(formRandom),groupVar,timeVar))
  vars_miss <- which(!(vars_used %in% names(data)))
  if (length(vars_miss)>0) stop(paste(vars_used[vars_miss],"not found in data"))
  data = data[,vars_used]
  #
  if (!is.factor(data[,groupVar])) data[,groupVar]<-as.factor(data[,groupVar])
  x <- model.matrix(formFixed,data=data)
  y <-data[,all.vars(formFixed)[1]]
  z<-model.matrix(formRandom,data=data)
  ind <-data[,groupVar]
  data$ind <-data[,groupVar]
  m<-nlevels(ind)#n_distinct(ind)
  if (m<=1) stop(paste(groupVar,"must have more than 1 level"))
  p<-ncol(x)
  q1<-ncol(z)
  if ((sum(is.na(x))+sum(is.na(z))+sum(is.na(y))+sum(is.na(ind)))>0) stop ("NAs not allowed")
  if (!is.null(timeVar) & sum(is.na(data[,timeVar]))) stop ("NAs not allowed")
  #
  if (!(distr %in% c("norm","t","sl","cn"))) stop("Accepted distributions: norm, t, sl, cn")
  if ((!is.null(lb))&distr!="norm") if((distr=="t"&(lb<=1))|(distr=="sl"&(lb<=.5))) stop("Invalid lb")
  if (is.null(lb)&distr!="norm") lb = ifelse(distr=="cn",rep(.01,2),ifelse(distr=="t",1.01,.51))
  if (is.null(lu)&distr!="norm") lu = ifelse(distr=="cn",rep(.99,2),ifelse(distr=="t",100,50))
  #
  if (depStruct=="ARp" & !is.null(timeVar) & ((sum(!is.wholenumber(data[,timeVar]))>0)|(sum(data[,timeVar]<=0)>0))) stop("timeVar must contain positive integer numbers when using ARp dependency")
  if (depStruct=="ARp" & !is.null(timeVar)) if (min(data[,timeVar])!=1) warning("consider using a transformation such that timeVar starts at 1")
  if (!(depStruct %in% c("CI","ARp","CS","DEC","CAR1"))) stop("accepted depStruct: CI, ARp, CS, DEC or CAR1")
  #
  if (is.null(initialValues$beta)|is.null(initialValues$sigma2)|is.null(initialValues$D)) {
    lmefit = try(lme(formFixed,random=formula(paste('~',as.character(formRandom)[length(formRandom)],
                                                    '|',"ind")),data=data),silent=T)
    if (class(lmefit)=="try-error") {
      lmefit = try(lme(formFixed,random=~1|ind,data=data),silent=TRUE)
      if (class(lmefit)=="try-error") {
        stop("error in calculating initial values")
      } else {
        D1init <- diag(q1)*as.numeric(var(random.effects(lmefit)))
      }
    } else {
      D1init <- (var(random.effects(lmefit)))
    }
  }
  if (!is.null(initialValues$beta)) {
    beta1 <- initialValues$beta
  } else beta1 <- as.numeric(lmefit$coefficients$fixed)
  if (!is.null(initialValues$sigma2)) {
    sigmae <- initialValues$sigma2
  } else sigmae <- as.numeric(lmefit$sigma^2)
  if (!is.null(initialValues$D)) {
    D1 <- initialValues$D
  } else D1 <- D1init
  #
  if (length(D1)==1 & !is.matrix(D1)) D1=as.matrix(D1)
  #
  if (length(beta1)!=p) stop ("wrong dimension of beta")
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
  if (distr=="t"&is.null(nu)) nu=10
  if (distr=="sl"&is.null(nu)) nu=5
  if (distr=="cn"&is.null(nu)) nu=c(.05,.8)
  #
  if (distr=="t"&length(nu)!=1) stop ("wrong dimension of nu")
  if (distr=="sl"&length(nu)!=1) stop ("wrong dimension of nu")
  if (distr=="cn"&length(nu)!=2) stop ("wrong dimension of nu")
  #
  if (distr=="norm") distrs="sn"
  if (distr=="t") distrs="st"
  if (distr=="sl") distrs="ss"
  if (distr=="cn") distrs="scn"
  ###
  if (depStruct=="CI") obj.out <- EM.sim(formFixed,formRandom,data,groupVar,distr=distrs,beta1,sigmae,D1,
                                         nu,lb,lu,precisao=tol,informa=calc.se,max.iter=max.iter,showiter=!quiet,showerroriter = (!quiet)&showCriterium)
  if (depStruct=="ARp") obj.out <- EM.AR(formFixed,formRandom,data,groupVar,pAR,timeVar,
                                         distr=distrs,beta1,sigmae,phiAR,D1,nu,lb,lu,
                                         precisao=tol,informa=calc.se,max.iter=max.iter,showiter=!quiet,showerroriter = (!quiet)&showCriterium)
  if (depStruct=="CS") obj.out <-EM.CS(formFixed,formRandom,data,groupVar,
                                       distr=distrs,beta1,sigmae,phiCS,D1,nu,lb,lu,
                                       precisao=tol,informa=calc.se,max.iter=max.iter,showiter=!quiet,showerroriter = (!quiet)&showCriterium)
  if (depStruct=="DEC") obj.out <-EM.DEC(formFixed,formRandom,data,groupVar,timeVar,
                                         beta1,sigmae,D1,distr=distrs,nu,parDEC,lb,lu,luDEC,
                                         precisao=tol,informa=calc.se,max.iter=max.iter,showiter=!quiet,showerroriter = (!quiet)&showCriterium)
  if (depStruct=="CAR1") obj.out <-EM.CAR1(formFixed = formFixed,formRandom = formRandom,
                                           data = data,groupVar = groupVar,timeVar = timeVar,
                                           distr=distrs,beta1 = beta1,sigmae = sigmae,
                                           phiCAR1 = phiCAR1,D1 = D1,nu = nu,lb = lb,lu = lu,
                                           precisao=tol,informa=calc.se,max.iter=max.iter,showiter=!quiet,showerroriter = (!quiet)&showCriterium)
  obj.out$call <- match.call()

  npar<-length(obj.out$theta);N<-nrow(data)
  obj.out$criteria$AIC <- 2*npar-2*obj.out$loglik
  obj.out$criteria$BIC <- log(N)*npar - 2*obj.out$loglik
  obj.out$data = data
  obj.out$formula$formFixed=formFixed
  obj.out$formula$formRandom=formRandom
  obj.out$depStruct = depStruct
  obj.out$distr=distr
  obj.out$N = N
  obj.out$n = m#n_distinct(ind)
  obj.out$groupVar = groupVar
  obj.out$timeVar = timeVar
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

  class(obj.out)<- c("SMN","list")
  obj.out
}

print.SMN <- function(x,...){
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
    colnames(tab) = names(x$theta)
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

summary.SMN <- function(object,confint.level=.95,...){
  cat("Linear mixed models with distribution", object$distr, "and dependency structure",object$depStruct,"\n")
  cat("Call:\n")
  print(object$call)
  cat("\nDistribution", object$distr)
  if (object$distr!="norm") cat(" with nu =", object$estimates$nu,"\n")
  cat("\nRandom effects: ")
  print(object$formula$formRandom)
  cat("  Estimated variance (D):\n")
  D1 = Dmatrix(object$estimates$dsqrt)%*%Dmatrix(object$estimates$dsqrt)
  colnames(D1)=row.names(D1)= colnames(model.matrix(object$formula$formRandom,data=object$data))
  print(D1)
  cat("\nFixed effects: ")
  print(object$formula$formFixed)
  if (!is.null(object$std.error)) cat("with approximate confidence intervals\n")
  else cat(" (std errors not estimated)\n")
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
  } else {
    tab = (rbind(object$estimates$beta))
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

fitted.SMN <- function(object,...) object$fitted
#ranef.SMN <- function(object,...) object$random.effects

#colocar if subj not in data
predict.SMN <- function(object,newData,...){
  if (missing(newData)) stop("newData must be a dataset containing the covariates, groupVar and timeVar (when used) from data that should be predicted")
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
  if (any(!(dataPred[,groupVar] %in% dataFit[,groupVar]))) stop("subjects for which future values should be predicted must also be at fitting data")
  if (!is.factor(dataFit[,groupVar])) dataFit[,groupVar]<-as.factor(dataFit[,groupVar])
  if (!is.factor(dataPred[,groupVar])) dataPred[,groupVar]<-factor(dataPred[,groupVar],levels=levels(dataFit[,groupVar]))
  #
  if (depStruct=="CI") obj.out <- predictf.sim(formFixed,formRandom,dataFit,dataPred,groupVar,distr=distrs,theta=object$theta)
  if (depStruct=="ARp") obj.out <- predictf.AR(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr=distrs,
                                                   pAR=length(object$estimates$phi),theta=object$theta)
  if (depStruct=="CS") obj.out <-predictf.CS(formFixed,formRandom,dataFit,dataPred,groupVar,distr=distrs,theta=object$theta)
  if (depStruct=="DEC") obj.out <-predictf.DEC(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr=distrs,theta=object$theta)
  if (depStruct=="CAR1") obj.out <-predictf.CAR1(formFixed,formRandom,dataFit,dataPred,groupVar,timeVar,distr=distrs,theta=object$theta)
  obj.out
}
