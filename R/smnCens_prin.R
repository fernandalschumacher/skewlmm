###########################################################
##   Linear Mixed-Effects Models with Censored Response  ##
###########################################################

# Censored mixed-effects models for irregularly observed repeated measures
# ------------------------------------------------------------------------------
smnCens.lmm = function(data, formFixed, groupVar, formRandom=~1, depStruct="UNC",
                       ci=NULL, lcl=NULL, ucl=NULL, timeVar=NULL, distr="norm", nufix=FALSE,
                       pAR=NULL, control=lmmControl()){

  if (!is(formFixed,"formula")) stop("formFixed must be a formula")
  if (!is(formRandom,"formula")) stop("formRandom must be a formula")
  if (!inherits(control,"lmmControl")) stop("control must be a list generated with lmmControl()")
  if (!is.logical(nufix)) stop ("nufix must be TRUE or FALSE")
  #
  if (!is.character(groupVar)) stop("groupVar must be a character containing the name of the grouping variable in data")
  if (!is.null(timeVar)&&!is.character(timeVar)) stop("timeVar must be a character containing the name of the time variable in data")
  if (!is.null(ci)&&!is.character(ci)) stop("ci must be a character containing the name of the censoring indicator in data")
  if (!is.null(lcl)&&!is.character(lcl)) stop("lcl must be a character containing the name of the lower censoring limit in data")
  if (!is.null(ucl)&&!is.character(ucl)) stop("ucl must be a character containing the name of the upper censoring limit in data")
  if (length(formFixed)!=3) stop("formFixed must be a two-sided linear formula object")
  if (!is.data.frame(data)) stop("data must be a data.frame")
  if (length(class(data))>1) data=as.data.frame(data)
  vars_used = unique(c(all.vars(formFixed), all.vars(formRandom), groupVar, timeVar, ci, lcl, ucl))
  vars_miss = which(!(vars_used %in% names(data)))
  if (length(vars_miss)>0) stop(paste(vars_used[vars_miss],"not found in data"))
  data = data[,vars_used]
  #
  if (!is.factor(data[,groupVar])) data[,groupVar] = haven::as_factor(data[,groupVar])
  x  = model.matrix(formFixed, data=model.frame(formFixed,data,na.action=NULL))
  y = data[,all.vars(formFixed)[1]]
  z = model.matrix(formRandom, data=data)
  ind = data[,groupVar]
  data$ind = ind
  #
  N = length(c(y))
  if (is.null(ci)) { cc = rep(0, N)} else { cc = data[,ci] }
  if (is.null(lcl)){ lcl = rep(-Inf, N)} else { lcl = data[,lcl] }
  if (is.null(ucl)){ ucl = rep(Inf, N) } else { ucl = data[,ucl] }
  yna = which(is.na(y))
  if (sum(cc%in%c(0,1)) < length(cc)) stop("The elements of ci must be 0 or 1")
  if (sum(cc[yna])!=length(yna)) stop("NA values in the response must be specified through arguments ci, lcl, and ucl")
  if (sum(cc) > 0){
    if (length(yna)>0){
      censor = (cc==1 & !is.na(y))
      if (any(is.infinite(lcl[censor]) & is.infinite(ucl[censor]))) stop("lcl or ucl must be finite for censored data")
    } else {
      if (any(is.infinite(lcl[cc==1]) & is.infinite(ucl[cc==1]))) stop("lcl or ucl must be finite for censored data")
    }
    if (sum(is.na(lcl))>0 | sum(is.na(ucl))>0) stop("There are some NA values in lcl or ucl")
    if (!all(lcl[cc==1]<ucl[cc==1])) stop ("lcl must be smaller than ucl")
  }
  #
  m = nlevels(ind)
  if (m<=1) stop(paste(groupVar,"must have more than 1 level"))
  if (all(table(ind)==1)) stop(paste(groupVar,"must have more than 1 observation by level"))
  p = ncol(x)
  q1 = ncol(z)
  if ((sum(is.na(x))+sum(is.na(z))+sum(is.na(y))+sum(is.na(ind)))>0) stop ("NAs not allowed")
  if (!is.null(timeVar) && sum(is.na(data[,timeVar]))) stop ("NAs not allowed")
  nj = tapply(ind, ind, length)[unique(ind)]
  if (is.null(timeVar)) {
    ttc = flatten_int(lapply(nj, function(x) seq(1, x)))
  } else { ttc = data[,timeVar] }
  #
  if (!(distr %in% c("norm","t"))) stop("Accepted distributions: norm and t")
  if (depStruct=="ARp" && !is.null(timeVar) &&
      ((sum(!is.wholenumber(data[,timeVar]))>0)||(sum(data[,timeVar]<=0)>0))) stop("timeVar must contain positive integer numbers when using ARp dependency")
  if (depStruct=="ARp" && !is.null(timeVar)) if (min(data[,timeVar])!=1) warning("consider using a transformation such that timeVar starts at 1")
  if (depStruct=="CI") depStruct = "UNC"
  if (!(depStruct %in% c("UNC","ARp","CS","DEC","MA1"))) stop("accepted depStruct: UNC, ARp, CS, DEC or MA1")
  #
  # Initial values
  if (is.null(control$initialValues$beta)||is.null(control$initialValues$sigma2)||is.null(control$initialValues$D)) {
    if (length(yna)>0){
      lmefit = try(lme(formFixed,random=formula(paste('~',as.character(formRandom)[length(formRandom)],'|',"ind")),data=data[-yna,]), silent=T)
      if (is(lmefit,"try-error")) {
        lmefit = try(lme(formFixed,random=~1|ind,data=data[-yna,]), silent=TRUE)
        if (is(lmefit,"try-error")) { stop("error in calculating initial values")
        } else { D1init = diag(q1)*as.numeric(var(random.effects(lmefit))) }
      } else { D1init = (var(random.effects(lmefit))) }

    } else {
      lmefit = try(lme(formFixed,random=formula(paste('~',as.character(formRandom)[length(formRandom)],'|',"ind")),data=data), silent=T)
      if (is(lmefit,"try-error")) {
        lmefit = try(lme(formFixed,random=~1|ind,data=data), silent=TRUE)
        if (is(lmefit,"try-error")) { stop("error in calculating initial values")
        } else { D1init = diag(q1)*as.numeric(var(random.effects(lmefit))) }
      } else { D1init = (var(random.effects(lmefit))) }
    }
  } # End if
  if (!is.null(control$initialValues$beta)){ beta1 = control$initialValues$beta } else { beta1 = as.numeric(lmefit$coefficients$fixed) }
  if (!is.null(control$initialValues$sigma2)){ sigmae = control$initialValues$sigma2 } else { sigmae = as.numeric(lmefit$sigma^2) }
  if (!is.null(control$initialValues$D)){ D1 = control$initialValues$D } else { D1 = D1init }
  if (length(D1)==1 && !is.matrix(D1)) D1 = as.matrix(D1)
  #
  if (length(beta1)!=p) stop ("wrong dimension of beta")
  if (!is.matrix(D1)) stop ("D must be a matrix")
  if ((ncol(D1)!=q1)|(nrow(D1)!=q1)) stop ("wrong dimension of D")
  if (length(sigmae)!=1) stop ("wrong dimension of sigma2")
  if (sigmae<=0) stop ("sigma2 must be positive")
  #
  if (distr=="t"){
    if (is.null(control$initialValues$nu)){ nu = 10 } else { nu = control$initialValues$nu }
    if (length(c(nu))>1) stop ("wrong dimension of nu")
    if (nu <=2 ) stop ("nu must be greater than 2")
  }
  #
  if (depStruct%in%c("ARp", "CS", "DEC", "MA1")){
    if (depStruct=="ARp"){
      if (!is.numeric(pAR)) stop ("pAR must be provided")
      if (pAR<=0 | pAR%%1!=0) stop ("pAR must be integer greater than 1")
    }

    if (!is.null(control$initialValues$phi)){
      phi = control$initialValues$phi
      if (depStruct=="ARp"){
        if (length(c(phi))!= pAR) stop ("initial value from phi must be in agreement with pAR")
        piis = tphitopi(phi)
        if (any(piis<(-1)|piis>1)) stop ("invalid initial value from phi")
      } else if (depStruct=="DEC"){
        if (length(c(phi))!=2) stop ("phi must be a bi-dimensional vector")
        if (any(phi<=0 | phi>=1)) stop ("phi elements must be in (0, 1)")
      } else if (depStruct=="CS"){
        if (length(c(phi))>1) stop ("phi must be a number in (0, 1)")
        if (phi<=0 | phi>=1) stop ("phi must be a number in (0, 1)")
      } else if (depStruct=="MA1"){
        if (length(c(phi))>1) stop ("phi must be a number in (-0.50, 0.50)")
        if (phi<=-0.5 | phi>=0.5) stop ("phi must be a number in (-0.50, 0.50)")
      }
    } else {
      if (length(yna)>0) ydif = (y-x%*%beta1)[-yna]
      else ydif = (y-x%*%beta1)

      if (depStruct=="ARp") phi = as.numeric(pacf(ydif, lag.max=pAR, plot=F)$acf) # pit
      if (depStruct=="CS")  phi = abs(as.numeric(pacf(ydif, lag.max=1, plot=F)$acf))
      if (depStruct=="MA1"){
        phit = seq(-0.4, 0.4, by=.05)
        logveroDECvs = function(phitheta, distri){
          if (distri=="norm") fun = loglikFuncN(y,cc,x,z,ttc,nj,lcl,ucl,beta1,sigmae,D1,phitheta,"MA1")
          if (distri=="t") fun = loglikFunct(nu,y,x,z,cc,ttc,nj,lcl,ucl,beta1,sigmae,D1,phitheta,"MA1")
          return (fun)
        }
        logverovec = unlist(lapply(phit, logveroDECvs, distri=distr))
        phi = as.numeric(phit[which.max(logverovec)])
      }
      if (depStruct=="DEC"){
        thetat = seq(0.1, .9, by=.05)
        phit = seq(0.1, .9, by=.05)
        vect = merge(phit, thetat, all=T)
        logveroDECvs = function(phitheta, distri){
          if (distri=="norm") fun = loglikFuncN(y,cc,x,z,ttc,nj,lcl,ucl,beta1,sigmae,D1,phitheta,"DEC")
          if (distri=="t") fun = loglikFunct(nu,y,x,z,cc,ttc,nj,lcl,ucl,beta1,sigmae,D1,phitheta,"DEC")
          return (fun)
        }
        logverovec = apply(vect, 1, logveroDECvs, distri=distr)
        phi = as.numeric(vect[which.max(logverovec),])
      }
    }
  } else { phi = NULL }

  # Estimation
  # ------------------------
  if (distr == "norm") obj.out = censEM.norm(y, x, z, cc, lcl, ucl, ind, ttc, data, beta1, sigmae, D1, phi,
                                             depStruct, pAR, precision=control$tol, informa=control$calc.se, itermax=control$max.iter,
                                             showiter=(!control$quiet), showerroriter=(!control$quiet)&&control$showCriterium)
  if (distr == "t")    obj.out = censEM.t(y, x, z, cc, lcl, ucl, ind, ttc, data, beta1, sigmae, D1, phi, nu, nufix,
                                          depStruct, pAR, precision=control$tol, informa=control$calc.se, itermax=control$max.iter,
                                          showiter=(!control$quiet), showerroriter=(!control$quiet)&&control$showCriterium)
  obj.out$call = match.call()
  obj.out$data = data
  obj.out$formula = list(formFixed=formFixed, formRandom=formRandom)
  obj.out$depStruct = depStruct
  obj.out$covRandom = 'pdSymm'
  obj.out$distr = distr
  obj.out$N = N
  obj.out$n = m
  obj.out$groupVar = groupVar
  obj.out$timeVar = timeVar
  #
  fitted = numeric(N)
  ind_levels = levels(ind)
  for (i in seq_along(ind_levels)) {
    seqi = ind==ind_levels[i]
    xfiti = matrix(x[seqi,], ncol=p)
    zfiti = matrix(z[seqi,], ncol=q1)
    fitted[seqi] = xfiti%*%obj.out$estimates$beta + zfiti%*%obj.out$random.effects[i]
  }
  obj.out$fitted = fitted
  #
  class(obj.out) = c("SMNCens","list")
  obj.out
}

# Fitted values
# ------------------------------------------------------------------------------
fitted.SMNCens = function(object,...) object$fitted

# Summary and print functions
# ------------------------------------------------------------------------------
print.SMNCens = function(x,...){
  cat("Linear mixed models with distribution", x$distr, "and dependency structure", x$depStruct,"\n")
  cat("Call:\n")
  print(x$call)
  cat("\nFixed: ")
  print(x$formula$formFixed)
  cat("Random:\n")
  cat("  Formula: ")
  print(x$formula$formRandom)
  cat("  Structure:", ifelse(x$covRandom=='pdSymm','General positive-definite',
                            'Diagonal'),'\n')
  cat("  Estimated variance (D):\n")
  D1 = x$estimates$D
  colnames(D1)=row.names(D1)= colnames(model.matrix(x$formula$formRandom,data=x$data))
  print(D1)
  cat("\nEstimated parameters:\n")
  if (!is.null(x$std.error)) {
    tab = round(rbind(x$theta, c(x$std.error, rep(NA, length(x$theta)-length(x$std.error)))),4)
    colnames(tab) = names(x$theta)
    rownames(tab) = c("","s.e.")
  } else {
    tab = round(rbind(x$theta),4)
    colnames(tab) = names(x$theta)
    rownames(tab) = c("")
  }
  print(tab)
  cat('\n')
  cat('Model selection criteria:\n')
  critFin = c(x$loglik, x$criteria$AIC, x$criteria$BIC, x$criteria$SIC)
  critFin = round(t(as.matrix(critFin)),digits=3)
  dimnames(critFin) = list(c(""),c("logLik", "AIC", "BIC", "SIC"))
  print(critFin)
  cat('\n')
  cat('Number of observations:',x$N,'\n')
  cat('Number of groups:',x$n,'\n')
}

summary.SMNCens = function(object, confint.level=0.95, ...){
  cat("Linear mixed models with distribution", object$distr, "and dependency structure", object$depStruct,"\n")
  cat("Call:\n")
  print(object$call)
  cat("\nDistribution", object$distr)
  if (object$distr!="norm") cat(" with nu =", object$estimates$nu, "\n")
  cat("Random effects:\n")
  cat("  Formula: ")
  print(object$formula$formRandom)
  cat("  Structure:", ifelse(object$covRandom=='pdSymm','General positive-definite',
                             'Diagonal'),'\n')
  cat("  Estimated variance (D):\n")
  D1 = object$estimates$D
  colnames(D1)=row.names(D1)= colnames(model.matrix(object$formula$formRandom,data=object$data))
  print(D1)
  cat("\nFixed effects: ")
  print(object$formula$formFixed)
  p = length(c(object$estimates$beta))
  if (!is.null(object$std.error)){
    cat("with approximate confidence intervals\n")
    qIC = qnorm(.5+confint.level/2)
    ICtab = cbind(object$estimates$beta-qIC*object$std.error[1:p],
                  object$estimates$beta+qIC*object$std.error[1:p])
    tab = (cbind(object$estimates$beta, object$std.error[1:p], ICtab))
    rownames(tab) = names(object$theta[1:p])
    colnames(tab) = c("Value", "Std.error", paste0("CI ",confint.level*100,"% lower"), paste0("CI ",confint.level*100,"% upper"))
  } else {
    cat(" (std errors not estimated)\n")
    tab = matrix(object$estimates$beta, nrow=1)
    colnames(tab) = names(object$theta[1:p])
    rownames(tab) = c("Value")
  }
  print(tab)
  cat("\nDependency structure: ", object$depStruct, "\n")
  cat("  Estimate(s):\n")
  covParam = c(object$estimates$sigma2, object$estimates$phi)
  if (object$depStruct=="UNC") names(covParam) = "sigma2"
  else names(covParam) = c("sigma2", paste0("phi",1:(length(covParam)-1)))
  print(covParam)
  cat("\nModel selection criteria:\n")
  criteria = c(object$loglik, object$criteria$AIC, object$criteria$BIC, object$criteria$SIC)
  criteria = round(t(as.matrix(criteria)), digits=3)
  dimnames(criteria) = list(c(""),c("logLik", "AIC", "BIC", "SIC"))
  print(criteria)
  cat('\n')
  cat('Number of observations:', object$N,'\n')
  cat('Number of groups:', object$n,'\n')
}

# Prediction
# ------------------------------------------------------------------------------




# Residuals
# ------------------------------------------------------------------------------
residuals.SMNCens = function(object, level="conditional", type="response",...){
  if (!(level %in% c("marginal","conditional"))) stop("Accepted levels: marginal, conditional")
  if (!(type %in% c("response","modified","normalized"))) stop("Accepted types: response, normalized or modified")
  data = object$data
  formFixed = object$formula$formFixed
  formRandom = object$formula$formRandom
  groupVar = object$groupVar
  timeVar = object$timeVar
  x = model.matrix(formFixed, data=model.frame(formFixed,data,na.action=NULL))
  y = object$yest
  z = model.matrix(formRandom, data=data)
  ind = data[,groupVar]
  if (!is.null(timeVar)){
    time = data[,timeVar]
  } else {
    time = numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] = seq_len(sum(ind==indi))
  }
  p  = ncol(x)
  q1 = ncol(z)
  N  = nrow(data)
  ind_levels = levels(ind)
  distr = object$distr
  sigmae = object$estimates$sigma2
  #
  res = numeric(N)
  if (type=="response") {
    if (level=="marginal") {
      lab = "marginal raw residuals"
      for (i in seq_along(ind_levels)) {
        seqi = ind==ind_levels[i]
        xfiti = matrix(x[seqi,], ncol=p)
        res[seqi] = y[seqi] - xfiti%*%object$estimates$beta
      }
    } else {
      lab = "conditional raw residuals"
      for (i in seq_along(ind_levels)) {
        seqi = ind==ind_levels[i]
        xfiti = matrix(x[seqi,],ncol=p)
        zfiti = matrix(z[seqi,],ncol=q1)
        res[seqi] = y[seqi] - (xfiti%*%object$estimates$beta + zfiti%*%object$random.effects[i])
      }
    }
  } # End response

  if (type=="modified"){
    if (level=="marginal") {
      lab = "marginal modified residuals"
      Dest = object$estimates$D
      for (i in seq_along(ind_levels)) {
        seqi = ind==ind_levels[i]
        xfiti = matrix(x[seqi,],ncol=p)
        zfiti = matrix(z[seqi,],ncol=q1)
        timei = time[seqi]
        Sigmaest = sigmae*MatDec(timei, object$estimates$phi, object$depStruct)
        vary = Sigmaest + (zfiti)%*%Dest%*%t(zfiti)
        sigFitinv = matrix.sqrt(solve(vary))
        res[seqi] = sigFitinv%*%(y[seqi] - xfiti%*%object$estimates$beta)
      }
    } else {
      lab = "conditional modified residuals"
      for (i in seq_along(ind_levels)) {
        seqi = ind==ind_levels[i]
        xfiti = matrix(x[seqi,],ncol=p)
        zfiti = matrix(z[seqi,],ncol=q1)
        timei = time[seqi]
        Sigmaest = sigmae*MatDec(timei, object$estimates$phi, object$depStruct)
        sigeFitinv = matrix.sqrt(solve(Sigmaest))
        res[seqi]  = sigeFitinv%*%(y[seqi] - (xfiti%*%object$estimates$beta + zfiti%*%object$random.effects[i]))
      }
    }
  } # End modified

  if (type=="normalized"){
    if (distr=="t") if (object$estimates$nu<=2) stop("normalized residual not defined for nu<=2")
    if (distr=="norm") { k2 = 1 }
    if (distr=="t") { k2 = object$estimates$nu/(object$estimates$nu-2) }

    if (level=="marginal") {
      Dest = object$estimates$D
      lab = "marginal standardized residuals"
      for (i in seq_along(ind_levels)) {
        seqi = ind==ind_levels[i]
        xfiti = matrix(x[seqi,],ncol=p)
        zfiti = matrix(z[seqi,],ncol=q1)
        timei = time[seqi]
        Sigmaest = sigmae*MatDec(timei, object$estimates$phi, object$depStruct)
        vary = Sigmaest + (zfiti)%*%Dest%*%t(zfiti)
        sigFitinv = matrix.sqrt(solve(k2*vary))
        res[seqi] = sigFitinv%*%(y[seqi]- xfiti%*%object$estimates$beta)
      }
    } else{
      lab = "conditional standardized residuals"
      for (i in seq_along(ind_levels)) {
        seqi = ind==ind_levels[i]
        xfiti = matrix(x[seqi,],ncol=p)
        zfiti = matrix(z[seqi,],ncol=q1)
        timei = time[seqi]
        Sigmaest = sigmae*MatDec(timei, object$estimates$phi, object$depStruct)
        sigeFitinv = matrix.sqrt(solve(k2*Sigmaest))
        res[seqi] = sigeFitinv%*%(y[seqi] - (xfiti%*%object$estimates$beta + zfiti%*%object$random.effects[i]))
      }
    }
  } # End normalized
  attr(res, "label") = lab
  res
}

## Plot residuals
plot.SMNCens = function(x, type="response", level="conditional", useweight=TRUE, alpha=.3,...) {
  resid = residuals(x, type=type, level=level)
  distrp = toupper(x$distr)
  if (distrp=='NORM') distrp = "N"
  depStructp = x$depStruct
  if (depStructp=="ARp") depStructp = paste0("AR(",length(x$estimates$phi),")")
  if (useweight){
    peso = data.frame(weight=x$uhat)
    peso$ind = row.names(peso)
    peso = left_join(x$data,peso,by='ind')
    peso$fitted = fitted(x)
    peso$resid = resid
    ggplot(peso, aes_string(x="fitted",y="resid",color="weight"))+geom_point()+theme_minimal() +
      geom_hline(yintercept = 0,linetype="dashed") + ylab(attr(resid,"label")) +
      xlab("fitted values") +
      scale_color_continuous(high = "#132B43", low = "#56B1F7") +
      ggtitle(paste0(depStructp,'-',distrp,'Cens-LMM')) +
      theme(plot.title = element_text( face="italic", size=10))
  } else {
    qplot(fitted(x),resid,alpha=I(alpha),...=...) +theme_minimal() +
      geom_hline(yintercept = 0,linetype="dashed") + ylab(attr(resid,"label")) +
      xlab("fitted values")+ ggtitle(paste0(depStructp,'-',distrp,'Cens-LMM')) +
      theme(plot.title = element_text( face="italic", size=10))
  }
}


###########################################################
##     Random generator from LMEM with Censored Data     ##
###########################################################
rsmnCens.lmm = function(time, ind, x, z, sigma2, D, beta, depStruct="UNC",
                        phi=NULL, distr="norm", nu=NULL, type="left", pcens=0.10,
                        LOD=NULL) {
  if (length(D)==1 && !is.matrix(D)) D = as.matrix(D)
  if (!is.matrix(D)) stop("D must be a matrix")
  if (!is.matrix(x)) x = as.matrix(x)
  if (!is.matrix(z)) z = as.matrix(z)
  if (!is.factor(ind)) ind = haven::as_factor(ind)
  N = length(c(time))
  p = length(c(beta))
  q = nrow(D)
  if (length(c(ind))!=N) stop("incompatible dimension of ind/time")
  if (nrow(x)!=N) stop("incompatible dimension of x/time")
  if (nrow(z)!=N) stop("incompatible dimension of z/time")
  if (ncol(x)!=p) stop("incompatible dimension of x/beta")
  if (ncol(z)!=q) stop("incompatible dimension of z/D")
  if (ncol(D)!=q) stop("D must be a square numeric matrix")
  if (!is.symmetric.matrix(D)) stop("D must be a symmetric matrix")
  if (!is.positive.definite(D)) stop("D must be a positive-definite matrix")
  if ((sum(is.na(x))+sum(is.na(z))+sum(is.na(ind)))>0) stop ("NAs not allowed")
  #
  if (length(sigma2)!=1) stop ("wrong dimension of sigma2")
  if (sigma2<=0) stop("sigma2 must be positive")
  if (!(distr %in% c("norm","t"))) stop("Accepted distributions: norm and t")
  if (distr=="t"){
    if (length(c(nu))>1) stop ("wrong dimension of nu")
    if (nu <=2 ) stop ("nu must be greater than 2")
  }
  #
  if (!(depStruct %in% c("UNC","ARp","CS", "CI","DEC","MA1"))) stop("accepted depStruct: UNC, ARp, CS, DEC or MA1")
  if (depStruct=="CI") depStruct = "UNC"
  if (depStruct=="ARp" && ((sum(!is.wholenumber(time))>0)||(sum(time<=0)>0))) stop("time must contain positive integer numbers when using ARp dependency")
  if (depStruct=="ARp") if (min(time)!=1) warning("consider using a transformation such that time starts at 1")
  if (depStruct%in%c("ARp", "CS", "DEC", "MA1")){
    if (is.null(phi)) stop("phi must be provided for ARp, CS, DEC, and MA1 dependency")
    if (depStruct=="ARp"){
      piis = tphitopi(phi)
      if (any(piis<(-1)|piis>1)) stop ("invalid value for phi")
    } else if (depStruct=="DEC"){
      if (length(c(phi))!=2) stop ("phi must be a bi-dimensional vector")
      if (any(phi<=0 | phi>=1)) stop ("phi elements must be in (0, 1)")
    } else if (depStruct=="CS"){
      if (length(c(phi))>1) stop ("phi must be a number in (0, 1)")
      if (phi<=0 | phi>=1) stop ("phi must be a number in (0, 1)")
    } else if (depStruct=="MA1"){
      if (length(c(phi))>1) stop ("phi must be a number in (-0.50, 0.50)")
      if (phi<=-0.5 | phi>=0.5) stop ("phi must be a number in (-0.50, 0.50)")
    }
  }
  #
  if (!(type%in%c("left","right"))) stop("accepted type: left or right")
  if (!is.null(LOD)){
    if (!is.numeric(LOD) | length(c(LOD))>1) stop("pcens must be a real number")
  } else if (!is.null(pcens)){
    if (!is.numeric(pcens) | length(c(pcens))>1) stop("pcens must be a number in [0,1]")
    if (pcens<0 | pcens>=1) stop("pcens must be a number in [0,1]")
  } else { stop("pcens or LOD must be provided") }
  #
  sample = randomCens.lmm(time,ind,x,z,sigma2,D,beta,depStruct,phi,distr,nu,pcens,LOD,type)
  return (sample)
}


# Update function: based on nlme update.lme
update.SMNCens = function (object, ..., evaluate=TRUE){
  call = object$call
  if (is.null(call))
    stop("need an object with call component")
  extras = match.call(expand.dots = FALSE)$...
  if(length(extras) > 0) {
    existing = !is.na(match(names(extras), names(call)))
    ## do these individually to allow NULL to remove entries.
    for (a in names(extras)[existing]) call[[a]] = extras[[a]]
    if(any(!existing)) {
      call = c(as.list(call), extras[!existing])
      call = as.call(call)
    }
  }
  if(evaluate) eval(call, parent.frame())
  else call
}
