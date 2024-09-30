###########################################################
##   Linear Mixed-Effects Models with Censored Response  ##
###########################################################

# Censored mixed-effects models for irregularly observed repeated measures
# ------------------------------------------------------------------------------
smn.clmm = function(data, formFixed, groupVar, formRandom=~1, depStruct="UNC",
                    ci, lcl, ucl, timeVar=NULL, distr="norm", nufix=FALSE,
                    pAR=1, control=lmmControl()){

  if (!is(formFixed,"formula")) stop("formFixed must be a formula")
  if (!is(formRandom,"formula")) stop("formRandom must be a formula")
  if (!inherits(control,"lmmControl")) stop("control must be a list generated with lmmControl()")
  if (!is.logical(nufix)) stop ("nufix must be TRUE or FALSE")
  #
  if (!is.character(groupVar)) stop("groupVar must be a character containing the name of the grouping variable in data")
  if (!is.null(timeVar)&&!is.character(timeVar)) stop("timeVar must be a character containing the name of the time variable in data")
  if (!missing(ci)&&!is.character(ci)) stop("ci must be a character containing the name of the censoring indicator in data")
  if (!missing(lcl)&&!is.character(lcl)) stop("lcl must be a character containing the name of the lower censoring limit in data")
  if (!missing(ucl)&&!is.character(ucl)) stop("ucl must be a character containing the name of the upper censoring limit in data")
  if (missing(ci)&&!missing(lcl)) stop("ci must be specified when lcl or ucl is used")
  if (missing(ci)&&!missing(ucl)) stop("ci must be specified when lcl or ucl is used")
  if (length(formFixed)!=3) stop("formFixed must be a two-sided linear formula object")
  if (!is.data.frame(data)) stop("data must be a data.frame")
  if (length(class(data))>1) data=as.data.frame(data)
  if (missing(ci)) ci = NULL
  if (missing(lcl)) lcl = NULL
  if (missing(ucl)) ucl = NULL
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
  if (sum(cc%in%c(0,1)) < length(cc)) stop("The elements of ci must contain only 0 or 1")
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
  if (!(depStruct %in% c("UNC","ARp","CS","DEC","MA1", "CAR1"))) stop("accepted depStruct: UNC, ARp, CS, DEC, CAR1 or MA1")
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
  if (depStruct%in%c("ARp", "CS", "DEC", "MA1", "CAR1")){
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
      } else if (depStruct=="CAR1"){
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
      if (depStruct=="CAR1") phi = abs(as.numeric(pacf(ydif, lag.max=1, plot=F)$acf))
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
  obj.out$ncens = sum(cc)
  obj.out$groupVar = groupVar
  obj.out$timeVar = timeVar
  #
  fitted = numeric(N)
  ind_levels = levels(ind)
  for (i in seq_along(ind_levels)) {
    seqi = ind==ind_levels[i]
    xfiti = matrix(x[seqi,], ncol=p)
    zfiti = matrix(z[seqi,], ncol=q1)
    fitted[seqi] = xfiti%*%obj.out$estimates$beta + zfiti%*%obj.out$random.effects[i,]
  }
  obj.out$fitted = fitted
  obj.out$estimates$beta <- as.numeric(obj.out$estimates$beta)
  names(obj.out$estimates$beta) <- colnames(x)
  #
  class(obj.out) = c("SMNclmm","list")
  obj.out
}

# Fitted values
# ------------------------------------------------------------------------------
fitted.SMNclmm = function(object,...) object$fitted

# Summary and print functions
# ------------------------------------------------------------------------------
print.SMNclmm = function(x, confint.level=0.95, ...){
  cat("Censored linear mixed models with distribution", x$distr, "and dependency structure", x$depStruct,"\n")
  cat("Call:\n")
  print(x$call)
  cat("\nDistribution", x$distr)
  if (x$distr!="norm") cat(" with nu =", x$estimates$nu)
  cat("\n")
  cat("\nRandom effects:\n")
  cat("  Formula: ")
  print(x$formula$formRandom)
  cat("  Structure:", ifelse(x$covRandom=='pdSymm','General positive-definite',
                             'Diagonal'),'\n')
  cat("  Estimated variance (D):\n")
  D1 = x$estimates$D
  colnames(D1)=row.names(D1)= colnames(model.matrix(x$formula$formRandom,data=x$data))
  print(D1)
  cat("\nFixed effects: ")
  print(x$formula$formFixed)
  p = length(c(x$estimates$beta))
  if (!is.null(x$std.error)){
    cat("with approximate confidence intervals\n")
    qIC = qnorm(.5+confint.level/2)
    ICtab = cbind(x$estimates$beta - qIC*x$std.error[1:p],
                  x$estimates$beta + qIC*x$std.error[1:p])
    tab = (cbind(x$estimates$beta, x$std.error[1:p], ICtab))
    rownames(tab) = names(x$theta[1:p])
    colnames(tab) = c("Value", "Std.error", paste0("CI ",confint.level*100,"% lower"), paste0("CI ",confint.level*100,"% upper"))
  } else {
    cat(" (std errors not estimated)\n")
    tab = matrix(x$estimates$beta, nrow=1)
    colnames(tab) = names(x$theta[1:p])
    rownames(tab) = c("Value")
  }
  print(tab)
  cat("\nDependency structure: ", x$depStruct, "\n")
  cat("  Estimate(s):\n")
  covParam = c(x$estimates$sigma2, x$estimates$phi)
  if (x$depStruct=="UNC") names(covParam) = "sigma2"
  else names(covParam) = c("sigma2", paste0("phi",1:(length(covParam)-1)))
  print(covParam)
  cat("\nModel selection criteria:\n")
  criteria = c(x$loglik, x$criteria$AIC, x$criteria$BIC)
  criteria = round(t(as.matrix(criteria)), digits=3)
  dimnames(criteria) = list(c(""),c("logLik", "AIC", "BIC"))
  print(criteria)
  cat('\n')
  cat('Number of observations:', x$N,'\n')
  cat('Number of censored/missing observations:', x$ncens,'\n')
  cat('Number of groups:', x$n,'\n')
}

# summary.SMNclmm = function(object, confint.level=0.95, ...){
#   cat("Censored linear mixed models with distribution", object$distr, "and dependency structure", object$depStruct,"\n")
#   cat("Call:\n")
#   print(object$call)
#   cat("\nDistribution", object$distr)
#   if (object$distr!="norm") cat(" with nu =", object$estimates$nu)
#   cat("\n")
#   cat("\nRandom effects:\n")
#   cat("  Formula: ")
#   print(object$formula$formRandom)
#   cat("  Structure:", ifelse(object$covRandom=='pdSymm','General positive-definite',
#                              'Diagonal'),'\n')
#   cat("  Estimated variance (D):\n")
#   D1 = object$estimates$D
#   colnames(D1)=row.names(D1)= colnames(model.matrix(object$formula$formRandom,data=object$data))
#   print(D1)
#   cat("\nFixed effects: ")
#   print(object$formula$formFixed)
#   p = length(c(object$estimates$beta))
#   if (!is.null(object$std.error)){
#     cat("with approximate confidence intervals\n")
#     qIC = qnorm(.5+confint.level/2)
#     ICtab = cbind(object$estimates$beta-qIC*object$std.error[1:p],
#                   object$estimates$beta+qIC*object$std.error[1:p])
#     tab = (cbind(object$estimates$beta, object$std.error[1:p], ICtab))
#     rownames(tab) = names(object$theta[1:p])
#     colnames(tab) = c("Value", "Std.error", paste0("CI ",confint.level*100,"% lower"), paste0("CI ",confint.level*100,"% upper"))
#   } else {
#     cat(" (std errors not estimated)\n")
#     tab = matrix(object$estimates$beta, nrow=1)
#     colnames(tab) = names(object$theta[1:p])
#     rownames(tab) = c("Value")
#   }
#   print(tab)
#   cat("\nDependency structure: ", object$depStruct, "\n")
#   cat("  Estimate(s):\n")
#   covParam = c(object$estimates$sigma2, object$estimates$phi)
#   if (object$depStruct=="UNC") names(covParam) = "sigma2"
#   else names(covParam) = c("sigma2", paste0("phi",1:(length(covParam)-1)))
#   print(covParam)
#   cat("\nModel selection criteria:\n")
#   criteria = c(object$loglik, object$criteria$AIC, object$criteria$BIC)
#   criteria = round(t(as.matrix(criteria)), digits=3)
#   dimnames(criteria) = list(c(""),c("logLik", "AIC", "BIC"))
#   print(criteria)
#   cat('\n')
#   cat('Number of observations:', object$N,'\n')
#   cat('Number of censored/missing observations:', object$ncens,'\n')
#   cat('Number of groups:', object$n,'\n')
# }
summary.SMNclmm <- function(object, confint.level=.95, ...){
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
                estimates = object$estimates, n = object$n, N = object$N, ncens = object$ncens,
                covParam = covParam)
  class(outobj) <- c("SMNcenssumm","list")
  outobj
}

print.SMNcenssumm <- function(x,...){
  cat("Censored linear mixed models with distribution", x$distr, "and dependency structure", x$depStruct,"\n")
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
  cat("\nDependency structure:", x$depStruct)
  cat("\n  Estimate(s):\n")
  print(x$covParam)
  cat("\nSkewness parameter estimate:", x$estimates$lambda)
  cat('\n')
  cat('\nModel selection criteria:\n')
  print(x$criteria)
  cat('\n')
  cat('Number of observations:',x$N,'\n')
  cat('Number of censored/missing observations:', x$ncens,'\n')
  cat('Number of groups:',x$n,'\n')
}


# Mahalanobis distance
# ------------------------------------------------------------------------------
mahalDistCens = function(object){
  if(!is(object, "SMNclmm")) stop("object must inherit from class SMNclmm")
  formFixed  = object$formula$formFixed
  formRandom = object$formula$formRandom
  groupVar = object$groupVar
  timeVar  = object$timeVar
  data = object$data
  x = model.matrix(formFixed, data=model.frame(formFixed,data,na.action=NULL))
  y = object$yest
  z = model.matrix(formRandom, data=data)
  ind = data[,groupVar]
  if (!is.null(timeVar)) {
    time = data[,timeVar]
  } else{
    time = numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] = seq_len(sum(ind==indi))
  }
  p  = ncol(x)
  q1 = ncol(z)
  N  = nrow(data)
  ind_levels = levels(ind)
  #
  distr <- object$distr
  mahaldist = distbi = distei = numeric(length(ind_levels))
  Dest = object$estimates$D
  for (i in seq_along(ind_levels)) {
    seqi = (ind==ind_levels[i])
    xfiti = x[seqi,]
    zfiti = z[seqi,]
    timei = time[seqi]
    Sigmaest = object$estimates$sigma2*MatDec(timei, object$estimates$phi, object$depStruct)
    Psiy = Sigmaest + (zfiti)%*%Dest%*%t(zfiti)
    #
    ytil = y[seqi] - xfiti%*%object$estimates$beta
    mahaldist[i] = t(ytil)%*%solve(Psiy)%*%ytil
  }
  #
  out = mahaldist
  names(out) = row.names(object$random.effects)
  class(out) = c("mahalDistCens","numeric")
  attr(out,'call') = match.call()
  out
}

plot.mahalDistCens = function(x, fitobject, level=.99, nlabels=3,...){
  if (missing(fitobject)) fitobject = eval(str2lang(as.character(attr(x,'call')[2])))
  if (!is(x, "mahalDistCens")) stop("x must inherit from class mahalDistCens")
  if (!is.data.frame(x)) x = data.frame(md=x)
  if (level>=1|level<=0) stop("0<level<1 needed")
  #
  if(!is(fitobject,"SMNclmm")) stop("fitobject must inherit from class SMNclmm")
  data = fitobject$data
  timeVar = fitobject$timeVar
  ind = data[,fitobject$groupVar]
  if (!is.null(timeVar)) {
    time = data[,timeVar]
  } else{
    time = numeric(length = length(ind))
    for (indi in levels(ind)) time[ind==indi] = seq_len(sum(ind==indi))
  }
  distr = fitobject$distr
  nu = fitobject$estimates$nu
  #
  x$nj = tapply(time, ind, length)
  x$index = seq_along(x$nj)
  x$ind = levels(ind)
  #
  distrp = toupper(distr)
  if (distrp=='NORM') distrp = "N"
  depStructp = fitobject$depStruct
  if (depStructp=="ARp") depStructp = paste0("AR(",length(fitobject$estimates$phi),")")
  #
  if (n_distinct(x$nj) == 1) {
    nj1 = x$nj[1]
    if (distr=="norm") mdquantile = qchisq(level, nj1)
    if (distr=="t") mdquantile = nj1*qf(level,nj1,nu)
    plotout = ggplot(x, aes_string("index","md")) +
      geom_point(shape=1) + ylab("Mahalanobis distance") +
      geom_text_repel(aes_string(label="ind"),data=subset(x,rank(x$md)>length(x$nj)-nlabels),
                      nudge_x=.5, nudge_y=.5, size=3) +
      geom_hline(yintercept=mdquantile, col=4, linetype="dashed")
    attr(plotout,"info") = data.frame(nj=nj1, quantile=c(mdquantile))
  } else {
    njvec = sort(unique(x$nj))
    if (distr=="norm") mdquantile = qchisq(level,njvec)
    if (distr=="t") mdquantile = njvec*qf(level,njvec,nu)
    datline = data.frame(nj=c(njvec,max(njvec)+1), quantile=c(mdquantile,max(mdquantile)))
    datline$nj2 = datline$nj-.5
    plotout = ggplot(x,aes_string("nj","md")) +
      geom_point(position=position_jitter(width=.25, height=0),
                 shape=1) + ylab("Mahalanobis distance") + xlab("number of observations") +
      geom_text(aes_string(label="ind"), data=subset(x,rank(x$md)>length(x$nj)-nlabels),
                nudge_x=0, nudge_y=0, size=3) +
      geom_step(aes_string(x="nj2",y="quantile"), data=datline, color=4, linetype="dashed") +
      scale_x_continuous(breaks=njvec)
    attr(plotout,"info") = data.frame(nj=c(njvec), quantile=c(mdquantile))
  }
  plotout + theme_minimal() + ggtitle(paste0(depStructp,'-',distrp,'-CLMM')) +
    theme(plot.title=element_text(face="italic", size=10))
}


# Prediction
# ------------------------------------------------------------------------------
predict.SMNclmm = function(object, newData,...){
  if (missing(newData)||is.null(newData)) return(fitted(object))
  if (!is.data.frame(newData)) stop("newData must be a data.frame object")
  if (nrow(newData)==0) stop("newData can not be an empty dataset")
  dataFit = object$data
  formFixed = object$formula$formFixed
  formRandom = object$formula$formRandom
  groupVar = object$groupVar
  timeVar = object$timeVar
  dataPred = newData
  #
  vars_used = unique(c(all.vars(formFixed)[-1],all.vars(formRandom),groupVar,timeVar))
  vars_miss = which(!(vars_used %in% names(newData)))
  if (length(vars_miss)>0) stop(paste(vars_used[vars_miss],"not found in newData"))
  depStruct = object$depStruct
  if (depStruct=="CI") depStruct = "UNC"
  if (any(!(dataPred[,groupVar] %in% dataFit[,groupVar]))) stop("subjects for which future values should be predicted must also be at fitting data")
  if (!is.factor(dataFit[,groupVar])) dataFit[,groupVar] = haven::as_factor(dataFit[,groupVar])
  if (!is.factor(dataPred[,groupVar])) dataPred[,groupVar]= factor(dataPred[,groupVar],levels=levels(dataFit[,groupVar]))
  #
  obj.out = predictionSMNCens(formFixed, formRandom, dataFit, dataPred, groupVar, timeVar, object$estimates, object$yest, depStruct)
  obj.out
}

# Residuals
# ------------------------------------------------------------------------------
residuals.SMNclmm = function(object, level="conditional",...){
  if (!(level %in% c("marginal","conditional"))) stop("Accepted levels: marginal, conditional")
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
      res[seqi] = y[seqi] - (xfiti%*%object$estimates$beta + zfiti%*%object$random.effects[i,])
    }
  }
  attr(res, "label") = lab
  res
}

## Plot residuals
plot.SMNclmm = function(x, level="conditional", useweight=TRUE, alpha=.3,...) {
  resid = residuals(x, level=level)
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
    ggplot(peso, aes_string(x="fitted",y="resid",color="weight")) + geom_point() + theme_minimal() +
      geom_hline(yintercept=0, linetype="dashed") + ylab(attr(resid,"label")) +
      xlab("fitted values") +
      scale_color_continuous(high="#132B43", low="#56B1F7") +
      ggtitle(paste0(depStructp, '-', distrp, '-CLMM')) +
      theme(plot.title=element_text( face="italic", size=10))
  } else {
    qplot(fitted(x),resid,alpha=I(alpha),...=...) + theme_minimal() +
      geom_hline(yintercept=0,linetype="dashed") + ylab(attr(resid,"label")) +
      xlab("fitted values") + ggtitle(paste0(depStructp, '-', distrp, '-CLMM')) +
      theme(plot.title=element_text( face="italic", size=10))
  }
}


# Update function: based on nlme update.lme
# ------------------------------------------------------------------------------
update.SMNclmm = function (object, ..., evaluate=TRUE){
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


###########################################################
##     Random generator from LMEM with Censored Data     ##
###########################################################
rsmsn.clmm = function(time, ind, x, z, sigma2, D, beta, lambda=rep(0, nrow(D)),
                      depStruct="UNC", phi=NULL, distr="norm", nu=NULL, type="left",
                      pcens=0.10, LOD=NULL) {
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
  if (length(lambda)!=q) stop ("incompatible dimension of lambda/D")
  if (ncol(D)!=q) stop("D must be a square numeric matrix")
  if (!is.symmetric.matrix(D)) stop("D must be a symmetric matrix")
  if (!is.positive.definite(D)) stop("D must be a positive-definite matrix")
  if ((sum(is.na(x))+sum(is.na(z))+sum(is.na(ind)))>0) stop ("NAs not allowed")
  #
  if (length(sigma2)!=1) stop ("wrong dimension of sigma2")
  if (sigma2<=0) stop("sigma2 must be positive")
  if (!(distr %in% c("norm","t","sn","st"))) stop("Accepted distributions: norm, t, sn, and st")
  if (distr%in%c("t", "st")){
    if (is.null(nu)) stop ("nu must be provided for t and st distributions")
    if (length(c(nu))>1) stop ("wrong dimension of nu")
    if (nu <=2 ) stop ("nu must be greater than 2")
  }
  #
  if (!(depStruct %in% c("UNC","ARp","CS", "CI","DEC","MA1", "CAR1"))) stop("accepted depStruct: UNC, ARp, CS, DEC, CAR1 or MA1")
  if (depStruct=="CI") depStruct = "UNC"
  if (depStruct=="ARp" && ((sum(!is.wholenumber(time))>0)||(sum(time<=0)>0))) stop("time must contain positive integer numbers when using ARp dependency")
  if (depStruct=="ARp") if (min(time)!=1) warning("consider using a transformation such that time starts at 1")
  if (depStruct%in%c("ARp", "CS", "DEC", "MA1", "CAR1")){
    if (is.null(phi)) stop("phi must be provided for ARp, CS, DEC, CAR1 and MA1 dependency")
    if (depStruct=="ARp"){
      piis = tphitopi(phi)
      if (any(piis<(-1)|piis>1)) stop ("invalid value for phi")
    } else if (depStruct=="DEC"){
      if (length(c(phi))!=2) stop ("phi must be a bi-dimensional vector")
      if (any(phi<=0 | phi>=1)) stop ("phi elements must be in (0, 1)")
    } else if (depStruct%in%c("CS", "CAR1")){
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
  sample = randomCens.lmm(time,ind,x,z,sigma2,D,beta,lambda,depStruct,phi,distr,nu,pcens,LOD,type)
  return (sample)
}
