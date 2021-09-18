#daarem codes from daarem package version 0.4.1 
#modified to show some information at each iteration
daarem <- function(par, fixptfn, objfn, ..., control=list(),showiter,showerroriter) {
  
  control.default <- list(maxiter=2000, order=5, tol=1.e-08, mon.tol=0.01, 
                          cycl.mon.tol=0.0, kappa=40, alpha=1.3)
  # control.default <- list(maxiter=2000, order=5, tol=1.e-08, mon.tol=0.01, 
  #                         cycl.mon.tol=0.0, kappa=25, alpha=1.2)
  namc <- names(control)
  if (!all(namc %in% names(control.default))) {
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  control <- modifyList(control.default, control)
  
  maxiter <- control$maxiter
  tol <- control$tol
  mon.tol <- control$mon.tol  ## monotonicity tolerance
  cycl.mon.tol <- control$cycl.mon.tol
  a1 <- control$alpha
  kappa <- control$kappa
  
  num.params <- length(par)
  order.default <- ifelse(num.params>20,10,floor(num.params/2))
  nlag <- min(control$order, order.default)
  #nlag <- min(control$order, ceiling(num.params/2))
  
  #if(!missing(objfn)) {
  ans <- daarem_base_objfn(par, fixptfn, objfn, maxiter, tol, mon.tol, 
                          cycl.mon.tol, a1, kappa, num.params, nlag,
                          showiter, showerroriter, ...)
  return(ans)
}
#

daarem_base_objfn <- function(par, fixptfn, objfn, maxiter, tol, mon.tol, 
                              cycl.mon.tol, a1, kappa, num.params, nlag,
                              showiter, showerroriter, ...) {
  
  #num.params <- length(par)
  #nlag <- min(control$order, ceiling(num.params/2))
  
  Fdiff <- Xdiff <- matrix(0.0, nrow=num.params, ncol=nlag)
  obj_funvals <- rep(NA, maxiter + 2)
  
  xold <- par
  xnew <- fixptfn(xold, ...)
  obj_funvals[1] <- objfn(xold, ...)
  obj_funvals[2] <- objfn(xnew, ...)
  likchg <- obj_funvals[2] - obj_funvals[1]
  obj.evals <- 2
  
  fold <- xnew - xold
  k <- 0
  count <- 0
  shrink.count <- 0
  shrink.target <- 1/(1 + a1^kappa)
  lambda.ridge <- 100000
  r.penalty <- 0
  conv <- FALSE
  ss.resids <- 10
  ell.star <- obj_funvals[2]
  while(k < maxiter) {
    k <- k+1
    count <- count + 1
    
    fnew <- fixptfn(xnew, ...) - xnew
    #ss.resids <- sqrt(crossprod(fnew))
    # if (k > 2){
    #   at<- (obj_funvals[k+1]-obj_funvals[k])/(obj_funvals[k]-obj_funvals[k-1])
    #   ss.resids <- abs((obj_funvals[k+1]-obj_funvals[k])/(1-at))
    # } #mudando para outro criterio
    if (k>1) ss.resids <- abs((obj_funvals[k+1]-obj_funvals[k])/obj_funvals[k])
    if (is.nan(ss.resids)) ss.resids=10
    #
    if (showiter&&!showerroriter) cat("Iteration ",k,"\r")#," of ",maxiter,"\r") 
    if (showerroriter) cat("Iteration ",k," of ",maxiter," - criterium =",ss.resids,
                           " - loglik =",obj_funvals[k+1],"\r")
    #
    if(ss.resids < tol) { #& count==nlag
      conv <- TRUE
      break
    }
    Fdiff[,count] <- fnew - fold
    Xdiff[,count] <- xnew - xold
    
    np <- count
    if(np==1) {
      Ftmp <- matrix(Fdiff[,1], nrow=num.params, ncol=np)
      Xtmp <- matrix(Xdiff[,1], nrow=num.params, ncol=np)  ## is this matrix function needed?
      
    } else {
      Ftmp <- Fdiff[,1:np]
      Xtmp <- Xdiff[,1:np]  
    }
    tmp <- La.svd(Ftmp)
    dvec <- tmp$d
    dvec.sq <- dvec*dvec
    uy <- crossprod(tmp$u, fnew)
    uy.sq <- uy*uy
    
    ### Still need to compute Ftf
    Ftf <- sqrt(sum(uy.sq*dvec.sq))
    tmp_lam <- DampingFind(uy.sq, dvec, a1, kappa, shrink.count, Ftf, lambda.start=lambda.ridge, r.start=r.penalty)
    lambda.ridge <- tmp_lam$lambda
    r.penalty <- tmp_lam$rr
    dd <- (dvec*uy)/(dvec.sq + lambda.ridge)
    gamma_vec <- crossprod(tmp$vt, dd)
    
    xbar <- xnew - drop(Xtmp%*%gamma_vec)
    fbar <- fnew - drop(Ftmp%*%gamma_vec)
    
    x.propose <- xbar + fbar
    new.objective.val <- try(objfn(x.propose, ...), silent=TRUE)
    obj.evals <- obj.evals + 1
    
    if(class(new.objective.val)[1] != "try-error" && !is.na(obj_funvals[k+1]) &&
       !is.nan(new.objective.val)) {
      if(new.objective.val >= obj_funvals[k+1] - mon.tol) {  ## just change this line in daarem_base_noobjfn
        ## Increase delta
        obj_funvals[k+2] <- new.objective.val
        fold <- fnew
        xold <- xnew
        
        xnew <- x.propose
        shrink.count <- shrink.count + 1
      } else {
        ## Keep delta the same
        fold <- fnew
        xold <- xnew
        
        xnew <- fold + xold
        obj_funvals[k+2] <- objfn(xnew, ...)
        obj.evals <- obj.evals + 1
      }
    } else {
      ## Keep delta the same
      fold <- fnew
      xold <- xnew
      
      xnew <- fold + xold
      obj_funvals[k+2] <- objfn(xnew, ...)
      obj.evals <- obj.evals + 1
      count <- 0
    }
    if(count==nlag) {
      count <- 0
      ## restart count
      ## make comparison here l.star vs. obj_funvals[k+2]
      if(obj_funvals[k+2] < ell.star - cycl.mon.tol) {
        ## Decrease delta
        shrink.count <- max(shrink.count - nlag, -2*kappa)
      }
      ell.star <- obj_funvals[k+2]
    }
    shrink.target <-  1/(1 + a1^(kappa - shrink.count))
  }
  if (showiter||showerroriter) cat("\n")
  obj_funvals <- obj_funvals[!is.na(obj_funvals)]
  value.obj <- objfn(xnew, ...)
  # if(k >= maxiter) {
  #   conv <- FALSE
  # }
  return(list(par=c(xnew), fpevals = k, value.objfn=value.obj, objfevals=obj.evals, 
              criterio=as.numeric(ss.resids), convergence=conv, objfn.track=obj_funvals[-1]))
}
#

fpiter <- function(par, fixptfn, objfn=NULL, control=list( ), showiter, showerroriter, ...){
  
  control.default <- list(tol=1.e-07, maxiter=5000)
  namc <- names(control)
  if (!all(namc %in% names(control.default))){
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  ctrl <- modifyList(control.default, control)
  
  tol <- ctrl$tol
  maxiter <- ctrl$maxiter
  
  iter <- 1
  #resid <- rep(NA,1)
  objeval <- 0
  conv <- FALSE
  obj_funvals <- rep(NA, maxiter)
  res <- 10
  
  while (iter <= maxiter) {
    #print(par)
    p.new <- fixptfn(par, ...)
    #p.new <- fixptfn(par, y=y,x=x,z=z,time=time,ind=ind,distr=distr,pAR=pAR,lb=lb,lu=lu)
    par <- p.new
    obj_funvals[iter] <- objfn(par, ...)
    #obj_funvals[iter] <- objfn(par, y=y,x=x,z=z,time=time,ind=ind,distr=distr,pAR=pAR,lb=lb,lu=lu)
    # if (iter > 3){
    #   at<- (obj_funvals[iter]-obj_funvals[iter-1])/(obj_funvals[iter]-obj_funvals[iter-2])
    #   res <- abs((obj_funvals[iter]-obj_funvals[iter-1])/(1-at))
    # } #mudar para abs((llji-llj1)/llj1) ?
    if (iter>2) res <- abs((obj_funvals[iter]-obj_funvals[iter-1])/obj_funvals[iter-1])
    if (is.nan(res)) res=10
    #res <- sqrt(crossprod(p.new - par))
    if (showiter&&!showerroriter) cat("Iteration ",iter,"\r")#," of ",maxiter,"\r") 
    #if (showerroriter) cat("Iteration ",iter," of ",maxiter," - criterium =",res,"\r")
    if (showerroriter) cat("Iteration ",iter," of ",maxiter," - criterium =",res,
                           " - loglik =",obj_funvals[iter],"\r")
    if ( res < tol) {conv <- TRUE; break}
    iter <- iter+1
  }
  if (showiter||showerroriter) cat("\n")
  
  loglik.best <- obj_funvals[iter-1]#objfn(par, ...)
  obj_funvals <- obj_funvals[!is.na(obj_funvals)]
  
  return(list(par=par, value.objfn=loglik.best, fpevals=iter-1, criterio=as.numeric(res), 
              convergence=conv,objfn.track=obj_funvals))
}
#

DampingFind <- function(uy.sq, dvec, aa, kappa, sk, Ftf, lambda.start=NULL,
                        r.start=NULL, maxit=10) {
  ## uy.sq is a vector whose jth component (U^t*y)^2
  ## d.sq is a vector whose jth component is d_j^2.
  ## (aa, kappa, sk) - parameters in determining delta_k
  ## Ftf <- F_{k}^T f_{k}
  ##
  ## Output: value of penalty term such that (approximately)
  ## .       || beta.ridge ||^2/||beta.ols||^2 = delta_{k}^{2}, with 0 < target < 1.
  
  if(is.null(lambda.start) | is.na(lambda.start)) {
    lambda.start <- 100
  }
  if(is.null(r.start) | is.na(r.start)) {
    r.start <- 0
  }
  if(sum(dvec==0.0) > 0) {
    ind <- dvec > 0
    dvec <- dvec[ind]
    uy.sq <- uy.sq[ind]
    ## what about if dvec has length 1?
  }
  pow <- kappa - sk
  target <- exp(-0.5*log1p(aa^pow))  ### This is sqrt(delta_k)
  #d.sq <- pmax(dvec*dvec, 1e-7)  ## in case the least-squares problem is very poorly conditioned
  d.sq <- dvec*dvec
  betahat.ls <- uy.sq/d.sq
  betahat.ls.norm <- sqrt(sum(betahat.ls))
  vk <- target*betahat.ls.norm
  
  if(vk == 0) {
    ## the norm of the betas is zero, so the value of lambda shouldn't matter
    return(list(lambda=lambda.start, rr=r.start))
  }
  ### Initialize lambda and lower and upper bounds
  lambda <- lambda.start - r.start/vk
  LL <- (betahat.ls.norm*(betahat.ls.norm - vk))/sum(uy.sq/(d.sq*d.sq))
  UU <- Ftf/vk
  
  ## Compute lstop and ustop
  pow.low <- pow + 0.5
  pow.up <- pow - 0.5
  l.stop <- exp(-0.5*log1p(aa^pow.low))
  u.stop <- exp(-0.5*log1p(aa^pow.up))
  
  ### Start iterations
  for(k in 1:maxit) {
    ### Check if lambda is inside lower and upper bounds
    if(lambda <= LL | lambda >= UU) {
      lambda <- max(.0001*UU, sqrt(LL*UU))
    }
    ### Evaluate ||s(lambda)|| and \phi(lambda)/phi'(lambda)
    
    d.lambda <- (dvec/(d.sq + lambda))^2
    d.prime <- d.lambda/(d.sq + lambda)
    s.norm <- sqrt(sum(uy.sq*d.lambda))
    phi.val <- s.norm - vk
    phi.der <- (-1)*sum(uy.sq*d.prime)/s.norm
    phi.ratio <- phi.val/phi.der
    
    d.u <- (dvec/(d.sq + LL))^2
    s.up <- sqrt(sum(uy.sq*d.u))
    ### Initial Lower bound is not correct
    
    ### Check convergence
    if(s.norm <= u.stop*betahat.ls.norm & s.norm >= l.stop*betahat.ls.norm) {
      break
    }
    
    ### If not converged, update lower and upper bounds
    UU <- ifelse(phi.val >= 0, UU, lambda)
    LL <- max(LL, lambda - phi.ratio)
    
    ### Now update lambda
    lambda <- lambda - (s.norm*phi.ratio)/vk
    bb <- (s.norm*phi.ratio)/vk
  }
  ans <- list(lambda=lambda, rr=s.norm*phi.ratio)
  return(ans)
}
