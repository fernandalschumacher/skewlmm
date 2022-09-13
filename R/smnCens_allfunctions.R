
# Auxiliary functions
# ------------------------------------------
# Correlation matrix for the error term
MatDec = function(tt, phi, struc){
  r = length(tt)

  if(struc=="UNC"){
    V = diag(1, nrow=r, ncol=r)

  } else if (struc=="DEC"){
    if (phi[2]<=0.0000001){
      V = matrix(phi[1], nrow=r, ncol=r)
      diag(V) = 1
    } else {
      H = (abs(outer(tt, tt, "-")))^phi[2]
      V = (phi[1]^H)
    }

  } else if (struc=="CS"){
    V = matrix(phi, nrow=r, ncol=r)
    diag(V) = 1

  } else if (struc=="ARp"){
    n = max(tt)
    if (n==1) Rn = matrix(1)
    else Rn = toeplitz(ARMAacf(ar=phi, ma=0, lag.max=n-1))
    rhos = ARMAacf(ar=phi, ma=0, lag.max=length(phi))[-1]
    Rn = Rn/(1 - sum(rhos*phi))
    V = as.matrix(Rn[tt,tt])

  } else if (struc=="MA1"){
    W = matrix(0, nrow=r, ncol=r)
    for (i in 1:r){
      W[i,i] = 1
      for(j in i:r){
        dif = abs(tt[i] - tt[j])
        if(dif==1){W[i,j] = W[j,i] = phi}
        }
      }
    V = W
  }
  return(V)
}

###########################################################
##     Random generator from LMEM with Censored Data     ##
###########################################################
randomCens.lmm = function(time, ind, x, z, sigmae, D1, beta, struc, phi, distr, nu,
                          pcens, lod, type){
  N  = length(c(time))
  nj = tapply(ind, ind, length)[unique(ind)]
  m  = length(c(nj))
  q1 = ncol(z)
  p  = ncol(x)
  yobs = numeric(N)
  #
  for (j in 1:m){
    a1 = sum(nj[1:j-1])+1; b1 = sum(nj[1:j])
    x1  = matrix(x[a1:b1, ], ncol=p)
    z1  = matrix(z[a1:b1, ], ncol=q1)
    tt1 = time[a1:b1]
    #
    Gamma = sigmae*MatDec(tt1, phi, struc)
    if (distr=="norm"){
      bi  = matrix(rmvnorm(n=1, mean=rep(0,q1), sigma=D1), nrow=q1)
      xii = matrix(rmvnorm(n=1, mean=rep(0,nj[j], sigma=Gamma)), nrow=nj[j])
    } else if (distr=="t"){
      bi  = matrix(rmvt(n=1, sigma=D1, df=nu, delta=rep(0,q1), type="shifted"), nrow=q1)
      xii = matrix(rmvt(n=1, sigma=Gamma, df=nu, delta=rep(0,nj[j]), type="shifted"), nrow=nj[j])
    }
    yobs[a1:b1] = x1%*%beta + z1%*%bi + xii
  }
  #
  if (all(x[,1]==1)) Xi = x[,-1]
  if (all(z[,1]==1)) Zi = z[,-1]
  #
  if (!is.null(lod)){
    if (type=="left"){
      cc  = (yobs<lod) + 0
      y   = yobs*(1-cc) + lod*cc
      lcl = rep(-Inf, N)
      ucl = rep(lod, N)
    } else if (type=="right"){
      cc  = (yobs>lod) + 0
      y   = yobs*(1-cc) + lod*cc
      lcl = rep(lod, N)
      ucl = rep(Inf, N)
    }
  } else if (!is.null(pcens)){
    if (pcens==0) return (data.frame(time=time, ind=ind, y=yobs, x=Xi, z=Zi))
    if (pcens>0){
      if (type=="left"){
        cte = as.numeric(quantile(yobs, probs=pcens))
        cc  = (yobs<cte) + 0
        y   = yobs*(1-cc) + cte*cc
        lcl = rep(-Inf, N)
        ucl = rep(cte, N)
      } else if (type=="right"){
        cte = as.numeric(quantile(yobs, probs=1-pcens))
        cc  = (yobs>cte) + 0
        y   = yobs*(1-cc) + cte*cc
        lcl = rep(cte, N)
        ucl = rep(Inf, N)
      }
    }
  }
 return (data.frame(time=time, ind=ind, y=y, ci=cc, lcl=lcl, ucl=ucl, x=Xi, z=Zi))
}

###########################################################
## Linear Normal Mixed-Effects Models with Censored Data ##
###########################################################

# Auxiliary functions
# ------------------------------------------

# Estimate phi1 and phi2 (DEC model)
FCi = function(phiG, beta1, sigmae, ttc, ubi, ubbi, uybi, uyyi, uyi, x, z, nj){
  m  = length(nj)[1]
  p  = dim(x)[2]
  q1 = dim(z)[2]
  soma = 0

  for (j in 1:m ){
    a1 = sum(nj[1:j-1])+1; b1 = sum(nj[1:j]); c1 = ((j-1)*q1)+1; d1 = j*q1
    x1 = matrix(x[a1:b1, ], ncol=p)
    z1 = matrix(z[a1:b1, ], ncol=q1)
    muii = x1%*%beta1
    tt1  = ttc[a1:b1]
    #
    ub = ubi[c1:d1, j]
    ubb = ubbi[c1:d1, c1:d1]
    uyb = uybi[a1:b1, c1:d1]
    uyy = uyyi[a1:b1, a1:b1]
    uy  = uyi[a1:b1, j]
    #
    Cii = MatDec(tt1, phiG, "DEC")
    Cii = (Cii + t(Cii))/2
    if (det(Cii)<=0){A = 1} else {A = det(Cii)}
    invCii = solve(Cii)
    #
    Ai = as.vector(sum(diag(uyy%*%invCii)) -t(uy)%*%invCii%*%muii - t(muii)%*%invCii%*%uy - sum(diag(invCii%*%((uyb)%*%t(z1)))) - sum(diag(invCii%*%((uyb)%*%t(z1))))
                     + t(muii)%*%invCii%*%z1%*%ub + t(ub)%*%t(z1)%*%invCii%*%muii + t(muii)%*%invCii%*%muii + sum(diag(ubb%*%t(z1)%*%invCii%*%z1)))
    #
    soma = soma - 0.5*log(A) - (0.5/sigmae)*Ai
  }
  return(-soma)
}

# Estimate phi (CS and MA(1) model)
FCiphi1 = function(phi1, beta1, sigmae, ttc, ubi, ubbi, uybi, uyyi, uyi, x, z, nj, struc){
  m  = length(nj)[1]
  p  = dim(x)[2]
  q1 = dim(z)[2]
  soma = 0

  for (j in 1:m ){
    a1 = sum(nj[1:j-1])+1; b1 = sum(nj[1:j]); c1 = ((j-1)*q1)+1; d1 = j*q1
    x1 = matrix(x[a1:b1, ], ncol=p)
    z1 = matrix(z[a1:b1, ], ncol=q1)
    muii = x1%*%beta1
    tt1  = ttc[a1:b1]
    #
    ub = ubi[c1:d1, j]
    ubb = ubbi[c1:d1, c1:d1]
    uyb = uybi[a1:b1, c1:d1]
    uyy = uyyi[a1:b1, a1:b1]
    uy  = uyi[a1:b1, j]
    #
    Cii = MatDec(tt1, phi1, struc)
    Cii = (Cii + t(Cii))/2
    if (det(Cii)<=0){A = 1} else {A = det(Cii)}
    invCii = solve(Cii)
    #
    Ai = as.vector(sum(diag(uyy%*%invCii)) -t(uy)%*%invCii%*%muii - t(muii)%*%invCii%*%uy - sum(diag(invCii%*%((uyb)%*%t(z1)))) - sum(diag(invCii%*%((uyb)%*%t(z1))))
                     + t(muii)%*%invCii%*%z1%*%ub + t(ub)%*%t(z1)%*%invCii%*%muii + t(muii)%*%invCii%*%muii + sum(diag(ubb%*%t(z1)%*%invCii%*%z1)))
    #
    soma = soma - 0.5*log(A) - (0.5/sigmae)*Ai
  }
  return(-soma)
}

# Estimate phi (AR(p) model)
FCiArp = function(piis, beta1, sigmae, ttc, ubi, ubbi, uybi, uyyi, uyi, x, z, nj){
  m  = length(nj)[1]
  p  = dim(x)[2]
  q1 = dim(z)[2]
  soma = 0

  for (j in 1:m ){
    a1 = sum(nj[1:j-1])+1; b1 = sum(nj[1:j]); c1 = ((j-1)*q1)+1; d1 = j*q1
    x1 = matrix(x[a1:b1, ], ncol=p)
    z1 = matrix(z[a1:b1, ], ncol=q1)
    muii = x1%*%beta1
    tt1  = ttc[a1:b1]
    #
    ub = ubi[c1:d1, j]
    ubb = ubbi[c1:d1, c1:d1]
    uyb = uybi[a1:b1, c1:d1]
    uyy = uyyi[a1:b1, a1:b1]
    uy  = uyi[a1:b1, j]
    #
    phi2 = estphit(piis)
    Cii = MatDec(tt1, phi2, "ARp")
    Cii = (Cii + t(Cii))/2
    if (det(Cii)<=0){A = 1} else {A = det(Cii)}
    invCii = solve(Cii)
    #
    Ai = as.vector(sum(diag(uyy%*%invCii)) -t(uy)%*%invCii%*%muii - t(muii)%*%invCii%*%uy - sum(diag(invCii%*%((uyb)%*%t(z1)))) - sum(diag(invCii%*%((uyb)%*%t(z1))))
                   + t(muii)%*%invCii%*%z1%*%ub + t(ub)%*%t(z1)%*%invCii%*%muii + t(muii)%*%invCii%*%muii + sum(diag(ubb%*%t(z1)%*%invCii%*%z1)))
    #
    soma = soma - 0.5*log(A) - (0.5/sigmae)*Ai
  }
  return(-soma)
}


# Log-likelihood function
# ------------------------------------------
loglikFuncN = function(y, cc, x, z, ttc, nj, LL, LU, betas, sigmae, D1, phi1, struc){
  m  = length(nj)
  p  = ncol(x)
  q1 = ncol(z)
  ver = matrix(0,m,1)

  for(j in 1:m){
    a1  = sum(nj[1:j-1])+1; b1  = sum(nj[1:j])
    y1  = y[a1:b1]
    x1  = matrix(x[a1:b1, ], ncol=p)
    z1  = matrix(z[a1:b1, ], ncol=q1)
    tt1 = ttc[a1:b1]
    cc1 = cc[a1:b1]
    LL1 = LL[a1:b1]
    LU1 = LU[a1:b1]

    muii = x1%*%betas
    Gama = MatDec(tt1, phi1, struc)
    SIGMA = (sigmae*Gama + (z1)%*%D1%*%t(z1))
    SIGMA = (SIGMA + t(SIGMA))/2

    if(sum(cc1)==0){ # No censored observations
      ver[j,] = dmvnorm(x=as.vector(y1), mean=as.vector(muii), sigma=SIGMA, log=TRUE)

    } else {
      if(sum(cc1)>=1){
        if(sum(cc1)==nj[j]){ # All censored observations
          Sc = (SIGMA + t(SIGMA))/2
          ver[j,] = log(pmvnorm(lower=LL1, upper=LU1, mean=as.vector(muii), sigma=Sc)[1])

        } else { # At least one censored observation
          isigma1 = solve(SIGMA[cc1==0,cc1==0])
          muiic = x1[cc1==1,]%*%betas + SIGMA[cc1==1,cc1==0]%*%isigma1%*%(y1[cc1==0]-x1[cc1==0,]%*%betas)
          Sc = SIGMA[cc1==1,cc1==1] - SIGMA[cc1==1,cc1==0]%*%isigma1%*%SIGMA[cc1==0,cc1==1]
          Sc = (Sc + t(Sc))/2
          LL1c = LL1[cc1==1]
          LU1c = LU1[cc1==1]
          ver[j,] = dmvnorm(x=as.vector(y1[cc1==0]), mean=as.vector(muii[cc1==0]), sigma=as.matrix(SIGMA[cc1==0,cc1==0]), log=TRUE) +
            log(as.numeric(pmvnorm(lower=as.vector(LL1c), upper=as.vector(LU1c), mean=as.vector(muiic), sigma=as.matrix(Sc))[1]))
        }
      }
    }
  } # End for
  logvero = sum(ver)
  return(logvero)
}

# Parameter estimation
# ------------------------------------------
censEM.norm = function(y, x, z, cc, lcl, ucl, ind, ttc, data, beta1, sigmae, D1, phi,
                       struc, pAR, precision, informa, itermax, showiter, showerroriter){
  t1 = Sys.time()
  nj = tapply(ind, ind, length)[unique(ind)]
  m  = nlevels(ind)
  N  = length(ind)
  p  = ncol(x)
  q1 = ncol(z)
  m2 = m*q1

  #Initial values
  if (struc=="UNC"){
    phi1 = NULL
    teta = c(beta1, sigmae, D1[upper.tri(D1, diag=T)])
  } else {
    if (struc=="ARp"){
      piis = phi
      phi1 = estphit(piis)
    } else {
      phi1 = phi
    }
    teta = c(beta1, sigmae, phi1, D1[upper.tri(D1, diag=T)])
  } # End if
  iD1 = solve(D1)
  iD1 = (iD1 + t(iD1))/2
  qr = length(D1[lower.tri(D1, diag=T)])

  ## EM algorithm
  criterio = 10
  count = 0
  loglik = loglikFuncN(y=y, cc=cc, x=x, z=z, ttc=ttc, nj=nj, LL=lcl, LU=ucl, betas=beta1, sigmae=sigmae, D1=D1, phi1=phi1, struc=struc)
  if (is.nan(loglik)|is.infinite(abs(loglik))) stop("NaN/infinity initial likelihood")
  loglikvec = NULL
  Theta = NULL

  while ((criterio>precision)&(count<itermax)){
    count = count + 1
    #
    soma1 = matrix(0, q1, q1)
    soma2 = 0
    soma3 = matrix(0, p, p)
    soma4 = matrix(0, p, 1)
    if (informa) soma5 = matrix(0, p, p)
    #
    uyi  = matrix(0, N, m)
    uyyi = matrix(0, N, N)
    ubi  = matrix(0, m2, m)
    ubbi = matrix(0, m2, m2)
    uybi = matrix(0, N,m2)
    yest = matrix(0, N, 1)     # yi = xibeta + times + zibi
    yhi  = matrix(0, N, 1)     # yi = E(yi|Ci, Vi)

    ## E-step: Compute expectations
    ## -----------------------------------------
    for (j in 1:m){
      a1  = sum(nj[1:j-1])+1; b1  = sum(nj[1:j])
      y1  = y[a1:b1]
      x1  = matrix(x[a1:b1, ], ncol=p)
      z1  = matrix(z[a1:b1, ], ncol=q1)
      tt1 = ttc[a1:b1]
      cc1 = cc[a1:b1]
      LL1 = lcl[a1:b1]
      LU1 = ucl[a1:b1]

      muii = x1%*%beta1
      Gama = MatDec(tt1, phi1, struc)
      invGama = solve(Gama)
      SIGMA = (sigmae*Gama + (z1)%*%D1%*%t(z1))
      SIGMA = (SIGMA + t(SIGMA))/2
      SIGMAinv = solve(SIGMA)
      Lambda1  = solve(iD1 + (t(z1)%*%invGama%*%z1)*(1/sigmae))
      Lambda1  = (Lambda1 + t(Lambda1))/2

      if (sum(cc1)==0){ # All variables observed
        uy = matrix(y1, nj[j], 1)
        uyy = y1%*%t(y1)
        ub = (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy-muii)
        ubb = Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy-uy%*%t(muii)-muii%*%t(uy)+muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
        ubb = (ubb + t(ubb))/2
        uyb = (uyy-uy%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)

      } else if (sum(cc1)>=1){ # At least one censored observation

        if (sum(cc1)==nj[j]){ # All observations are censored
          Sc  = (SIGMA + t(SIGMA))/2
          aux = meanvarTMD(lower=LL1, upper=LU1, mu = muii, Sigma=Sc, dist="normal")
          uy  = aux$mean
          uyy = aux$EYY
          ub  = (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy-muii)
          ubb = Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy-uy%*%t(muii)-muii%*%t(uy)+muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
          uyb = (uyy-uy%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)

        } else {
          inSigma00 = solve(SIGMA[cc1==0,cc1==0])
          muiic = x1[cc1==1,]%*%beta1 + SIGMA[cc1==1,cc1==0]%*%inSigma00%*%(y1[cc1==0] - x1[cc1==0,]%*%beta1)
          Sc = SIGMA[cc1==1,cc1==1] - SIGMA[cc1==1,cc1==0]%*%inSigma00%*%SIGMA[cc1==0,cc1==1]
          Sc = (Sc + t(Sc))/2
          LL1c = LL1[cc1==1]
          LU1c = LU1[cc1==1]

          aux = meanvarTMD(lower=LL1c, upper=LU1c, mu=as.vector(muiic), Sigma=Sc, dist="normal")
          w1aux = aux$mean
          w2aux = aux$EYY
          uy = matrix(y1, nj[j], 1)
          uy[cc1==1] = w1aux
          uyy = y1%*%t(y1)
          uyy[cc1==0,cc1==1] = y1[cc1==0]%*%t(w1aux)
          uyy[cc1==1,cc1==0] = w1aux%*%t(y1[cc1==0])
          uyy[cc1==1,cc1==1] = w2aux
          uyy = (uyy + t(uyy))/2
          ub  = (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy-muii)
          ubb = Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy-uy%*%t(muii)-muii%*%t(uy)+muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
          ubb = (ubb + t(ubb))/2
          uyb = (uyy-uy%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
        }
      } # end if sum(cc1)>=1

      soma1 = soma1 + ubb
      soma2 = soma2 + (sum(diag(uyy%*%invGama)) - t(uy)%*%invGama%*%muii - t(muii)%*%invGama%*%uy - sum(diag(t(uyb)%*%invGama%*%z1)) - sum(diag(uyb%*%t(z1)%*%invGama))
                        + t(muii)%*%invGama%*%z1%*%ub + t(ub)%*%t(z1)%*%invGama%*%muii + t(muii)%*%invGama%*%muii + sum(diag(ubb%*%t(z1)%*%invGama%*%z1)))
      soma3 = soma3 + (t(x1)%*%invGama%*%x1)
      soma4 = soma4 + (t(x1)%*%invGama%*%(uy-z1%*%ub))
      if (informa){ soma5 = soma5 + (t(x1)%*%SIGMAinv%*%x1 - t(x1)%*%SIGMAinv%*%(uyy-uy%*%t(uy))%*%SIGMAinv%*%x1) }

      c1 = ((j-1)*q1)+1; d1 = j*q1
      uyi[a1:b1, j] = uy
      uyyi[a1:b1, a1:b1] = uyy
      ubi[c1:d1, j] = ub
      ubbi[c1:d1, c1:d1] = ubb
      uybi[a1:b1, c1:d1] = uyb
      yest[a1:b1] = z1%*%ub + muii
      yhi[a1:b1]  = uy
    } # end for

    yhatorg = apply(yhi, 1, sum) # y original + imputado
    yfit = apply(yest, 1, sum)   # y fitted
    yfit[cc==1] = yhatorg[cc==1]
    bi = matrix(apply(ubi, 1, sum), nrow=m, ncol=q1, byrow=TRUE)

    ## M-step: Uptade parameters
    ## -----------------------------------------
    beta1  = solve(soma3)%*%soma4
    sigmae = (1/N)*(soma2)
    sigmae = as.numeric(sigmae)
    D1 = (1/m)*(soma1)
    iD1 = solve(D1)

    # Update phi1
    if (struc=="UNC"){
      teta1 = c(beta1, sigmae, D1[upper.tri(D1, diag = T)])
    } else {

      if (struc=="DEC"){
        phi1 = optim(phi1, FCi, lower=c(0.01,0.01), upper=c(0.9,0.9), method="L-BFGS-B", hessian=TRUE, beta1=beta1, sigmae=sigmae,
                     ttc=ttc, ubi=ubi, ubbi=ubbi, uybi=uybi, uyyi=uyyi, uyi=uyi, x=x, z=z, nj=nj)$par

      } else if (struc=="ARp"){
        piis = optim(piis, FCiArp, lower=rep(-.999,pAR), upper=rep(.999,pAR), method="L-BFGS-B", hessian=TRUE, beta1=beta1, sigmae=sigmae,
                     ttc=ttc, ubi=ubi, ubbi=ubbi, uybi=uybi, uyyi=uyyi, uyi=uyi, x=x, z=z, nj=nj)$par
        phi1 = estphit(piis)

      } else if(struc=="CS"){
        phi1 = optimize(f=FCiphi1, lower=0.001, upper=0.99, beta1=beta1, sigmae=sigmae,
                        ttc=ttc, ubi=ubi, ubbi=ubbi, uybi=uybi, uyyi=uyyi, uyi=uyi, x=x, z=z, nj=nj, struc=struc)$minimum

      } else if (struc=="MA1"){
        phi1 = optimize(f=FCiphi1, lower=-0.49, upper=0.49, beta1=beta1, sigmae=sigmae,
                        ttc=ttc, ubi=ubi, ubbi=ubbi, uybi=uybi, uyyi=uyyi, uyi=uyi, x=x, z=z, nj=nj, struc=struc)$minimum
      }
      teta1 = c(beta1, sigmae, phi1, D1[upper.tri(D1,diag=T)])
    } # End estimation phi1

    loglik1 = loglikFuncN(y=y, cc=cc, x=x, z=z, ttc=ttc, nj=nj, LL=lcl, LU=ucl, betas=beta1, sigmae=sigmae, D1=D1, phi1=phi1, struc=struc)
    criterio = sqrt(((loglik1/loglik)-1)%*%((loglik1/loglik)-1))
    loglik = loglik1
    loglikvec = c(loglikvec, loglik)
    Theta = rbind(Theta, teta1)
    #
    if (showiter&&!showerroriter) cat("Iteration ",count,"\r")
    if (showerroriter) cat("Iteration ",count," of ",itermax," - criterium =",criterio," - loglik =",loglik,"\r")
  } # End while

  t2 = Sys.time()
  timediff = t2 - t1
  #
  dd = matrix.sqrt(D1)[upper.tri(D1, diag=T)]

  # The information matrix
  if (informa){
    Infbetas = soma5
    Infbetas = (Infbetas + t(Infbetas))/2
  }

  # Selection criteria
  npar = length(c(teta1))
  AICc = -2*loglik + 2*npar
  BICc = -2*loglik + log(N)*npar
  SIC  = -2*loglik + npar*log(N)

  namesD = NULL
  for (k in 1:nrow(D1)) namesD = c(namesD, paste0("D",1:k,k))
  if (struc!="UNC") names(teta1) = c(colnames(x), "sigma2", paste0("phi", 1:length(phi1)), namesD)
  else names(teta1) = c(colnames(x), "sigma2", namesD)

  if (struc!="UNC") estimates = list(beta=beta1, sigma2=sigmae, phi=phi1, dsqrt=dd, D=D1)
  else estimates = list(beta=beta1, sigma2=sigmae, dsqrt=dd, D=D1)

  uhat = rep(1, m);  names(uhat) = names(nj)
  colnames(bi) = colnames(z)

  if (informa) obj.out = list(theta=teta1, iter=count, estimates=estimates, yest=yhatorg, uhat=uhat, loglik.track=loglikvec,
                              random.effects=bi, std.error=sqrt(diag(solve(Infbetas))), loglik=loglik, elapsedTime=timediff,
                              error=criterio, criteria=list(AIC=AICc, BIC=BICc, SIC=SIC))
  else obj.out = list(theta=teta1, iter=count, estimates=estimates, yest=yhatorg, uhat=uhat, loglik.track=loglikvec, random.effects=bi,
                      std.error=NULL, loglik=loglik, elapsedTime=timediff, error=criterio, criteria=list(AIC=AICc, BIC=BICc, SIC=SIC))
  return(obj.out)
}


####################################################################
## Linear Student-t Mixed-Effects Models with Censored Data       ##
####################################################################

# Auxiliary functions
# ------------------------------------------

# Estimate phi1 and phi2 (DEC model)
FCit = function(phiG, beta1, sigmae, ttc, ubi, ubbi, uybi, uyyi, uyi, ui, x, z, nj){
  m = length(nj)[1]
  p = dim(x)[2]
  q1 = dim(z)[2]
  soma = 0

  for (j in 1:m ){
    a1 = sum(nj[1:j-1])+1; b1 = sum(nj[1:j]); c1 = ((j-1)*q1)+1; d1 = j*q1
    x1   = matrix(x[a1:b1, ], ncol=p)
    z1   = matrix(z[a1:b1, ], ncol=q1)
    muii = x1%*%beta1
    tt1  = ttc[a1:b1]
    #
    ub  = matrix(ubi[c1:d1, j], nrow=q1, ncol=1)
    ubb = as.matrix(ubbi[c1:d1, c1:d1])
    uyb = matrix(uybi[a1:b1, c1:d1], ncol=q1)
    uyy = as.matrix(uyyi[a1:b1, a1:b1])
    uy  = matrix(uyi[a1:b1, j], ncol=1)
    u   = ui[j]
    #
    Cii = MatDec(tt1, phiG, "DEC")
    Cii = (Cii + t(Cii))/2
    if(det(Cii)<=0){A = 1} else {A = det(Cii)}
    invCii = solve(Cii)
    #
    Ai = as.vector(sum(diag(uyy%*%invCii)) - sum(diag(invCii%*%((uyb)%*%t(z1)))) - sum(diag(invCii%*%(z1%*%t(uyb)))) + sum(diag(ubb%*%t(z1)%*%invCii%*%z1))
                     - t(uy)%*%invCii%*%muii - t(muii)%*%invCii%*%uy  + t(muii)%*%invCii%*%z1%*%ub + t(ub)%*%t(z1)%*%invCii%*%muii + u*t(muii)%*%invCii%*%muii)
    #
    soma = soma - 0.5*log(A) - (0.5/sigmae)*Ai
  }
  return(-soma)
}

# Estimate phi (CS and MA(1) model)
FCiphi1t = function(phi1, beta1, sigmae, ttc, ubi, ubbi, uybi, uyyi, uyi, ui, x, z, nj, struc){
  m  = length(nj)[1]
  p  = dim(x)[2]
  q1 = dim(z)[2]
  soma = 0

  for (j in 1:m ){
    a1 = sum(nj[1:j-1])+1; b1 = sum(nj[1:j]); c1 = ((j-1)*q1)+1; d1 = j*q1
    x1   = matrix(x[a1:b1, ], ncol=p)
    z1   = matrix(z[a1:b1, ], ncol=q1)
    muii = x1%*%beta1
    tt1  = ttc[a1:b1]
    #
    ub  = matrix(ubi[c1:d1, j], nrow=q1, ncol=1)
    ubb = as.matrix(ubbi[c1:d1, c1:d1])
    uyb = matrix(uybi[a1:b1, c1:d1], ncol=q1)
    uyy = as.matrix(uyyi[a1:b1, a1:b1])
    uy  = matrix(uyi[a1:b1, j], ncol=1)
    u   = ui[j]
    #
    Cii = MatDec(tt1, phi1, struc)
    Cii = (Cii + t(Cii))/2
    if(det(Cii)<=0){A = 1} else {A = det(Cii)}
    invCii = solve(Cii)
    #
    Ai =  as.vector(sum(diag(uyy%*%invCii)) - sum(diag(invCii%*%((uyb)%*%t(z1)))) - sum(diag(invCii%*%(z1%*%t(uyb)))) + sum(diag(ubb%*%t(z1)%*%invCii%*%z1))
                     - t(uy)%*%invCii%*%muii - t(muii)%*%invCii%*%uy  + t(muii)%*%invCii%*%z1%*%ub + t(ub)%*%t(z1)%*%invCii%*%muii + u*t(muii)%*%invCii%*%muii)
    #
    soma = soma - 0.5*log(A) - (0.5/sigmae)*Ai
  }
  return(-soma)
}

# Estimate phi (AR(p) model)
FCiArpt = function(piis, beta1, sigmae, ttc, ubi, ubbi, uybi, uyyi, uyi, ui, x, z, nj){
  m  = length(nj)[1]
  p  = dim(x)[2]
  q1 = dim(z)[2]
  soma = 0

  for (j in 1:m ){
    a1 = sum(nj[1:j-1])+1; b1 = sum(nj[1:j]); c1 = ((j-1)*q1)+1; d1 = j*q1
    x1   = matrix(x[a1:b1, ], ncol=p)
    z1   = matrix(z[a1:b1, ], ncol=q1)
    muii = x1%*%beta1
    tt1  = ttc[a1:b1]
    #
    ub  = matrix(ubi[c1:d1, j], nrow=q1, ncol=1)
    ubb = as.matrix(ubbi[c1:d1, c1:d1])
    uyb = matrix(uybi[a1:b1, c1:d1], ncol=q1)
    uyy = as.matrix(uyyi[a1:b1, a1:b1])
    uy  = matrix(uyi[a1:b1, j], ncol=1)
    u   = ui[j]
    #
    phi2 = estphit(piis)
    Cii = MatDec(tt1, phi2, "ARp")
    Cii = (Cii + t(Cii))/2
    if (det(Cii)<=0){A = 1} else {A = det(Cii)}
    invCii = solve(Cii)
    #
    Ai = as.vector(sum(diag(uyy%*%invCii)) - sum(diag(invCii%*%((uyb)%*%t(z1)))) - sum(diag(invCii%*%(z1%*%t(uyb)))) + sum(diag(ubb%*%t(z1)%*%invCii%*%z1))
                    - t(uy)%*%invCii%*%muii - t(muii)%*%invCii%*%uy  + t(muii)%*%invCii%*%z1%*%ub + t(ub)%*%t(z1)%*%invCii%*%muii + u*t(muii)%*%invCii%*%muii)
    #
    soma = soma - 0.5*log(A) - (0.5/sigmae)*Ai
  }
  return(-soma)
}

# Probability of multivariate Student-t
acumTs = function(y, mu, sigma, nu){
  lower = rep(-Inf, length(y))
  resp = pmvt(lower=lower, upper=y, delta=mu, df=nu, sigma=sigma, type="shifted")[1]
  return(resp)
}


# Log-likelihood function (Student-t)
# ------------------------------------------
loglikFunct = function(nu, y, x, z, cc, ttc, nj, LL, LU, betas, sigmae, D1, phi1, struc){
  m = length(nj)[1]
  p = dim(x)[2]
  q1  = dim(z)[2]
  ver = matrix(0, m, 1)

  for(j in 1:m){
    a1 = sum(nj[1:j-1])+1; b1 = sum(nj[1:j])
    y1  = y[a1:b1]
    x1   = matrix(x[a1:b1, ], ncol=p)
    z1   = matrix(z[a1:b1, ], ncol=q1)
    tt1  = ttc[a1:b1]
    cc1 = cc[a1:b1]
    LL1  = LL[a1:b1]
    LU1  = LU[a1:b1]

    muii = as.vector(x1%*%betas)
    Gama = MatDec(tt1, phi1, struc)
    SIGMA = (sigmae*Gama + (z1)%*%D1%*%t(z1))
    SIGMA = (SIGMA + t(SIGMA))/2

    if (sum(cc1)==0){ # No censored observations
      ver[j,] = dmvt(x=y1, delta=muii, sigma=SIGMA, df=nu, log=TRUE, type="shifted") #dTs(y1,as.vector(muii),SIGMA,nu)

    } else {
      if (sum(cc1)>=1){
        if (sum(cc1)==nj[j]){ # All censored observations
          ver[j,] = log(pmvt(lower=LL1, upper=LU1, delta=muii, df=nu, sigma=SIGMA, type="shifted")[1]) #acumTs(as.vector(y1-muii),rep(0,sum(cc1)),SIGMA,nu)

        } else { # At least one censored observation
          if (sum(cc1)==1){
            isigma00 = solve(SIGMA[cc1==0,cc1==0])
            muiic = muii[cc1==1] + SIGMA[cc1==1,cc1==0]%*%isigma00%*%(y1[cc1==0]-muii[cc1==0])
            muiic = as.numeric(muiic)
            Si = SIGMA[cc1==1,cc1==1] - SIGMA[cc1==1,cc1==0]%*%isigma00%*%SIGMA[cc1==0,cc1==1]
            Si = as.numeric(Si)
            coef1 = as.numeric((mahalanobis(y1[cc1==0], muii[cc1==0], SIGMA[cc1==0,cc1==0]) + nu)/(nu + length(cc1[cc1==0])))
            Sic = sqrt(coef1*Si)
            #
            ver[j,] = dmvt(x=y1[cc1==0], delta=muii[cc1==0], sigma=SIGMA[cc1==0,cc1==0], df=nu, log=TRUE, type="shifted") +
              log(pt(q=(LU1[cc1==1]-muiic)/Sic, df=(nu+length(cc1==0))) - pt(q=(LL1[cc1==1]-muiic)/Sic, df=(nu+length(cc1==0))))

          } else {
            isigma00 = solve(SIGMA[cc1==0,cc1==0])
            muiic = muii[cc1==1] + SIGMA[cc1==1,cc1==0]%*%isigma00%*%(y1[cc1==0]-muii[cc1==0])
            muiic = as.vector(muiic)
            Si = SIGMA[cc1==1,cc1==1] - SIGMA[cc1==1,cc1==0]%*%isigma00%*%SIGMA[cc1==0,cc1==1]
            Si = (Si + t(Si))/2
            coef1 = as.numeric((mahalanobis(y1[cc1==0], muii[cc1==0], SIGMA[cc1==0,cc1==0]) + nu)/(nu + length(cc1[cc1==0])))
            Sic = coef1*Si
            Sic = (Sic + t(Sic))/2
            #
            ver[j,] = dmvt(x=y1[cc1==0], delta=muii[cc1==0], sigma=SIGMA[cc1==0,cc1==0], df=nu, log=TRUE, type="shifted") +
              log(pmvt(lower=LL1[cc1==1], upper=LU1[cc1==1], delta=muiic, df=(nu+length(cc1==0)), sigma=Sic, type="shifted")[1])
          }
        }
      } # End if (sum(cc1)>=1)
    } # End else
  } #end for
  logvero = sum(ver)
  return(logvero)
}


# Parameter estimation
# ------------------------------------------
censEM.t = function(y, x, z, cc, lcl, ucl, ind, ttc, data, beta1, sigmae, D1, phi, nu,
                    nu.fixed, struc, pAR, precision, informa, itermax, showiter, showerroriter){
  t1 = Sys.time()
  nj = tapply(ind, ind, length)[unique(ind)]
  m  = length(nj)
  N  = length(ind)
  p  = ncol(x)
  q1 = ncol(z)
  m2 = m*q1

  #  Initial values
  if (struc=="UNC"){
    phi1 = NULL
    teta = c(beta1, sigmae, D1[upper.tri(D1, diag=T)], nu)
  } else {
    if (struc=="ARp"){
      piis = phi
      phi1 = estphit(piis)
    } else {
      phi1 = phi
    }
    teta = c(beta1, sigmae, phi1, D1[upper.tri(D1, diag=T)], nu)
  } # End if
  iD1 = solve(D1)
  iD1 = (iD1 + t(iD1))/2
  qr = length(D1[lower.tri(D1, diag=T)])

  ## EM algorithm
  criterio = 10
  count = 0
  loglik = loglikFunct(nu=nu, y=y, x=x, z=z, cc=cc, ttc=ttc, nj=nj, LL=lcl, LU=ucl, betas=beta1, sigmae=sigmae, D1=D1, phi1=phi1, struc=struc)
  if (is.nan(loglik)|is.infinite(abs(loglik))) stop("NaN/infinity initial likelihood")
  loglikvec = NULL

  while ((criterio>precision)&(count<itermax)){
    count = count + 1
    #
    soma1 = matrix(0, q1, q1)
    soma2 = 0
    soma3 = matrix(0, p, p)
    soma4 = matrix(0, p, 1)
    if (informa) soma5 = matrix(0, p, p)
    #
    ui   = rep(0, m)
    uyi  = matrix(0, N, m)
    uyyi = matrix(0, N, N)
    ubi  = matrix(0, m2, m)
    ubbi = matrix(0, m2, m2)
    uybi = matrix(0, N,m2)
    yest = matrix(0, N, 1)     # yi = xibeta + zibi
    yhi  = matrix(0, N, 1)     # yi = E(yi|Ci, Vi)
    biest = matrix(0, m2, m)

    ## E-step: Compute expectations
    ## -----------------------------------------
    for (j in 1:m){
      a1  = sum(nj[1:j-1])+1; b1  = sum(nj[1:j])
      y1  = y[a1:b1]
      x1  = matrix(x[a1:b1, ], ncol=p)
      z1  = matrix(z[a1:b1, ], ncol=q1)
      tt1 = ttc[a1:b1]
      cc1 = cc[a1:b1]
      LL1 = lcl[a1:b1]
      LU1 = ucl[a1:b1]

      muii = x1%*%beta1
      Gama = MatDec(tt1, phi1, struc)
      invGama = solve(Gama)
      SIGMA = (sigmae*Gama + (z1)%*%D1%*%t(z1))
      SIGMA = (SIGMA + t(SIGMA))/2
      SIGMAinv = solve(SIGMA)
      Lambda1  = solve(iD1 + (t(z1)%*%invGama%*%z1)*(1/sigmae))
      Lambda1  = (Lambda1 + t(Lambda1))/2

      dm  = as.numeric(t(y1 - muii)%*%SIGMAinv%*%(y1-muii))
      cdm = as.numeric((nu+nj[j])/(nu+dm))

      if (sum(cc1)==0){ # All variables observed
        u = cdm
        uy  = matrix(y1,nj[j],1)*cdm
        uyy = (y1%*%t(y1))*cdm
        ub  = (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy - u*muii)
        ubb = Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy - uy%*%t(muii) - muii%*%t(uy) + u*muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
        uyb = (uyy - uy%*%t(muii))%*%(invGama%*%(z1*(1/sigmae))%*%Lambda1)
        yh = matrix(y1,nj[j],1)
        best = (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(yh - muii)

        if (informa){# Information matrix
          Eu2yy = (cdm^2)*(y1%*%t(y1))
          Eu2y  = (cdm^2)*(matrix(y1,nj[j],1))
          Eu2   = cdm^2
          E2 = Eu2yy - Eu2y%*%t(muii) - muii%*%t(Eu2y) + Eu2*(muii)%*%t(muii)
          E1 = (uy - u*(muii))%*%t(uy - u*(muii))
        }
      } else if (sum(cc1)>=1){ # At least one censored observation

        if(sum(cc1)==nj[j]){
          aux1U = acumTs(as.vector(y1-muii), rep(0,sum(cc1)), as.matrix((nu/(nu+2))*SIGMA), (nu + 2))
          aux2U = acumTs(as.vector(y1-muii), rep(0,sum(cc1)), as.matrix(SIGMA), nu)
          u = as.numeric(aux1U/aux2U)
          auxy = mvtelliptical(as.vector(muii), as.matrix((nu/(nu+2))*SIGMA), lower=as.vector(LL1), upper=as.vector(LU1), "t", nu=(nu+2), n=1e5)

          uy = u*auxy$EY
          uyy = u*auxy$EYY
          ub  = (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy - u*muii)
          ubb = Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy - uy%*%t(muii) - muii%*%t(uy) + u*muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
          uyb = (uyy - uy%*%t(muii))%*%(invGama%*%(z1*(1/sigmae))%*%Lambda1)
          auxb = mvtelliptical(as.vector(muii), as.matrix(SIGMA), lower=as.vector(LL1), upper=as.vector(LU1), "t", nu=(nu), n=1e5)
          yh   = auxb$EY
          best = (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(yh - muii)

          if (informa){## Information matrix
            cp = (((nu+nj[j])/nu)^2)*((gamma((nu+nj[j])/2)*gamma((nu+4)/2))/(gamma(nu/2)*gamma((nu+nj[j]+4)/2)))
            auxw = mvtelliptical(as.vector(muii), as.matrix((nu/(nu + 4))*SIGMA), lower=as.vector(LL1), upper=as.vector(LU1), "t", nu=(nu+4), n=1e5)
            auxEU = acumTs(as.vector(y1-muii), rep(0,sum(cc1)), as.matrix((nu/(nu+4))*SIGMA), (nu+4))
            Eu2yy = cp*(auxEU/aux2U)*auxw$EYY
            Eu2y  = cp*(auxEU/aux2U)*auxw$EY
            Eu2 = cp*(auxEU/aux2U)
            E2  = Eu2yy - Eu2y%*%t(muii) - muii%*%t(Eu2y) + Eu2*(muii)%*%t(muii)
            E1  = (uy - u*(muii))%*%t(uy - u*(muii))
          }

        } else {
          isigma00 = solve(SIGMA[cc1==0,cc1==0])
          n2 = length(cc1[cc1==0])
          muiic = x1[cc1==1,]%*%beta1 + SIGMA[cc1==1,cc1==0]%*%isigma00%*%(y1[cc1==0]-x1[cc1==0,]%*%beta1)
          Si = SIGMA[cc1==1,cc1==1] - SIGMA[cc1==1,cc1==0]%*%isigma00%*%SIGMA[cc1==0,cc1==1]
          Si = (Si+t(Si))/2
          Qy0 = as.numeric(t(y1[cc1==0]-x1[cc1==0,]%*%beta1)%*%isigma00%*%(y1[cc1==0]-x1[cc1==0,]%*%beta1))
          auxQy0  = as.numeric((nu + Qy0)/(nu + n2))
          auxQy02 = as.numeric((nu + Qy0)/(nu + 2 + n2))
          Sc0 = auxQy0*Si
          Sc0til = auxQy02*Si
          LL1c = LL1[cc1==1]
          LU1c = LU1[cc1==1]

          aux1U = acumTs(as.vector(y1[cc1==1]-muiic),rep(0,sum(cc1)),Sc0til,(nu + 2 + n2))
          aux2U = acumTs(as.vector(y1[cc1==1]-muiic),rep(0,sum(cc1)),Sc0,(nu + n2))
          u = as.numeric(aux1U/aux2U)*(1/auxQy0)
          auxy = mvtelliptical(as.vector(muiic), as.matrix(Sc0til), lower=as.vector(LL1c), upper=as.vector(LU1c), "t", nu=(nu + 2 + n2), n=1e5)
          w1aux = auxy$EY
          w2aux = auxy$EYY

          uy = matrix(y1,nj[j],1)*u
          uy[cc1==1] = w1aux*u
          uyy = y1%*%t(y1)*u
          uyy[cc1==0,cc1==1] = u*y1[cc1==0]%*%t(w1aux)
          uyy[cc1==1,cc1==0] = u*w1aux%*%t(y1[cc1==0])
          uyy[cc1==1,cc1==1] = u*w2aux
          uyy = (uyy + t(uyy))/2

          ub = (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uy - u*muii)
          ubb = Lambda1 + (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(uyy - uy%*%t(muii) - muii%*%t(uy) + u*muii%*%t(muii))%*%t(Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)
          uyb = (uyy - uy%*%t(muii))%*%(invGama%*%(z1*(1/sigmae))%*%Lambda1)

          auxb = mvtelliptical(as.vector(muiic), as.matrix(Sc0), lower=as.vector(LL1c), upper=as.vector(LU1c), "t", nu=(nu + n2), n=1e5)$EY
          yh = matrix(y1,nj[j],1)
          yh[cc1==1] = auxb
          best = (Lambda1%*%(t(z1)*(1/sigmae))%*%invGama)%*%(yh - muii)

          if (informa) {# Information matrix
            dp = ((nu+nj[j])^2)*((gamma((nu+nj[j])/2)*gamma((nu+4+n2)/2))/(gamma((nu+n2)/2)*gamma((nu+4+nj[j])/2)))
            auxEU = acumTs(as.vector(y1[cc1==1]-muiic), rep(0,sum(cc1)), as.matrix(as.numeric((nu+Qy0)/(nu+4+n2))*Si), (nu+4+n2))
            auxEw = mvtelliptical(as.vector(muiic), as.matrix(as.numeric((nu+Qy0)/(nu+4+n2))*Si), lower=as.vector(LL1c), upper=as.vector(LU1c), "t", nu=(nu+4+n2), n=1e5)
            Ew1aux = auxEw$EY
            Ew2aux = auxEw$EYY

            Eu2yy = (dp/((nu + Qy0)^2))*(auxEU/aux2U)*y1%*%t(y1)
            Eu2yy[cc1==0,cc1==1] = (dp/((nu + Qy0)^2))*(auxEU/aux2U)*y1[cc1==0]%*%t(Ew1aux)
            Eu2yy[cc1==1,cc1==0] = (dp/((nu + Qy0)^2))*(auxEU/aux2U)* Ew1aux%*%t(y1[cc1==0])
            Eu2yy[cc1==1,cc1==1] = (dp/((nu + Qy0)^2))*(auxEU/aux2U)*Ew2aux

            Eu2y = (dp/((nu + Qy0)^2))*(auxEU/aux2U)*matrix(y1,nj[j],1)
            Eu2y[cc1==1] = (dp/((nu + Qy0)^2))*(auxEU/aux2U)*Ew1aux
            Eu2 = (dp/((nu + Qy0)^2))*(auxEU/aux2U)
            E2 = Eu2yy - Eu2y%*%t(muii) - muii%*%t(Eu2y) + Eu2*(muii)%*%t(muii)
            E1 = (uy - u*(muii))%*%t(uy - u*(muii))
          }
        }
      } # End if

      soma1 = soma1 + ubb
      soma2 = soma2 + (sum(diag(uyy%*%invGama)) - t(uy)%*%invGama%*%muii - t(muii)%*%invGama%*%uy - sum(diag(t(uyb)%*%invGama%*%z1)) - sum(diag(uyb%*%t(z1)%*%invGama))
                       + t(muii)%*%invGama%*%z1%*%ub  + t(ub)%*%t(z1)%*%invGama%*%muii + u*t(muii)%*%invGama%*%muii + sum(diag(ubb%*%t(z1)%*%invGama%*%z1)))
      soma3 = soma3 + (u*t(x1)%*%invGama%*%x1)
      soma4 = soma4 + (t(x1)%*%invGama%*%(uy - z1%*%ub))
      if (informa) soma5 = soma5 + (((nu+nj[j])/(nu+nj[j]+2))*t(x1)%*%SIGMAinv%*%x1 - ((nu+nj[j]+2)/(nu+nj[j]))*t(x1)%*%SIGMAinv%*%(E2)%*%SIGMAinv%*%x1 + (t(x1)%*%SIGMAinv%*%(E1)%*%SIGMAinv%*%x1))

      c1 = ((j-1)*q1)+1; d1 = j*q1
      ui[j] = u
      uyi[a1:b1, j] = uy
      uyyi[a1:b1, a1:b1] = uyy
      ubi[c1:d1, j] = ub
      ubbi[c1:d1, c1:d1] = ubb
      uybi[a1:b1, c1:d1] = uyb
      yest[a1:b1] = z1%*%best + muii
      biest[c1:d1, j] = best
      yhi[a1:b1] = yh
    } # end for

    yhatorg = apply(yhi, 1, sum) # y original + imputado
    yfit = apply(yest, 1, sum)   # y fitted
    yfit[cc==1] = yhatorg[cc==1]
    bi = matrix(apply(biest, 1, sum), nrow=m, ncol=q1, byrow=TRUE)

    ## M-step: Uptade parameters
    ## -----------------------------------------
    beta1 = solve(soma3)%*%soma4
    sigmae = (1/N)*(soma2)
    sigmae = as.numeric(sigmae)
    D1  = (1/m)*(soma1)
    iD1 = solve(D1)

    # Update nu
    if (nu.fixed==FALSE){
      nu = optimize(f=loglikFunct, interval=c(2.1,50), tol=0.001, maximum=TRUE, y=y, x=x, z=z, cc=cc, ttc=ttc, nj=nj,
                    LL=lcl, LU=ucl, betas=beta1, sigmae=sigmae, D1=D1, phi1=phi1, struc=struc)$maximum
    }

    # Update phi1
    if (struc=="UNC"){
      if (nu.fixed==FALSE) teta1 = c(beta1, sigmae, D1[upper.tri(D1, diag = T)], nu)
      else teta1 = c(beta1, sigmae, D1[upper.tri(D1, diag = T)])

    } else {
      if (struc=="DEC"){
        phi1 = optim(phi1, FCit, lower=c(0.01,0.01), upper=c(0.9,0.9), method="L-BFGS-B", hessian=TRUE, beta1=beta1, sigmae=sigmae,
                     ttc=ttc, ubi=ubi, ubbi=ubbi, uybi=uybi, uyyi=uyyi, uyi=uyi, ui=ui, x=x, z=z, nj=nj)$par

      } else if (struc=="ARp"){
        piis = optim(piis, FCiArpt, lower=rep(-.999,pAR), upper=rep(.999,pAR), method="L-BFGS-B", hessian=TRUE, beta1=beta1, sigmae=sigmae,
                     ttc=ttc, ubi=ubi, ubbi=ubbi, uybi=uybi, uyyi=uyyi, uyi=uyi, ui=ui, x=x, z=z, nj=nj)$par
        phi1 = estphit(piis)

      } else if(struc=="CS"){
        phi1 = optimize(f=FCiphi1t, lower=0.001, upper=0.99, tol=0.001, beta1=beta1, sigmae=sigmae,
                        ttc=ttc, ubi=ubi, ubbi=ubbi, uybi=uybi, uyyi=uyyi, uyi=uyi, ui=ui, x=x, z=z, nj=nj, struc=struc)$minimum

      } else if (struc=="MA1"){
        phi1 = optimize(f=FCiphi1t, lower=-0.49, upper=0.49, tol=0.001, beta1=beta1, sigmae=sigmae,
                        ttc=ttc, ubi=ubi, ubbi=ubbi, uybi=uybi, uyyi=uyyi, uyi=uyi, ui=ui, x=x, z=z, nj=nj, struc=struc)$minimum
      }
      if (nu.fixed==FALSE) teta1 = c(beta1, sigmae, phi1, D1[upper.tri(D1,diag=T)], nu)
      else teta1 = c(beta1, sigmae, phi1, D1[upper.tri(D1,diag=T)])
    } # End estimation phi1

    loglik1 = loglikFunct(nu=nu, y=y, x=x, z=z, cc=cc, ttc=ttc, nj=nj, LL=lcl, LU=ucl, betas=beta1, sigmae=sigmae, D1=D1, phi1=phi1, struc=struc)
    criterio = sqrt(((loglik1/loglik)-1)%*%((loglik1/loglik)-1))
    loglik = loglik1
    loglikvec = c(loglikvec, loglik)
    #
    if (showiter&&!showerroriter) cat("Iteration ",count,"\r")
    if (showerroriter) cat("Iteration ",count," of ",itermax," - criterium =",criterio," - loglik =",loglik,"\r")
  } # End while

  t2 = Sys.time()
  timediff = t2 - t1
  #
  dd = matrix.sqrt(D1)[upper.tri(D1, diag=T)]

  # The information matrix
  if (informa){
    Infbetas = soma5
    Infbetas = (Infbetas + t(Infbetas))/2
  }

  # Selection criteria
  npar = length(c(teta1))
  AICc = -2*loglik + 2*npar
  BICc = -2*loglik + log(N)*npar
  SIC  = -2*loglik + npar*log(N)
  #
  namesD = NULL
  for (k in 1:nrow(D1)) namesD = c(namesD, paste0("D",1:k,k))
  if (nu.fixed) {
    if (struc!="UNC") names(teta1) = c(colnames(x), "sigma2", paste0("phi", 1:length(phi1)), namesD)
    else names(teta1) = c(colnames(x), "sigma2", namesD)
  } else {
    if (struc!="UNC") names(teta1) = c(colnames(x), "sigma2", paste0("phi", 1:length(phi1)), namesD, "nu")
    else names(teta1) = c(colnames(x), "sigma2", namesD, "nu")
  }

  if (struc!="UNC") estimates = list(beta=beta1, sigma2=sigmae, phi=phi1, dsqrt=dd, D=D1, nu=nu)
  else estimates = list(beta=beta1, sigma2=sigmae, dsqrt=dd, D=D1, nu=nu)

  uhat = ui
  names(uhat) = names(nj)
  colnames(bi) = colnames(z)

  if (informa) obj.out = list(theta=teta1, iter=count, estimates=estimates, yest=yhatorg, uhat=uhat, loglik.track=loglikvec,
                              random.effects=bi, std.error=sqrt(diag(solve(Infbetas))), loglik=loglik, elapsedTime=timediff,
                              error=criterio, criteria=list(AIC=AICc, BIC=BICc, SIC=SIC))
  else obj.out = list(theta=teta1, iter=count, estimates=estimates, yest=yhatorg, uhat=uhat, loglik.track=loglikvec,
                      random.effects=bi, std.error=NULL, loglik=loglik, elapsedTime=timediff, error=criterio,
                      criteria=list(AIC=AICc, BIC=BICc, SIC=SIC))
  return(obj.out)
}

