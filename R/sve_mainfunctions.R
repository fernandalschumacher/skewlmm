# Sandwich Variance Estimator --------------------------------------------------
# Main functions

# SVE for all parameters -------------------------------------------------------
gen_der <- function(object, ...){
  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<- object$groupVar
  depStruct <- object$depStruct
  distr <- object$distr
  if (distr=="norm") distr="sn"
  if (distr=="t") distr="st"
  if (distr=="sl"|| distr =="ssl") distr="ss"
  if (distr=="cn") distr="scn"

  data <- object$data
  x <- model.matrix(formFixed,data=data)
  y <- data[,all.vars(formFixed)[1]]
  z <- model.matrix(formRandom,data=data)
  p = ncol(x); q1=ncol(z)
  q2 = q1*(q1+1)/2
  ind <- data[,groupVar]
  nj <- as.vector(table(ind))
  n <- object$n
  N <- length(y)

  if(depStruct == "ARp"){
    beta1 <- object$estimates$beta
    sigmae <- object$estimates$sigma2
    phiAR <- object$estimates$phi
    D1 <- object$estimates$D
    if(class(object)[1] == "SMN") {lambda <- rep(0,nrow(D1))
    } else{lambda <- object$estimates$lambda}
    Dsqrti <- matrix.sqrt(D1)
    pAR <- length(phiAR)
    timeVar <- object$timeVar

    if (is.null(timeVar)) {
      time = flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
    } else time = data[,timeVar]

    if(distr == "sn"){nu = 1} else {nu <- object$estimates$nu}
    # Gerando uma amostra
    dadosi <- tapply(1:N,ind,gerar_smsn,x=x,z=z,sigma2=sigmae,Dsqrti=Dsqrti,
                     beta1=beta1,lambda=lambda,distr=distr,nu=nu,ind=ind,time=time,
                     depStruct=depStruct,phi=phiAR) %>% do.call("rbind",.) #key: bind_rows()
    names(dadosi)[1] <- all.vars(object$formula$formFixed)[1]
    yi <- dadosi[,all.vars(formFixed)[1]]

    # Calculando as derivadas
    if(class(object)[1] == "SMSN"){
      derivates_listi <- tapply(1:N, ind, derivatesARi, y=yi, x=x, z=z, time = time, beta1=beta1,
                                sigmae=sigmae, phiAR = phiAR, D1=D1, lambda=lambda, distr=distr, nu=nu)
    } else{
      derivates_listi <- tapply(1:N, ind, derivatesARis, y=yi, x=x, z=z, time = time, beta1=beta1,
                                sigmae=sigmae, phiAR = phiAR, D1=D1, lambda=lambda, distr=distr, nu=nu)
    }
  } else if (depStruct == "UNC"){
    beta1 <- object$estimates$beta
    sigmae <- object$estimates$sigma2
    D1 <- object$estimates$D
    if(class(object)[1] == "SMN") {lambda <- rep(0,nrow(D1))
    } else{lambda <- object$estimates$lambda}
    Dsqrti <- matrix.sqrt(D1)
    timeVar <- object$timeVar

    if (is.null(timeVar)) {
      time = flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
    } else time = data[,timeVar]

    if(distr == "sn"){nu = 1}else{nu <- object$estimates$nu}
    # Gerando uma amostra
    dadosi <- tapply(1:N,ind,gerar_smsn,x=x,z=z,sigma2=sigmae,Dsqrti=Dsqrti,
                     beta1=beta1,lambda=lambda,distr=distr,nu=nu,ind=ind,time=time,
                     depStruct=depStruct, phi = NULL) %>% do.call("rbind",.) #bind_rows()
    names(dadosi)[1] <- all.vars(object$formula$formFixed)[1]
    yi <- dadosi[,all.vars(formFixed)[1]]

    # Calculando as derivadas
    if(class(object)[1] == "SMSN"){
      derivates_listi <- tapply(1:N, ind, derivatesUNC, y=yi, x=x, z=z, time = time, beta1=beta1,
                                sigmae=sigmae, D1=D1, lambda=lambda, distr=distr, nu=nu)
    } else {
      derivates_listi <- tapply(1:N, ind, derivatesUNCs, y=yi, x=x, z=z, time = time, beta1=beta1,
                                sigmae=sigmae, D1=D1, lambda=lambda, distr=distr, nu=nu)
    }
  }
  if (depStruct == "DEC"){
    beta1 <- object$estimates$beta
    sigmae <- object$estimates$sigma2
    phiDEC <- object$estimates$phi[1]
    thetaDEC <- object$estimates$phi[2]
    D1 <- object$estimates$D
    if(class(object)[1] == "SMN") {lambda <- rep(0,nrow(D1))
    } else{lambda <- object$estimates$lambda}
    Dsqrti <- matrix.sqrt(D1)
    timeVar <- object$timeVar

    if (is.null(timeVar)) {
      time = flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
    } else time = data[,timeVar]

    if(distr == "sn"){nu = 1}else{nu <- object$estimates$nu}
    # Gerando uma amostra
    dadosi <- tapply(1:N,ind,gerar_smsn,x=x,z=z,sigma2=sigmae,Dsqrti=Dsqrti,
                     beta1=beta1,lambda=lambda,distr=distr,nu=nu,ind=ind,time=time,
                     depStruct=depStruct, phi= c(phiDEC, thetaDEC)) %>% do.call("rbind",.) #key: bind_rows()
    names(dadosi)[1] <- all.vars(object$formula$formFixed)[1]
    yi <- dadosi[,all.vars(formFixed)[1]]

    # Calculando as derivadas
    if(class(object)[1] == "SMSN"){
      derivates_listi <- tapply(1:N, ind, derivatesDEC, y=yi, x=x, z=z, time = time, beta1=beta1,
                                sigmae=sigmae, phiDEC = phiDEC, thetaDEC = thetaDEC, D1=D1, lambda=lambda, distr=distr, nu=nu)
    } else {
      derivates_listi <- tapply(1:N, ind, derivatesDECs, y=yi, x=x, z=z, time = time, beta1=beta1,
                                sigmae=sigmae, phiDEC = phiDEC, thetaDEC = thetaDEC, D1=D1, lambda=lambda, distr=distr, nu=nu)
    }
  }

  score_list <- lapply(derivates_listi, `[[`, 1) # Primeiras derivadas
  hessian_list <- lapply(derivates_listi, `[[`, 2) # Segundas derivadas

  #sumscorei <- Reduce("+", score_list)
  #prodscorei <- sumscorei%*%t(sumscorei)
  prodscorei <- Reduce("+", lapply(score_list, function(s) s %*% t(s)))
  hessiani <- Reduce("+", hessian_list)

  return(list(prodscorei = prodscorei, hessiani = hessiani))
}

sandwichvar <- function(object, MCiter = 100,  parallel = TRUE, seed = 123){
  if(!inherits(object,c("SMSN","SMN"))) stop("object must inherit from class SMSN or SMN")
  #if (is.null(ncores)){ncores <- availableCores()} # - 1}

  # Calculando as esperanças usando monte carlo
  #plan(multisession, workers = ncores)
  # mc_list <- future_map(.x = seq_len(MCiter), .f = gen_der, object = object,
  #                       .options = furrr_options(seed = 123))
  #plan(sequential)
  if (parallel) {
    with(plan(multisession, workers = availableCores(omit = 1)), local = TRUE)
    mc_list <- suppressMessages(future_map(.x = seq_len(MCiter), .f = gen_der, object = object,
                                    .options = furrr_options(seed = seed)))
  } else{
    mc_list <- suppressMessages(map(.x = seq_len(MCiter), .f = gen_der, object = object,
                             .options = furrr_options(seed = seed)))
  }

  # Organizando as listas mc
  resultmc <- list_transpose(mc_list)
  auxresult <- list(hessian = resultmc$hessiani, prodscore = resultmc$prodscorei)
  sums <- map(auxresult, function(x) Reduce('+', x))

  # Esperancas
  Atheta <- sums$hessian/MCiter
  Btheta = sums$prodscore/MCiter

  # Diagnostic 1: 1/kappa — W5 near-singularity check
  eigs        <- abs(eigen(Atheta, only.values = TRUE)$values)
  rel_min_eig <- min(eigs) / max(eigs)

  if (rel_min_eig < 1e-5) {
    warning(
      "Atheta is near-singular (1/kappa = ", formatC(rel_min_eig, format = "e", digits = 2), "). ",
      "SVE standard errors may be unreliable. ",
      "Consider checking model convergence or using method = 'asymptotic' or 'bootstrap'."
    )
  }

  sAtheta <- solve(Atheta)

  # Variancia sanduiche
  Ctheta = sAtheta%*%Btheta%*%sAtheta
  # Erro Padrao
  #stderror <- if(object$distr %in% c("st", "ssl", "scn")){
  #  c(sqrt(diag(Ctheta)), NA)
  #}else{c(sqrt(diag(Ctheta)))}

  # Diagnostic 2: SVE/ASE ratio — output-level reliability check
  ase   <- object$std.error[seq_along(diag(Ctheta))]
  ratio <- sqrt(abs(diag(Ctheta))) / ase
  bad   <- names(which(ratio > 10 & diag(Ctheta) > 0))

  if (length(bad) > 0) {
    warning(
      "SVE standard errors are more than 10x the asymptotic SE for: ",
      paste(bad, collapse = ", "), ". ",
      "These estimates are likely unreliable for this dataset."
    )
  }

  diag_C <- diag(Ctheta)
  bad_params <- which(diag_C < 0)
  if (length(bad_params) > 0) {
    warning("Sandwich variance estimator produced indefinite Ctheta for parameter(s): ",
            paste(names(object$theta)[bad_params], collapse = ", "),
            ".\nThis typically occurs with small n or weak skewness signal. ",
            "Standard errors for these parameters are set to NA.")
    diag_C[bad_params] <- NA
  }
  n_nu <- sum(grepl("^nu", names(object$theta)))
  stderror <- c(sqrt(diag_C), rep(NA, n_nu))
  names(stderror) <- names(object$theta)
  return(list(std.error = stderror, a_kappa = 1/rel_min_eig))
}

# SMSN - SVE for fixed effects -------------------------------------------------
gen_derBetas <- function(object, ...){
  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<- object$groupVar
  depStruct <- object$depStruct
  distr <- object$distr

  data <- object$data
  x <- model.matrix(formFixed,data=data)
  y <- data[,all.vars(formFixed)[1]]
  z <- model.matrix(formRandom,data=data)
  p = ncol(x); q1=ncol(z)
  q2 = q1*(q1+1)/2
  ind <- data[,groupVar]
  nj <- as.vector(table(ind))
  n <- object$n
  N <- length(y)

  if(depStruct == "ARp"){
    beta1 <- object$estimates$beta
    sigmae <- object$estimates$sigma2
    phiAR <- object$estimates$phi
    D1 <- object$estimates$D
    if(class(object)[1] == "SMN") {lambda <- rep(0,nrow(D1))
    } else{lambda <- object$estimates$lambda}
    Dsqrti <- matrix.sqrt(D1)
    pAR <- length(phiAR)
    timeVar <- object$timeVar

    if (is.null(timeVar)) {
      time = flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
    } else time = data[,timeVar]

    if(distr == "sn"){ nu = 1}else{nu <- object$estimates$nu}
    # Gerando uma amostra
    dadosi <- tapply(1:N,ind,gerar_smsn,x=x,z=z,sigma2=sigmae,Dsqrti=Dsqrti,
                     beta1=beta1,lambda=lambda,distr=distr,nu=nu,ind=ind,time=time,
                     depStruct=depStruct,phi=phiAR) %>% do.call("rbind",.) #bind_rows()
    names(dadosi)[1] <- all.vars(object$formula$formFixed)[1]
    yi <- dadosi[,all.vars(formFixed)[1]]

    # Calculando as derivadas
    derivates_listi <- tapply(1:N, ind, derBetasARi, y=yi, x=x, z=z, time = time, beta1=beta1,
                              sigmae=sigmae, phiAR = phiAR, D1=D1, lambda=lambda, distr=distr, nu=nu)
  }
  if(depStruct == "UNC"){
    beta1 <- object$estimates$beta
    sigmae <- object$estimates$sigma2
    D1 <- object$estimates$D

    if(class(object)[1] == "SMN") {lambda <- rep(0,nrow(D1))
    } else{lambda <- object$estimates$lambda}
    Dsqrti <- matrix.sqrt(D1)
    timeVar <- object$timeVar

    if (is.null(timeVar)) {
      time = flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
    } else time = data[,timeVar]

    if(distr == "sn"){ nu = 1}else{nu <- object$estimates$nu}
    # Gerando uma amostra
    dadosi <- tapply(1:N,ind,gerar_smsn,x=x,z=z,sigma2=sigmae,Dsqrti=Dsqrti,
                     beta1=beta1,lambda=lambda,distr=distr,nu=nu,ind=ind,time=time,
                     depStruct=depStruct, phi = NULL) %>% do.call("rbind",.) #bind_rows()
    names(dadosi)[1] <- all.vars(object$formula$formFixed)[1]
    yi <- dadosi[,all.vars(formFixed)[1]]

    # Calculando as derivadas
    derivates_listi <- tapply(1:N, ind, derBetasUNC, y=yi, x=x, z=z, time = time, beta1=beta1,
                              sigmae=sigmae, D1=D1, lambda=lambda, distr=distr, nu=nu)
  }
  if(depStruct == "DEC"){
    beta1 <- object$estimates$beta
    sigmae <- object$estimates$sigma2
    phiDEC <- object$estimates$phi[1]
    thetaDEC <- object$estimates$phi[2]
    D1 <- object$estimates$D
    if(class(object)[1] == "SMN") {lambda <- rep(0,nrow(D1))
    } else{lambda <- object$estimates$lambda}
    Dsqrti <- matrix.sqrt(D1)
    # pAR <- length(phiAR)
    timeVar <- object$timeVar

    if (is.null(timeVar)) {
      time = flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
    } else time = data[,timeVar]

    if(distr == "sn"){ nu = 1}else{nu <- object$estimates$nu}
    # Gerando uma amostra
    dadosi <- tapply(1:N,ind,gerar_smsn,x=x,z=z,sigma2=sigmae,Dsqrti=Dsqrti,
                     beta1=beta1,lambda=lambda,distr=distr,nu=nu,ind=ind,time=time,
                     depStruct=depStruct, phi = c(phiDEC,thetaDEC)) %>% do.call("rbind",.) #bind_rows()
    names(dadosi)[1] <- all.vars(object$formula$formFixed)[1]
    yi <- dadosi[,all.vars(formFixed)[1]]

    # Calculando as derivadas
    derivates_listi <- tapply(1:N, ind, derBetasDEC, y=yi, x=x, z=z, time = time, beta1=beta1,
                              sigmae=sigmae, phiDEC = phiDEC, thetaDEC=thetaDEC, D1=D1, lambda=lambda, distr=distr, nu=nu)
  }


  score_list <- lapply(derivates_listi, `[[`, 1) # Primeiras derivadas
  hessian_list <- lapply(derivates_listi, `[[`, 2) # Segundas derivadas

  # sumscorei <- Reduce("+", score_list)
  # prodscorei <- sumscorei%*%t(sumscorei)
  prodscorei <- Reduce("+", lapply(score_list, function(s) s %*% t(s)))
  hessiani <- Reduce("+", hessian_list)

  return(list(prodscorei = prodscorei, hessiani = hessiani))
}

sandwichvarBetas <- function(object, MCiter = 100,  parallel = TRUE, seed = 123){
  if(!inherits(object,c("SMSN","SMN"))) stop("object must inherit from class SMSN or SMN")
  #if (is.null(ncores)){ncores <- availableCores() - 1}
  if (object$distr=="norm") object$distr="sn"
  if (object$distr=="t") object$distr="st"
  if (object$distr=="sl") object$distr="ss"
  if (object$distr=="cn") object$distr="scn"
  # Calculando as esperanças usando monte carlo
  #plan(multisession, workers = ncores)
  # mc_list <- future_map(.x = seq_len(MCiter), .f = gen_derBetas, object = object,
  #                       .options = furrr_options(seed = seed))
  if (parallel) {
    with(plan(multisession, workers = availableCores(omit = 1)), local = TRUE)
    mc_list<-suppressMessages(future_map(.x = seq_len(MCiter), .f = gen_derBetas, object = object,
                                         .options = furrr_options(seed = seed)))
  } else{
    mc_list<-suppressMessages(map(.x = seq_len(MCiter), .f = gen_derBetas, object = object,
                                  .options = furrr_options(seed = seed)))
  }
  #plan(sequential)

  # Organizando as listas mc
  resultmc <- list_transpose(mc_list)
  auxresult <- list(hessian = resultmc$hessiani, prodscore = resultmc$prodscorei)
  sums <- map(auxresult, function(x) Reduce('+', x))

  # Esperancas
  Atheta <- sums$hessian/MCiter
  sAtheta <- solve(Atheta)
  Btheta = sums$prodscore/MCiter
  # Variancia sanduiche
  Ctheta = sAtheta%*%Btheta%*%sAtheta
  # Erro Padrao
  stderror <- sqrt(diag(Ctheta))
  names(stderror) <- names(object$theta[1:(length(object$estimates$beta))])

  return(list(std.error = stderror))
  #return(list(score = resultmc$scorei, hessian = resultmc$hessiani, std.error = stderror))
}

# SMN - SVE for fixed effects using exact expectation --------------------------
sandwichvarBetasExpec <- function(object){
  if(!inherits(object,c("SMN"))) stop("object must inherit from class SMN")

  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<- object$groupVar
  depStruct <- object$depStruct
  distr <- object$distr
  if (distr=="norm") distr="sn"
  if (distr=="t") distr="st"
  if (distr=="sl") distr="ss"
  if (distr=="cn") distr="scn"

  data <- object$data
  x <- model.matrix(formFixed,data=data)
  y <- data[,all.vars(formFixed)[1]]
  z <- model.matrix(formRandom,data=data)
  p = ncol(x); q1=ncol(z)
  q2 = q1*(q1+1)/2
  ind <- data[,groupVar]
  nj <- as.vector(table(ind))
  n <- object$n
  N <- length(y)

  if(depStruct == "ARp"){
    beta1 <- object$estimates$beta
    sigmae <- object$estimates$sigma2
    phiAR <- object$estimates$phi
    D1 <- object$estimates$D
    lambda <- rep(0,nrow(D1))
    Dsqrti <- matrix.sqrt(D1)
    pAR <- length(phiAR)
    timeVar <- object$timeVar
    if(distr == "sn"){nu = 1}else{nu <- object$estimates$nu}

    if (is.null(timeVar)) {
      time = flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
    } else time = data[,timeVar]

    # Calculando as esperancas
    expects_listi <- tapply(1:N, ind, expectBetasARi, y=y, x=x, z=z, time = time, beta1=beta1,
                            sigmae=sigmae, phiAR=phiAR, D1=D1,lambda =lambda, distr=distr, nu=nu)
  } else if (depStruct == "UNC"){
    beta1 <- object$estimates$beta
    sigmae <- object$estimates$sigma2
    D1 <- object$estimates$D
    lambda <- rep(0,nrow(D1))
    Dsqrti <- matrix.sqrt(D1)
    timeVar <- object$timeVar
    if(distr == "sn"){nu = 1}else{nu <- object$estimates$nu}

    if (is.null(timeVar)) {
      time = flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
    } else time = data[,timeVar]

    # Calculando as esperancas
    expects_listi <- tapply(1:N, ind, expectBetasUNC, y=y, x=x, z=z, time = time, beta1=beta1,
                            sigmae=sigmae, D1=D1,lambda =lambda, distr=distr, nu=nu)
  } else if(depStruct == "DEC"){
    beta1 <- object$estimates$beta
    sigmae <- object$estimates$sigma2
    phiDEC <- object$estimates$phi[1]
    thetaDEC <- object$estimates$phi[2]
    D1 <- object$estimates$D
    lambda <- rep(0,nrow(D1))
    Dsqrti <- matrix.sqrt(D1)
    #pAR <- length(phiAR)
    timeVar <- object$timeVar
    if(distr == "sn"){nu = 1}else{nu <- object$estimates$nu}

    if (is.null(timeVar)) {
      time = flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
    } else time = data[,timeVar]

    # Calculando as esperancas
    expects_listi <- tapply(1:N, ind, expectBetasDEC, y=y, x=x, z=z, time = time, beta1=beta1,
                            sigmae=sigmae, phiDEC=phiDEC, thetaDEC=thetaDEC, D1=D1,
                            lambda =lambda, distr=distr, nu=nu)
  }
  E2_list <- lapply(expects_listi, `[[`, 1)
  E3_list <- lapply(expects_listi, `[[`, 2)

  Btheta <- Reduce("+", E3_list)
  Atheta <- Reduce("+", E2_list) - Btheta
  sAtheta <- solve(Atheta)

  # Variancia sanduiche
  Ctheta = sAtheta%*%Btheta%*%sAtheta
  # Erro Padrao
  stderror <- sqrt(diag(Ctheta))
  names(stderror) <- names(object$theta[1:(length(object$estimates$beta))])

  return(std.error = stderror)
}

# SMN - Inverse of Fisher informations - for fixed effects ---------------------
IFisherInverse <- function(object){
  if(!inherits(object,c("SMN"))) stop("object must inherit from class SMN")

  formFixed <- object$formula$formFixed
  formRandom <- object$formula$formRandom
  groupVar<- object$groupVar
  depStruct <- object$depStruct
  distr <- object$distr
  if (distr=="norm") distr="sn"
  if (distr=="t") distr="st"
  if (distr=="sl") distr="ss"
  if (distr=="cn") distr="scn"

  data <- object$data
  x <- model.matrix(formFixed,data=data)
  y <- data[,all.vars(formFixed)[1]]
  z <- model.matrix(formRandom,data=data)
  p = ncol(x); q1=ncol(z)
  q2 = q1*(q1+1)/2
  ind <- data[,groupVar]
  nj <- as.vector(table(ind))
  n <- object$n
  N <- length(y)

  if(depStruct == "ARp"){
    beta1 <- object$estimates$beta
    sigmae <- object$estimates$sigma2
    phiAR <- object$estimates$phi
    D1 <- object$estimates$D
    lambda <- rep(0,nrow(D1))
    Dsqrti <- matrix.sqrt(D1)
    pAR <- length(phiAR)
    timeVar <- object$timeVar
    if(distr == "sn"){nu = 1}else{nu <- object$estimates$nu}

    if (is.null(timeVar)) {
      time = flatten_int(tapply(ind,ind,function(x.) seq_along(x.)))
    } else time = data[,timeVar]

    # Calculando as esperancas
    expects_listi <- tapply(1:N, ind, expectBetasARi, y=y, x=x, z=z, time = time, beta1=beta1,
                            sigmae=sigmae, phiAR=phiAR, D1=D1,lambda =lambda, distr=distr, nu=nu)
  }
  E2_list <- lapply(expects_listi, `[[`, 1)
  E3_list <- lapply(expects_listi, `[[`, 2)

  infisher <- -Reduce("+", E2_list) + Reduce("+", E3_list)
  sinfisher <- solve(infisher) # Variancia aproximada
  # Erro Padrao
  stderror <- sqrt(diag(infisher))
  names(stderror) <- names(object$theta[1:(length(object$estimates$beta))])

  return(std.error = stderror)
}
