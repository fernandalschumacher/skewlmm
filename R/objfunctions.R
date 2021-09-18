#### objfunctions for SMSN
objfn.skewAR <- function(thetav,y, x, z,time,ind,distr,pAR,...) {
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2#;N<-length(y)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  Gammab <- Dmatrix(thetav[(p+2):(p+1+q2)])
  Deltab<-thetav[(p+2+q2):(p+1+q2+q1)]
  piAR<-thetav[(p+2+q2+q1):(p+1+q2+q1+pAR)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+1+q2+q1+pAR))]
  #
  phiAR <- estphit(piAR)
  N <-length(ind)
  #logveroARpi function
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalAR,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtAR,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsAR,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnAR,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiAR=phiAR))
  lv
}
objfn.skewUNC <- function(thetav,y, x, z,ind,distr,...) {
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2#;N<-length(y)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  Gammab <- Dmatrix(thetav[(p+2):(p+1+q2)])
  Deltab<-thetav[(p+2+q2):(p+1+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+1+q2+q1))]
  #
  N <-length(ind)
  #logvero function
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormal,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljt,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljs,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcn,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  lv
}
objfn.skewCS <- function(thetav,y, x, z,ind,distr,...) {
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2#;N<-length(y)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  Gammab <- Dmatrix(thetav[(p+2):(p+1+q2)])
  Deltab<-thetav[(p+2+q2):(p+1+q2+q1)]
  phiCS<-thetav[(p+2+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+2+q2+q1))]
  #
  N <-length(ind)
  #logveroCS function
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalCS,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiCS=phiCS))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtCS,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiCS=phiCS))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsCS,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiCS=phiCS))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnCS,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiCS=phiCS))
  lv
}
objfn.skewDEC <- function(thetav,y, x, z,time,ind,distr,...) {
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2#;N<-length(y)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  Gammab <- Dmatrix(thetav[(p+2):(p+1+q2)])
  Deltab<-thetav[(p+2+q2):(p+1+q2+q1)]
  phiDEC<-thetav[(p+2+q2+q1)]
  thetaDEC<-thetav[(p+3+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+3+q2+q1))]
  #
  N <-length(ind)
  #logveroDEC function
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalDEC,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtDEC,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsDEC,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnDEC,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  lv
}
objfn.skewCAR1 <- function(thetav,y, x, z,time,ind,distr,...) {
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2#;N<-length(y)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  Gammab <- Dmatrix(thetav[(p+2):(p+1+q2)])
  Deltab<-thetav[(p+2+q2):(p+1+q2+q1)]
  phiDEC<-thetav[(p+2+q2+q1)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+2+q2+q1))]
  #
  N <-length(ind)
  #logveroCAR1 function
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalCAR1,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtCAR1,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsCAR1,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnCAR1,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae,phiDEC=phiDEC))
  lv
}

#### objfunctions for SMN
objfn.AR <- function(thetav,y, x, z,time,ind,distr,pAR,...) {
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2#;N<-length(y)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  D1 <- Dmatrix(thetav[(p+2):(p+1+q2)])
  piAR<-thetav[(p+2+q2):(p+1+q2+pAR)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+1+q2+pAR))]
  #
  phiAR <- estphit(piAR)
  N <-length(ind)
  #logveroARpis function
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalARs,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtARs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsARs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiAR=phiAR))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnARs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiAR=phiAR))
  lv
}
objfn.UNC <- function(thetav,y, x, z,ind,distr,...) {
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2#;N<-length(y)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  D1 <- Dmatrix(thetav[(p+2):(p+1+q2)])
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+1+q2))]
  #
  N <-length(ind)
  #logveros function
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormals,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljts,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljss,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcns,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  lv
}
objfn.CS <- function(thetav,y, x, z,ind,distr,...) {
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2#;N<-length(y)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  D1 <- Dmatrix(thetav[(p+2):(p+1+q2)])
  phiCS<-thetav[(p+2+q2)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+2+q2))]
  #
  N <-length(ind)
  #logveroCSs function
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalCSs,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae,phiCS=phiCS))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtCSs,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae,phiCS=phiCS))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsCSs,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae,phiCS=phiCS))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnCSs,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae,phiCS=phiCS))
  lv
}
objfn.DEC <- function(thetav,y, x, z,time,ind,distr,...) {
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2#;N<-length(y)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  D1 <- Dmatrix(thetav[(p+2):(p+1+q2)])
  phiDEC<-thetav[(p+2+q2)]
  thetaDEC<-thetav[(p+3+q2)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+3+q2))]
  #
  N <-length(ind)
  #logveroDECs function
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalDECs,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtDECs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsDECs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnDECs,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC,thetaDEC=thetaDEC))
  lv
}
objfn.CAR1 <- function(thetav,y, x, z,time,ind,distr,...) {
  p<-ncol(x);q1<-ncol(z);q2 <- q1*(q1+1)/2#;N<-length(y)
  beta1<-matrix(thetav[1:p],ncol=1)#thetav[1:p]
  sigmae<-as.numeric(thetav[p+1])#thetav[p+1]
  D1 <- Dmatrix(thetav[(p+2):(p+1+q2)])
  phiDEC<-thetav[(p+2+q2)]
  if (distr=="sn") {
    nu <- NULL
  } else nu<-thetav[-(1:(p+2+q2))]
  #
  N <-length(ind)
  #logveroCAR1s function
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalCAR1s,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtCAR1s,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsCAR1s,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnCAR1s,nu=nu,y=y,x=x,z=z,time=time,beta1=beta1,D1=D1,sigmae=sigmae,phiDEC=phiDEC))
  lv
}