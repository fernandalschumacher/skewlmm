library(skewlmm)

nj1 = 5
m = 30
time = rep(1:nj1, times=m)
groups = as.factor(rep(1:m, each=nj1))
dat1 = rsmsn.clmm(time = time, ind = groups, x = cbind(1,time), z = rep(1,m*nj1),
                  sigma2=0.6, D=0.5*diag(1), beta=c(1,2), depStruct="ARp",
                  phi=0.4, pcens = .1, type = "left")

# Estimation using UNC
fm0 = smn.clmm(dat1, formFixed=y~x, groupVar="ind", formRandom = ~x, ci="ci", ucl="ucl")
fm0$estimates$beta
summary(fm0)

ranef(fm0)

coef(fm0)

fm1 = smn.lmm(dat1, formFixed=y~x, groupVar="ind")
fm1$estimates$beta
sfm1<-summary(fm1)
sfm1$tableFixed

library(tidyverse)
lm(y~x, data = dat1) %>% broom::tidy()

confint(fm1)

