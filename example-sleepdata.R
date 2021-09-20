library(tidyverse)
library(gridExtra)
library(skewlmm)
library(scales)
library(knitr)

# loading the dataset
data("sleepstudy", package = "lme4")
sleepstudy <- subset(sleepstudy,Days>=2)
sleepstudy$Days <- sleepstudy$Days-2
sleepstudy$Dayst <- sleepstudy$Days-3.5

ggplot(sleepstudy,aes(x=Days,y=Reaction,group=Subject)) + geom_line(alpha=.4) +
  stat_summary(aes(group = 1), geom = "line", fun= mean, colour=1, size=1) +
  scale_x_continuous(breaks = pretty_breaks()) + ylab("reaction time") + xlab("days") +
  theme_minimal()

# fitting a N-LMM using package nlme
fitlme <- nlme::lme(Reaction~Dayst, data = sleepstudy, random = ~Dayst|Subject)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m, ymin=ymin, ymax=ymax))
}
g1 <- nlme::ranef(fitlme) %>% dplyr::rename(`intercepts`=`(Intercept)`,
                                          `slopes`=Dayst) %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(x = value)) + geom_histogram(bins = 7,aes(y=..density..)) +
  theme_minimal() + facet_wrap(~ name, scales = "free") + xlab('') + ylab('density')
g2 <- nlme::ranef(fitlme) %>% dplyr::rename(`intercepts`=`(Intercept)`,
                                          `slopes`=Dayst) %>%
  pivot_longer(cols = everything()) %>%
  ggplot(aes(sample = value)) + geom_qq() + geom_qq_line()+
  facet_wrap(~name,scales = 'free') + theme_minimal()
gridExtra::grid.arrange(g1, g2, ncol = 1)


# fitting other models using the skewlmm package
fit_norm <- smn.lmm(data = sleepstudy, formFixed = Reaction ~ Dayst,
                    formRandom = ~Dayst, groupVar = "Subject")

fit_sl <- update(object = fit_norm, distr = "sl")

fit_ssl <- smsn.lmm(data = sleepstudy, formFixed = Reaction ~ Dayst,
                    formRandom = ~Dayst, groupVar = "Subject",
                    distr = "ssl")
fit_ssl

# testing if lambda = 0
lr.test(fit_sl, fit_ssl)

# assessing the goodness of fit using a Healy-type plot
grid.arrange(healy.plot(fit_norm, calcCI = TRUE),
             healy.plot(fit_sl, calcCI = TRUE), nrow = 1)


# considering different dependence structures
fit_sl_ar1 <- update(fit_sl, depStruct = "ARp", pAR = 1)
fit_sl_ar2 <- update(fit_sl, depStruct = "ARp", pAR = 2)
fit_sl_DEC <- update(fit_sl, depStruct = "DEC", timeVar = "Days",
                     control = lmmControl(parallelphi = TRUE))
fit_sl_DEC_seq <- update(fit_sl_DEC, control = lmmControl(parallelphi = FALSE))

# comparing sequential and parallel times for DEC
fit_sl_DEC$elapsedTime
fit_sl_DEC_seq$elapsedTime

# comparing the fitted models
criteria(list(UNC = fit_sl, AR1 = fit_sl_ar1, AR2 = fit_sl_ar2,
              DEC = fit_sl_DEC)) %>% kable(digits = 3)

# estimates for the AR(2)-SL-LMM
rbind(fit_sl_ar2$theta,fit_sl_ar2$std.error) %>% kable(digits = 3)

# computing ACFs
grid.arrange(plot(acfresid(fit_sl, calcCI = TRUE, maxLag = 6)),
             plot(acfresid(fit_sl_ar1, calcCI = TRUE, maxLag = 6)), nrow = 1)

# summary for the AR(1)-SL-LMM
summary(fit_sl_ar1)

# Mahalanobis distance for the AR(1)-SMN-LMM
grid.arrange(plot(mahalDist(fit_sl_ar1), nlabels = 0) + ggtitle(NULL),
             qplot(mahalDist(fit_sl_ar1), fit_sl_ar1$uhat, shape=I(1)) +
               geom_point(shape=1)+ theme_minimal() +
               ylab("weight") + xlab("Mahalanobis distance"), ncol=2)

# plotting the AR(1)-SMN-LMM
plot(fit_sl_ar1,type = "normalized")

# additional options: diagonal random effects scale matrix
fit_sl_ar1D <- update(fit_sl_ar1, covRandom = "pdDiag")
lr.test(fit_sl_ar1, fit_sl_ar1D)

# additional options: forcing a subset of lambda to zero; using EM algorithm
fit_ssl1 <- update(fit_ssl, skewind = c(1, 0), control = lmmControl(algorithm = "EM"))
lr.test(fit_ssl1, fit_sl)

# additional options: changing initial values
fit_sl2 <- update(fit_sl, control = lmmControl(initialValues = list(nu = 1)))

# computing parametric bootstrap (takes a while to run)
boot_sl_ar1<-boot_par(fit_sl_ar1, B = 100)
boot_ci(boot_sl_ar1) %>% kable(digits = 2)
