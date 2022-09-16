\name{predict.SMNCens}
\alias{predict.SMNCens}

\title{Prediction of future observations from an smn.clmm object
}
\description{
  Predicted values are obtained through conditional expectation. For details, see Schumacher, Lachos, and Matos (2021).
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  \method{predict}{SMNCens}(object, newData, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object inheriting from class \code{SMNCens}, representing
    a fitted scale mixture normal linear mixed model.}
  \item{newData}{A data frame for which response variable should be predicted, including covariates, groupVar and possibly timeVar. If missing or \code{NULL}, fitted values are returned.}
  \item{\dots}{Additional arguments.}
}

\value{
  A data frame with covariates, groupVar and ypred, where ypred contains the predicted values from the response variable.
}

\references{
  Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Statistics in Medicine 40(7), 1790-1810.
}
\author{
  Fernanda L. Schumacher, Larissa A. Matos, Victor H. Lachos and Katherine L. Valeriano
}


\seealso{
  \code{\link{smn.clmm}}, \code{\link{fitted.SMNCens}}
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Simulated example: 20 individuals with 5 times for estimation and
# 1 time for prediction
set.seed(963)
nj1 = 6
m = 20
time = rep(1:nj1, times=m)
groups = as.factor(rep(1:m, each=nj1))
dat1 = rsmn.clmm(time, groups, cbind(1,time), rep(1,m*nj1), sigma2=0.4,
                 D=diag(1), beta=c(1,2), depStruct="CS", phi=0.25)
# Estimation
fm1 = smn.clmm(subset(dat1,time<6), formFixed=y~x, groupVar="ind",
               depStruct="CS", ci="ci", lcl="lcl", ucl="ucl",
               control=lmmControl(max.iter=30, tol=1e-4))
# Prediction
pred = predict(fm1, subset(dat1,time==6))
# MSPE
mean((subset(dat1,time==6)$y - pred$ypred)^2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
