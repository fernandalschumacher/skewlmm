\name{summary.SMNCens}
\alias{summary.SMNCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of a smn.clmm object
  %%  ~~function to do ... ~~
}
\description{
  \code{summary} method for class "SMNCens".
}

\usage{
  \method{summary}{SMNCens}(object, confint.level = 0.95, \dots)
}

\arguments{
  \item{object}{An object inheriting from class \code{SMNCens},
  representing a fitted scale mixture normal linear mixed model with censored responses.}
  \item{confint.level}{Level of the approximate confidence intervals presented.}
  \item{\dots}{Additional arguments.}
}
\author{
  Fernanda L. Schumacher, Larissa A. Matos, Victor H. Lachos and Katherine L. Valeriano
}
\seealso{
  \code{\link{smn.clmm}}, \code{\link{plot.SMNCens}}, \code{\link{residuals.SMNCens}}
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
nj1 = 5; m = 30
time = rep(1:nj1, times=m)
groups = as.factor(rep(1:m, each=nj1))
dat1 = rsmn.clmm(time, groups, cbind(1,time), rep(1,m*nj1), sigma2=0.7,
                 D=diag(1), beta=c(1,2), depStruct="UNC")
fm1 = smn.clmm(dat1, formFixed=y~x, groupVar="ind", formRandom=~1,
               depStruct="UNC", ci="ci", lcl="lcl", ucl="ucl")
summary(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
