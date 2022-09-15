\name{print.SMNCens}
\alias{print.SMNCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print a smnCens.lmm object
  %%  ~~function to do ... ~~
}
\description{
  Print a smnCens.lmm object.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  \method{print}{SMNCens}(x, \dots)
}

\arguments{
  \item{x}{An object inheriting from class \code{SMNCens}, representing
    a fitted scale mixture normal linear mixed model with censored responses.}
  \item{\dots}{Additional print arguments.}
}

\author{
  Fernanda L. Schumacher, Larissa A. Matos, Katherine L. Valeriano, and Victor H. Lachos
}
\seealso{
  \code{\link{smnCens.lmm}}, \code{\link{summary.SMNCens}}
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
nj1 = 5; m = 30
time = rep(1:nj1, times=m)
groups = as.factor(rep(1:m, each=nj1))
dat1 = rsmnCens.lmm(time, groups, cbind(1,time), rep(1,m*nj1), sigma2=0.7,
                    D=diag(1), beta=c(1,2), depStruct="UNC")
fm1 = smnCens.lmm(dat1, formFixed=y~x, groupVar="ind", formRandom=~1,
                  depStruct="UNC", ci="ci", lcl="lcl", ucl="ucl")
fm1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
