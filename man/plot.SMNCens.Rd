\name{plot.SMNclmm}
\alias{plot.SMNclmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a smn.clmm object
  %%  ~~function to do ... ~~
}
\description{
  Fitted values versus residuals plot. Censored values are imputed using their conditional expectation from the fitting algorithm.
}

\usage{
  \method{plot}{SMNclmm}(x, level = "conditional", useweight = TRUE,
                         alpha = 0.3, \dots)
}

\arguments{
  \item{x}{An object inheriting from class \code{SMNclmm}, representing
    a fitted scale mixture of normal censored linear mixed model.}
  \item{level}{Level of residual that should be used. For details see \code{\link{residuals.SMNclmm}}. Default is "conditional".}
  \item{useweight}{A logical value indicating if the estimated weights should be used as color in the plot.}
  \item{alpha}{Transparency parameter to be used (\code{0<alpha<1}). Meaningless if \code{useweight = TRUE}.}
  \item{\dots}{Additional arguments.}
}

\value{
  A ggplot object.
}

\author{
  Fernanda L. Schumacher, Larissa A. Matos, Victor H. Lachos and Katherine L. Valeriano
}


\seealso{
\link{ggplot}, \link{smn.clmm}, \link{fitted.SMNclmm}, \link{residuals.SMNclmm}
}
\examples{\donttest{nj1 = 5; m = 30
time = rep(1:nj1, times=m)
groups = as.factor(rep(1:m, each=nj1))
dat1 = rsmsn.clmm(time, groups, cbind(1,time), rep(1,m*nj1), sigma2=0.7,
                  D=0.5*diag(1), beta=c(1,2), depStruct="UNC", distr="t", nu=4)
# Estimation
fm1 = smn.clmm(dat1, formFixed=y~x, groupVar="ind", depStruct="UNC", ci="ci",
               lcl="lcl", ucl="ucl", distr="t", control=lmmControl(max.iter=30))
plot(fm1)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
