\name{fitted.SMNCens}
\alias{fitted.SMNCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract smn.clmm fitted values
  %%  ~~function to do ... ~~
}
\description{
  The fitted values are obtained by adding together the
  population fitted values (based only on the fixed effects estimates)
  and the estimated contributions of the random effects to the fitted
  values at grouping levels. }


\usage{
  \method{fitted}{SMNCens}(object, \dots)
}

\arguments{
  \item{object}{An object inheriting from class \code{SMNCens}, representing a fitted scale mixture normal linear mixed model with censored responses.}
  \item{\dots}{Additional arguments.}
}

\value{
  Vector of fitted values with length equal to \code{nrow(data)}.
}

\author{
  Fernanda L. Schumacher, Larissa A. Matos, Victor H. Lachos and Katherine L. Valeriano
}

\seealso{
  \code{\link{smn.clmm}}, \code{\link{predict.SMNCens}}
}
\examples{
nj1 = 5
m = 30
time = rep(1:nj1, times=m)
groups = as.factor(rep(1:m, each=nj1))
dat1 = rsmn.clmm(time, groups, cbind(1,time), rep(1,m*nj1), sigma2=0.7,
                 D=diag(1), beta=c(1,2), depStruct="UNC")

fm1 = smn.clmm(dat1, formFixed=y~x, groupVar="ind", depStruct="UNC",
               ci="ci", lcl="lcl", ucl="ucl", control=lmmControl(max.iter=30))
fitted(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
