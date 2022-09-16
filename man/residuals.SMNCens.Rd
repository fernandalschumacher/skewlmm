\name{residuals.SMNCens}
\alias{residuals.SMNCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract model residuals from smn.clmm objects
  %%  ~~function to do ... ~~
}
\description{
  The conditional residuals are obtained by subtracting the fitted values from the response vector, while the marginal residuals are obtained by subtracting only the fixed effects from the response vector.}

\usage{
  \method{residuals}{SMNCens}(object, level = "conditional", \dots)
}

\arguments{
  \item{object}{An object inheriting from class \code{SMNCens}, representing
    a fitted scale mixture of normal censored linear mixed model.}
  \item{level}{Either \code{"conditional"}, for obtaining conditional residuals, or \code{"marginal"}, for marginal residuals.}
  \item{\dots}{Additional arguments.}
}

\value{
  Vector with the residuals of length equal to \code{nrow(data)}.
}
\note{
The residuals are computed after imputing the censored observations.}
\author{
  Fernanda L. Schumacher, Larissa A. Matos, Victor H. Lachos and Katherine L. Valeriano
}

\seealso{
  \code{\link{smn.clmm}}, \code{\link{mahalDist}}
}
\examples{
nj1 = 5
m = 30
time = rep(1:nj1, times=m)
groups = as.factor(rep(1:m, each=nj1))
dat1 = rsmn.clmm(time, groups, cbind(1,time), rep(1,m*nj1), sigma2=0.7,
                 D=diag(1), beta=c(1,2), depStruct="CS", phi=0.4)

fm1 = smn.clmm(dat1, formFixed=y~x, groupVar="ind", depStruct="CS", ci="ci",
               lcl="lcl", ucl="ucl", control=lmmControl(max.iter=30))
residuals(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
