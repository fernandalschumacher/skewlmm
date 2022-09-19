\name{mahalDist.SMNclmm}
\alias{mahalDist.SMNclmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mahalanobis distance from a smn.clmm object
  %%  ~~function to do ... ~~
}
\description{
  Returns the squared Mahalanobis distance from a fitted SMN-CLMM.}

\usage{
  mahalDist.SMNclmm(object)
}

\arguments{
  \item{object}{An object inheriting from class \code{SMNclmm}, representing
    a fitted scale mixture of normal censored linear mixed model.}
}

\value{
  An object of class \code{mahalDist.SMNclmm} containing the Mahalanobis distance.
}

\references{
Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Statistics in Medicine 40(7), 1790-1810.

Zeller, C. B., Labra, F. V., Lachos, V. H. & Balakrishnan, N. (2010). Influence analyses of skew-normal/independent linear mixed models. Computational Statistics & Data Analysis, 54(5).
}
\author{
Fernanda L. Schumacher, Larissa A. Matos, Victor H. Lachos and Katherine L. Valeriano
}


\seealso{
\code{\link{smn.clmm}}, \code{\link{plot.mahalDist.SMNclmm}}
}
\examples{
nj1 = 5; m = 30
time = rep(1:nj1, times=m)
groups = as.factor(rep(1:m, each=nj1))
dat1 = rsmsn.clmm(time, groups, cbind(1,time), rep(1,m*nj1), sigma2=0.7,
                  D=0.5*diag(1), beta=c(1,2), depStruct="CS", phi=0.4)
# Estimation
fm1 = smn.clmm(dat1, formFixed=y~x, groupVar="ind", depStruct="CS", ci="ci",
               lcl="lcl", ucl="ucl", control=lmmControl(max.iter=30))
distance = mahalDist.SMNclmm(fm1)
plot(distance, level=0.95, nlabels=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
