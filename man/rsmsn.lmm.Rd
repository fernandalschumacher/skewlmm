\name{rsmsn.lmm}
\alias{rsmsn.lmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate data from SMSN-LMM
%%  ~~function to do ... ~~
}
\description{
It creates a simulated data set from SMSN-LMM (or from SMN-LMM, if \code{lambda=0}) with several possible dependence structures, for one subject.
}
\usage{
rsmsn.lmm(time1, x1, z1, sigma2, D1, beta, lambda, depStruct = "UNC",
    phi = NULL, distr = "sn", nu = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time1}{Vector containing times that should be used in data generation.}
  \item{x1}{Design matrix for fixed effects.}
  \item{z1}{Design matrix for random effects.}
  \item{sigma2}{Common variance parameter, such that \eqn{\Sigma=sigma2*R}.}
  \item{D1}{Variance matrix for random effects.}
  \item{beta}{Vector of fixed effects parameter.}
  \item{lambda}{Skewness parameter of random effects.}
  \item{depStruct}{Dependence structure. \code{"UNC"} for conditionally uncorrelated (\code{"CI"} is also accepted), \code{"ARp"} for AR(p) -- p is \code{length(phi)}--, \code{"CS"} for compound symmetry, \code{"DEC"} for DEC, and \code{"CAR1"} for  continuous-time AR(1).}
  \item{phi}{Parameter vector indexing the dependence structure.}
  \item{distr}{Distribution that should be used. \code{"sn"} for skew-normal, \code{"st"} for skew-t, \code{"ss"} for skew-slash, and \code{"scn"} for skew-contaminated normal.}
  \item{nu}{Parameter vector indexing \code{distr}. Should be \code{NULL} for \code{"sn"}, be a vector of length 1 for \code{"st"} and \code{"ss"}, and of length 2 for \code{"scn"}.}
}

\value{
  A data frame containing time, the generated response variable (\code{y}), and possible covariates.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Lachos, V. H., P. Ghosh, and R. B. Arellano-Valle (2010). Likelihood based inference for skew-normal independent linear mixed models. Statistica Sinica 20, 303-322.

Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Statistics in Medicine. \doi{10.1002/sim.8870}.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smsn.lmm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#generating a sample for 1 individual at 5 times
nj1<-5
rsmsn.lmm(1:nj1,cbind(1,1:nj1),rep(1,nj1),sigma2=.25,
          D1=diag(1),beta = c(1,2),lambda = 2,depStruct = "ARp",phi=.5,
          distr="st",nu=5)

#generating a sample for m=20 individuals at 5 times
library(dplyr)
library(purrr)
library(ggplot2)
nj1<-5
m<-50
gendatList<- map(rep(nj1,m), function(nj) rsmsn.lmm(1:nj,cbind(1,1:nj),rep(1,nj),
  sigma2=.25,D1=.5*diag(1),beta = c(1,2),lambda = 2,depStruct = "ARp",phi=.5))
gendat<- bind_rows(gendatList,.id="ind")
ggplot(gendat,aes(x=x,y=y,group=ind)) + geom_line() +
  stat_summary(aes(group = 1),geom = "line", fun.y = mean, col="blue",size=2)
#
\donttest{fm1 <- smsn.lmm(y ~ x, data = gendat, groupVar="ind", depStruct="ARp", pAR=1)
summary(fm1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}% use one of  RShowDoc("KEYWORDS")
