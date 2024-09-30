\name{coef}
\alias{coef}
\alias{coef.SMN}
\alias{coef.SMSN}
\alias{coef.SMNclmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract coefficients from smsn.lmm, smn.lmm and smn.clmm objects
%%  ~~function to do ... ~~
}
\description{
  It extracts estimated coefficients from smsn.lmm, smn.lmm and smn.clmm objects. The estimated coefficients are obtained by adding together the fixed and random effects estimates.}

\usage{
\method{coef}{SMN}(object, \dots)
\method{coef}{SMSN}(object, \dots)
\method{coef}{SMNclmm}(object, \dots)
}

\arguments{
 \item{object}{An object inheriting from class \code{SMN}, \code{SMSN}, or \code{SMNclmm}, representing
   a fitted scale mixture skew-normal linear mixed model.}
  \item{\dots}{Additional arguments}
   }

\value{
  Matrix of coefficients.
}

\references{
Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Statistics in Medicine 40(7), 1790-1810.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smsn.lmm}}, \code{\link{smn.lmm}}, \code{\link{smn.clmm}}, \code{\link{fitted.SMSN}}, \code{\link{fitted.SMN}}, \code{\link{fitted.SMNclmm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
fm1 = smn.lmm(distance ~ age+Sex, data=nlme::Orthodont, groupVar="Subject")
coef(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
