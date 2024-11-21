\name{confint}
\alias{confint}
\alias{confint.SMN}
\alias{confint.SMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes confidence intervals from smn.lmm and smsn.lmm fitted models
%%  ~~function to do ... ~~
}
\description{
Computes either asymptotical (based on normality from maximum likelihood estimators) or parametric bootstraped confidence intervals from a model fit.}


\usage{
\method{confint}{SMN}(object, param, level = 0.95, method, \dots)
\method{confint}{SMSN}(object, param, level = 0.95, method, \dots)
}

\arguments{
 \item{object}{An object inheriting from class \code{SMN} or \code{SMSN}}
 \item{param}{A character indicating for which parameter the intervals should be returned. Available options: "beta" for fixed effects, or "all" for all parameters. Default is "all".}
 \item{level}{Confidence level to be used.}
  \item{\dots}{Additional arguments to be passed to \code{boot_par}.}
   }

\value{
  A table containing the estimate and the respective confidence interval.
}


\seealso{
\code{\link{smn.lmm}},
\code{\link{smsn.lmm}},
\code{\link{boot_par}},
\code{\link{boot_ci}}
}
\examples{
fm1 = smn.lmm(distance ~ age+Sex, data=nlme::Orthodont, groupVar="Subject")
confint(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
