\name{logLik}
\alias{logLik.SMN}
\alias{logLik.SMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Log-likelihood of an smn.lmm and smsn.lmm models
%%  ~~function to do ... ~~
}
\description{
It returns the log-likelihood value of the linear mixed model represented by object evaluated at the estimated coefficients.
}


\usage{
\method{logLik}{SMN}(object, \dots)
\method{logLik}{SMSN}(object, \dots)
}

\arguments{
 \item{object}{An object inheriting from class \code{SMN} or \code{SMSN}}
  \item{\dots}{Additional arguments.}
   }

\value{
The log-likelihood value of the model represented by object evaluated at the estimated coefficients.
}


\seealso{
\code{\link{smn.lmm}},
\code{\link{smsn.lmm}},
\code{\link{criteria}}
}
\examples{
fm1 = smn.lmm(distance ~ age+Sex, data=nlme::Orthodont, groupVar="Subject")
logLik(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
