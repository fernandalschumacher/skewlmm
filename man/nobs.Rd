\name{nobs}
\alias{nobs.SMN}
\alias{nobs.SMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract the number of observations from smn.lmm and smsn.lmm fitted models
%%  ~~function to do ... ~~
}
\description{
 Extract the total number of observations from a model fit (considering all repeated measurements from all subjects or clusters). }


\usage{
\method{nobs}{SMN}(object, \dots)
\method{nobs}{SMSN}(object, \dots)
}

\arguments{
 \item{object}{An object inheriting from class \code{SMN} or \code{SMSN}}
  \item{\dots}{Additional arguments.}
   }

\value{
  A single integer, expected to be equal to \code{nrow(data)}.
}


\seealso{
\code{\link{smn.lmm}},
\code{\link{smsn.lmm}},
\code{\link{criteria}}
}
\examples{
fm1 = smn.lmm(distance ~ age+Sex, data=nlme::Orthodont, groupVar="Subject")
nobs(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
