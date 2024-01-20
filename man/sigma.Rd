\name{sigma}
\alias{sigma}
\alias{sigma.SMN}
\alias{sigma.SMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Residual standard deviation from smn.lmm and smsn.lmm objects
%%  ~~function to do ... ~~
}
\description{
Extract the estimated residual standard deviation from smn.lmm and smsn.lmm objects.
}


\usage{
\method{sigma}{SMN}(object, \dots)
\method{sigma}{SMSN}(object, \dots)
}

\arguments{
 \item{object}{An object inheriting from class \code{SMN} or \code{SMSN}}
  \item{\dots}{Additional arguments.}
   }

\value{
  A positive number.
}


\seealso{
\code{\link{smn.lmm}},
\code{\link{smsn.lmm}},
\code{\link{criteria}}
}
\examples{
fm1 = smn.lmm(distance ~ age+Sex, data=nlme::Orthodont, groupVar="Subject")
sigma(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
