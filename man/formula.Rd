\name{formula}
\alias{formula.SMN}
\alias{formula.SMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Formula from an smn.lmm and smsn.lmm models
%%  ~~function to do ... ~~
}
\description{
It returns the formula used for both fixed and random terms of the linear mixed model represented by object.
}


\usage{
\method{formula}{SMN}(x, \dots)
\method{formula}{SMSN}(x, \dots)
}

\arguments{
 \item{x}{An object inheriting from class \code{SMN} or \code{SMSN}}
  \item{\dots}{Additional arguments}
   }

\value{
\item{formFixed}{Fixed effects formula}
\item{formRandom}{Random effects formula}
\item{groupVar}{Variable identified subjects or clusters}
}


\seealso{
\code{\link{smn.lmm}},
\code{\link{smsn.lmm}},
\code{\link{criteria}}
}
\examples{
fm1 = smn.lmm(distance ~ age+Sex, data=nlme::Orthodont, groupVar="Subject")
formula(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
