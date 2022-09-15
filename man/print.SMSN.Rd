\name{print.SMSN}
\alias{print.SMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print a smsn.lmm object
%%  ~~function to do ... ~~
}
\description{
Print a smsn.lmm object.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{print}{SMSN}(x, \dots)
}

\arguments{
 \item{x}{An object inheriting from class \code{SMSN}, representing
   a fitted scale mixture skew-normal linear mixed model.}
  \item{\dots}{Additional print arguments.}
   }

\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smsn.lmm}}, \code{\link{summary.SMSN}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{\donttest{
fm1 = smsn.lmm(distance ~ age+Sex, data=nlme::Orthodont, groupVar="Subject")
fm1
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
