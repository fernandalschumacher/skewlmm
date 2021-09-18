\name{boot_ci}
\alias{boot_ci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract confidence intervals from \code{lmmBoot} object
%%  ~~function to do ... ~~
}
\description{
  It  extracts confidence intervals from parametric bootstrap results obtained using the \code{boot_par()} function.
}
\usage{
boot_ci(object, conf = 0.95)
}
\arguments{
  \item{object}{An object containing the results of the \code{boot_par()} function.}
  \item{conf}{Confidence level to be considered.}
}
\value{
  A matrix containing the confidence intervals.
  }
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{boot_par}}, \code{\link{smsn.lmm}}, \code{\link{smn.lmm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{\donttest{
fm1<- smn.lmm(nlme::Orthodont, formFixed = distance ~ age+Sex, groupVar = "Subject")
b1<- boot_par(fm1, B = 100)
boot_ci(b1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}% use one of  RShowDoc("KEYWORDS")
