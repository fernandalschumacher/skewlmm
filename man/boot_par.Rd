\name{boot_par}
\alias{boot_par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parametric bootstrap for SMSN/SMN objects
%%  ~~function to do ... ~~
}
\description{
  It generates and estimates \code{B} Monte Carlo samples identical to the fitted model.
}
\usage{
boot_par(object, B = 100, seed = 123)
}
\arguments{
  \item{object}{A smsn.lmm or smn.lmm object containing the fitted model to be updated.}
  \item{B}{Number of samples to be used.}
  \item{seed}{Seed to be used.}
}
\value{
  A tibble of class \code{lmmBoot} with \code{B} lines, where each line contains the estimated parameters from a simulated sample.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}

\details{
This function provides an alternative for the asymptotic standard errors and confidence intervals given in \code{summary}, which may be helpful for small samples.
Nevertheless, the computational cost is higher and it may take several minutes to get the results.
}
\seealso{
\code{\link{boot_ci}}, \code{\link{smsn.lmm}}, \code{\link{smn.lmm}}
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
