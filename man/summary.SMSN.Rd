\name{summary.SMSN}
\alias{summary.SMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of a smsn.lmm object
%%  ~~function to do ... ~~
}
\description{
  \code{summary} method for class "SMSN".
}

\usage{
\method{summary}{SMSN}(object, confint.level = 0.95, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{SMSN}, representing
   a fitted scale mixture skew-normal linear mixed model.}
  \item{confint.level}{Level of the approximate confidence intervals presented. Default=0.95.}
  \item{\dots}{Additional arguments.}
   }

\value{
  \item{varRandom}{estimated variance matrix from random effects (\eqn{D}).}
  \item{varFixed}{parameter estimates of variance from random errors (\eqn{\Sigma}). For recovering the error variance matrix use \code{\link{errorVar}} function.}
  \item{tableFixed}{estimated fixed effects, their standard errors and approximated confidence intervals.}
  \item{criteria}{Log-likelihood value from MLE, AIC and BIC criteria.}
}

\references{
Schumacher, F.L., Matos, L.A., and Lachos, V.H. (2020+) "Scale mixture of skew-normal linear mixed models with within-subject serial dependence". Submitted.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smsn.lmm}}, \code{\link{errorVar}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
    fm1 <- smsn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
    summary(fm1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
