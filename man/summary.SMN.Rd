\name{summary.SMN}
\alias{summary.SMN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of a smn.lmm object
%%  ~~function to do ... ~~
}
\description{
  \code{summary} method for class "SMN".
}

\usage{
\method{summary}{SMN}(object, confint.level = 0.95, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{SMN}, representing
   a fitted scale mixture normal linear mixed model.}
  \item{confint.level}{Level of the approximate confidence intervals presented. Default=0.95.}
  \item{\dots}{Additional arguments.}
   }

\value{
  \item{varRandom}{estimated variance matrix from random effects (\eqn{D}).}
  \item{varFixed}{parameter estimates of variance from random errors (\eqn{\Sigma}). For recovering the error variance matrix use \code{\link{errorVar}} function.}
  \item{tableFixed}{estimated fixed effects, their standard errors and approximated confidence intervals.}
  \item{criteria}{Log-likelihood value from MLE, AIC and BIC criteria.}
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smn.lmm}}, \code{\link{errorVar}}, \code{\link{plot.SMN}}, \code{\link{residuals.SMN}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
summary(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
