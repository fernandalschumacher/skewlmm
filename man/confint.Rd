\name{confint}
\alias{confint}
\alias{confint.SMN}
\alias{confint.SMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes confidence intervals from smn.lmm and smsn.lmm fitted models
%%  ~~function to do ... ~~
}
\description{
Computes either asymptotical (based on normality from maximum likelihood estimators), using a sandwich method, or parametric bootstraped confidence intervals from a model fit.}


\usage{
\method{confint}{SMN}(object, parm, level = 0.95,  method = "asymptotic",
        parallel = NULL, seed = 123, \dots)
\method{confint}{SMSN}(object, parm, level = 0.95,  method = "asymptotic",
        parallel = NULL, seed = 123, \dots)
}

\arguments{
 \item{object}{An object inheriting from class \code{SMN} or \code{SMSN}}
 \item{parm}{A character indicating for which parameter the intervals should be returned. Available options: "beta" for fixed effects, or "all" for all parameters. Default is "all" for asymptotic and bootstrap methods, and "beta" for sandwich method.}
 \item{level}{Confidence level to be used.}
 \item{method}{A character indicating with method should be used. "asymptotic" refers to traditional confidence intervals based on asymptotical normality from maximum likelihood estimators; "bootstrap" performs a parametric bootstrap method based on B samples (100 by default), and is only recommended to small to moderate sample sizes; "sandwich" uses a sandwich variance estimator to compute the standard errors and is only available for the fixed-effects (\code{parm = "beta"}).}
  \item{\dots}{Additional arguments to be passed to \code{boot_par}.}
   }

\value{
  A table containing the estimate and the respective confidence interval.
}
\references{
Travassos, K., Matos, L. A., and Schumacher, F. L. (2026+). Sandwich variance estimation in scale mixture of skew-normal linear mixed models. Under review.

Schumacher, F. L., Matos, L. A., and Lachos, V. H. (2025). "skewlmm: An R Package for Fitting Skewed and Heavy-Tailed Linear Mixed Models." \emph{Journal of Statistical Software}, \strong{115}(7), 1--32.

Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. \emph{Statistics in Medicine} 40(7), 1790-1810.
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

\donttest{#Sandwich variance estimator
confint(fm1, method = "sandwich")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
