\name{lmmControl}
\alias{lmmControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Control options for \code{smsn.lmm()}, \code{smn.lmm()} and \code{smn.clmm()}
%%  ~~function to do ... ~~
}
\description{
The values supplied in the function call replace the defaults and a list with all possible arguments is returned. The returned list has class \code{"lmmControl"} and is used as the control argument to the \code{smsn.lmm()}, \code{smn.lmm()} and \code{smn.clmm()} functions.
}
\usage{
lmmControl(tol = 1e-06, max.iter = 300, calc.se = TRUE, lb = NULL,
           lu = NULL, luDEC = 10,
           initialValues = list(beta = NULL, sigma2 = NULL, D = NULL,
                                lambda = NULL, phi = NULL, nu = NULL),
           quiet = FALSE, showCriterium = FALSE, algorithm = "DAAREM",
           parallelphi = NULL, parallelnu = NULL, ncores = NULL,
           control.daarem = list())
}
\arguments{
  \item{tol}{Tolerance for the convergence criterion. Default = \code{1e-6}.}
  \item{max.iter}{Maximum number of iterations for the EM algorithm. Default = \code{200}.}
  \item{calc.se}{A logical value indicating if standard errors should be calculated.}
  \item{lb}{Optional. Bottom limit for estimating \code{nu}.}
  \item{lu}{Optional. Upper limit for estimating \code{nu}.}
  \item{luDEC}{Optional. Upper limit for estimating the "damping" parameter for DEC covariance. If \code{luDEC<=1}, only attenuation of the exponential decay can be obtained.}
  \item{initialValues}{Optional. A named list containing initial parameter values, with at most the following elements: \code{beta}, \code{sigma2}, \code{D}, \code{lambda}, \code{phi}, \code{nu}.}
  \item{quiet}{A logical value indicating if the iteration message should be suppressed. Useful when calling the function in R Markdown.}
  \item{showCriterium}{A logical value indicating if the criterium should be shown at each iteration.}
  \item{algorithm}{Algorithm to be used for estimation, either \code{"DAAREM"} (default) or \code{"EM"}. DAAREM is an acceleration method and usually converges with fewer iterations, but it may result in numerical errors (in this case, please use the \code{"EM"} option).}
  \item{parallelphi}{A logical value indicating if parallel optimization should be used in the numerical update of the parameters related to the within-subject dependence structure. Default is \code{TRUE} if the data contains more than \code{30} subjects, and \code{FALSE} otherwise. Meaningless if \code{depStruct = "UNC"}.}
  \item{parallelnu}{A logical value indicating if parallel optimization should be used in the numerical update of nu. Meaningless if \code{distr="norm"} or \code{distr="sn"}.}
  \item{ncores}{Number of cores to be used for the parallel optimization. Meaningless if \code{parallelphi = FALSE} and \code{parallelnu = FALSE}.}
  \item{control.daarem}{List of control for the daarem algorithm. See \code{help(daarem, package = "daarem")} for details. Meaningless if \code{algorithm = "EM"}}
}
\references{
Henderson, N.C. and Varadhan, R. (2019) Damped Anderson acceleration with restarts and monotonicity control for accelerating EM and EM-like algorithms, Journal of Computational and Graphical Statistics, Vol. 28(4), 834-846.

Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Statistics in Medicine 40(7), 1790-1810.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smsn.lmm}}, \code{\link{smn.lmm}}, \code{\link{smn.clmm}}, \code{\link{update}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
lmmControl(algorithm = "EM")

fm1 = smn.lmm(nlme::Orthodont, formFixed=distance ~ age+Sex,
              groupVar="Subject", control=lmmControl(tol=1e-7))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
