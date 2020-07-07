\name{smsn.lmm}
\alias{smsn.lmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ML estimation of scale mixture of skew-normal linear mixed models
}
\description{
It fits a scale mixture of skew-normal linear mixed model with possible within-subject dependence structure, using an EM-type algorithm.
It provides estimates and standard errors of parameters.
}
\usage{
smsn.lmm(data, formFixed, groupVar, formRandom = ~1, depStruct = "CI",
    timeVar = NULL, distr = "sn", pAR = 1, luDEC = 10,
    tol = 1e-06, max.iter = 200, calc.se = TRUE, lb = NULL,
    lu = NULL, initialValues = list(beta = NULL, sigma2 = NULL,
      D = NULL, lambda = NULL, phi = NULL, nu = NULL), quiet = FALSE, showCriterium=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data frame containing the variables named in \code{formFixed}, \code{formRandom}, \code{groupVar}, and \code{timeVar}.}
  \item{formFixed}{A two-sided linear formula object describing the fixed effects part of the model, with the response on the left of a \code{~} operator and the covariates, separated by \code{+} operators, on the right.}
  \item{groupVar}{A character containing the name of the variable which represents the subjects or groups in \code{data}.}
  \item{formRandom}{A one-sided linear formula object describing the random effects part of the model, with the covariates, separated by \code{+} operators, on the right of a \code{~} operator. By default, a model with random intercept is considered.}
  \item{depStruct}{A character indicating which dependence structure should be used. \code{"CI"} for condicional independence, \code{"ARp"} for AR(p) -- p is \code{length(phi)}--, \code{"CS"} for compound symmetry, \code{"DEC"} for DEC, and \code{"CAR1"} for  continuous-time AR(1).}
  \item{timeVar}{A character containing the name of the variable which represents the time in \code{data}. Meaningless if \code{depStruct="CI"} or \code{depStruct="CS"}. For other structures, if \code{is.null(timeVar)} the observations are considered equally spaced and ordered. If \code{depStruct="ARp"}, timeVar must the an index, preferably starting at 1.}
  \item{distr}{A character indicating which distribution should be used. \code{"sn"} for skew-normal, \code{"st"} for skew-t, \code{"ssl"} for skew-slash, and \code{"scn"} for skew-contaminated normal.}
  \item{pAR}{If \code{depStruct="ARp"}, pAR indicates the order of the autoregressive process that should be used. Otherwise, it is meaningless.}
  \item{luDEC}{Optional. Upper limit for estimating the "damping" parameter for DEC covariance. If \code{luDEC<=1}, only attenuation of the exponential decay can be obtained.}
  \item{tol}{Tolerance for the convergence criterion. Default=\code{1e-6}.}
  \item{max.iter}{Maximum number of iterations for the EM algorithm. Default=\code{200}.}
  \item{calc.se}{A logical value indicating if standard errors should be calculated.}
  \item{lb}{Optional. Bottom limit for estimating \code{nu}.}
  \item{lu}{Optional. Upper limit for estimating \code{nu}.}
  \item{initialValues}{Optional. A named list containing initial parameter values, with at most the following elements: \code{beta}, \code{sigma2}, \code{D}, \code{lambda}, \code{phi}, \code{nu}.}
  \item{quiet}{A logical value indicating if the iteration message should be suppressed. Useful when calling the function in R Markdown.}
  \item{showCriterium}{A logical value indicating if the criterium should be shown at each iteration.}
  }
\details{
It fits the model \eqn{Y_i = X_i \beta + Z_i b_i + \epsilon_i}, for \eqn{i=1,\dots,n}, where \eqn{Y_i} is a vector with \eqn{n_i} observed continuous responses, \eqn{b_i ~ SMSN(c\Delta, D,\lambda;H)} and \eqn{\epsilon_i ~ SMN(0, \Sigma_i; H)}, indexed by the same mixing distribution.

For details see Schumacher, Matos & Lachos (2020+).
}
\value{
An object of class "SMSN" representing the SMSN-LMM fit. Generic functions such as print and summary have methods to show the results of the fit. The functions fitted and ranef can be used to extract some of its components.

Specifically, the following components are returned:
  \item{theta}{Named vector with parameter estimates.}
  \item{iter}{Number of iterations runned.}
  \item{estimates}{A named list containing parameter estimates.}
  \item{random.effects}{Estimated random effects.}
  \item{std.error}{A vector with standard errors.}
  \item{loglik}{Value of the log-likelihood at last iteration.}
  \item{elapsedTime}{Time elapsed in processing, in seconds.}
  \item{error}{Convergence criterion at last iteration.}
  \item{call}{The \code{smsn.lmm} call that produced the object.}
  \item{criteria}{A list with AIC and BIC criterion.}
  \item{data}{The data frame used on \code{smsn.lmm} call.}
  \item{forumla}{A list containing the formulas used on \code{smsn.lmm} call.}
  \item{depStruct}{A character indicating which depedence structure was used.}
  \item{distr}{A character indicating which distribution was used.}
  \item{N}{The number of observations used.}
  \item{n}{The number of individuals/groups used.}
  \item{groupVar}{A character indicating the name of the grouping variable.}
  \item{timeVar}{A character indicating the name of the time variable, if any.}
  \item{fitted}{A vector of fitted values, if \code{calc.bi=TRUE}.}
}
\references{
Lachos, V. H., P. Ghosh, and R. B. Arellano-Valle (2010). Likelihood based inference for skew-normal independent linear mixed models. Statistica Sinica 20, 303-322.

Schumacher, F.L., V.H. Lachos, and L.A. Matos (2020+). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Submitted.
Preprint available at \url{https://arxiv.org/abs/2002.01040}.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{print.SMSN}}, \code{\link{summary.SMSN}}, \code{\link{predict.SMSN}}, \code{\link{residuals.SMSN}}, \code{\link{plot.SMSN}}, \code{\link{smn.lmm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#simple example
  dat1 <- as.data.frame(nlme::Orthodont)
  fm1 <- smsn.lmm(dat1,formFixed=distance ~ age,groupVar="Subject",max.iter = 30)
  fm1
\donttest{#fitting for several distributions / dependence structures
  fm1 <- smsn.lmm(dat1,formFixed=distance ~ age+Sex,groupVar="Subject")
  fm2 <- smsn.lmm(dat1,formFixed=distance ~ age+Sex,groupVar="Subject", distr="st")
  fm3 <- smsn.lmm(dat1,formFixed=distance ~ age+Sex,groupVar="Subject", distr="ssl")
  fm4 <- smsn.lmm(dat1,formFixed=distance ~ age+Sex,groupVar="Subject", depStruct="ARp",pAR=1)
  rbind(fm1$criteria,fm2$criteria,fm3$criteria,fm4$criteria)
  summary(fm3)

  #some diagnostic tools
  plot(fm1)
  acf1<- acfresid(fm1,calcCI=TRUE,MCiter=100)
  plot(acf1)
  plot(mahalDist(fm1),fm1,nlabels=2)
  healy.plot(fm1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
