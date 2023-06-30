\name{smn.clmm}
\alias{smn.clmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ML estimation of scale mixture of normal linear mixed models with censored responses
}
\description{
It fits left, right, or interval censored scale mixture of normal linear mixed model with possible within-subject dependence structure, using the EM algorithm.
It provides estimates and standard errors of parameters.
}
\usage{
smn.clmm(data, formFixed, groupVar, formRandom = ~1, depStruct = "UNC",
         ci, lcl, ucl, timeVar = NULL, distr = "norm",
         nufix = FALSE, pAR = 1, control = lmmControl())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data frame containing the variables named in \code{formFixed}, \code{formRandom}, \code{groupVar}, \code{timeVar}, \code{ci}, \code{lcl}, and \code{ucl}.}
  \item{formFixed}{A two-sided linear formula object describing the fixed effects part of the model, with the response on the left of a \code{~} operator and the covariates, separated by \code{+} operators, on the right.}
  \item{groupVar}{A character containing the name of the variable which represents the subjects or groups in \code{data}.}
  \item{formRandom}{A one-sided linear formula object describing the random effects part of the model, with the covariates, separated by \code{+} operators, on the right of a \code{~} operator. By default, a model with random intercept is considered.}
  \item{depStruct}{A character indicating which dependence structure should be used. \code{"UNC"} for conditionally uncorrelated (\code{"CI"} is also accepted), \code{"ARp"} for AR(p) -- p is \code{length(phi)}--, \code{"CS"} for compound symmetry, \code{"DEC"} for DEC, \code{"CAR1"} for continuous-time AR(1), and \code{"MA1"} for moving average of order 1.}
  \item{ci}{A character containing the name of the censoring indicator variable in \code{data}, which should be 1 if the respective observation is censored or missing, and 1 otherwise. If missing, it is assumed that none of the observations is censored.}
  \item{lcl}{A character containing the name of the lower censoring limit in \code{data}. If missing, it is assumed \code{lcl=-Inf}, i.e., no left limit.}
  \item{ucl}{A character containing the name of the upper censoring limit in \code{data}. If missing, it is assumed \code{ucl=Inf}, i.e., no right limit.}
  \item{timeVar}{A character containing the name of the variable which represents the time in \code{data}. Meaningless if \code{depStruct="UNC"} or \code{depStruct="CS"}. For other structures, if \code{is.null(timeVar)} the observations are considered equally spaced and ordered. If \code{depStruct="ARp"}, timeVar must be an index, preferably starting at 1.}
  \item{distr}{A character indicating which distribution should be used. \code{"norm"} for normal and \code{"t"} for Student-t.}
  \item{nufix}{\code{TRUE} or \code{FALSE} indicating if \code{nu} should be estimated for t distribution. If \code{nufix=TRUE}, nu must be especified through \code{lmmControl()}.}
  \item{pAR}{If \code{depStruct="ARp"}, \code{pAR} indicates the order of the autoregressive process that should be used (1 by default). Otherwise, it is meaningless.}
  \item{control}{An object resulting from the function \code{lmmControl()}, containing additional options for the estimation algorithm.}
}
\details{
  It fits the model \eqn{Y_i = X_i \beta + Z_i b_i + \epsilon_i}, for \eqn{i=1,\dots,n}, where \eqn{Y_i} is a vector with \eqn{n_i} observed continuous responses, \eqn{b_i ~ SMN(0, D; H)} and \eqn{\epsilon_i ~ SMN(0, \Sigma_i; H)}, indexed by the same mixing distribution.

  For details see Matos et al. (2013) and Lachos et al. (2019).
}
\value{
  An object of class "SMNclmm" representing the SMN-CLMM fit. Generic functions such as print and summary have methods to show the results of the fit. The functions fitted and ranef can be used to extract some of its components.

  Specifically, the following components are returned:
  \item{theta}{Named vector with parameter estimates.}
  \item{iter}{Number of iterations runned.}
  \item{estimates}{A named list containing parameter estimates.}
  \item{yest}{Imputed values in the response variable.}
  \item{uhat}{Estimated weights.}
  \item{loglik.track}{Vector containing the log-likelihood value from each iteration of the estimation procedure.}
  \item{random.effects}{Estimated random effects.}
  \item{std.error}{A vector with standard errors.}
  \item{loglik}{Value of the log-likelihood at last iteration.}
  \item{elapsedTime}{Time elapsed in processing, in seconds.}
  \item{error}{Convergence criterion at last iteration.}
  \item{criteria}{A list with AIC, BIC, and SIC criterion.}
  \item{call}{The \code{smn.clmm} call that produced the object.}
  \item{data}{The data frame used on \code{smn.clmm} call.}
  \item{formula}{A list containing the formulas used on \code{smn.clmm} call.}
  \item{depStruct}{A character indicating which dependence structure was used.}
  \item{covRandom}{A character indicating which structure was used for the random effects scale matrix.}
  \item{distr}{A character indicating which distribution was used.}
  \item{N}{The number of observations used.}
  \item{n}{The number of individuals/groups used.}
  \item{groupVar}{A character indicating the name of the grouping variable.}
  \item{timeVar}{A character indicating the name of the time variable, if any.}
  \item{fitted}{A vector of fitted values.}
}
\references{
Henderson, N. C. and R. Varadhan (2019). Damped anderson acceleration with restarts
and monotonicity control for accelerating EM and EM-like algorithms. Journal of
Computational and Graphical Statistics 28(4), 834-846.

Matos, L. A., Prates, M. O., Chen, M. H., and Lachos, V. H. (2013). Likelihood-based inference for mixed-effects models with censored response using the multivariate-t distribution. Statistica Sinica 23(3), 1323-1345.

Lachos, V. H., A. Matos, L., Castro, L. M., and Chen, M. H. (2019). Flexible longitudinal linear mixed models for multiple censored responses data. Statistics in medicine, 38(6), 1074-1102.
}
\author{
Larissa A. Matos, Victor H. Lachos, Katherine L. Valeriano and Fernanda L. Schumacher
}

\seealso{
\code{\link{lmmControl}}, \code{\link{update}}, \code{\link{predict.SMNclmm}}, \code{\link{residuals.SMNclmm}}, \code{\link{plot.SMNclmm}},
\code{\link{smn.lmm}}, \code{\link{smsn.lmm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Generating a sample for m=30 individuals with 5 times
# Considering 10\% of left censoring responses
# AR(1) and normal distribution
nj1 = 5
m = 30
time = rep(1:nj1, times=m)
groups = as.factor(rep(1:m, each=nj1))
dat1 = rsmsn.clmm(time = time, ind = groups, x = cbind(1,time), z = rep(1,m*nj1),
                  sigma2=0.6, D=0.5*diag(1), beta=c(1,2), depStruct="ARp",
                  phi=0.4, pcens = .1, type = "left")

# Estimation using UNC
fm0 = smn.clmm(dat1, formFixed=y~x, groupVar="ind", ci="ci", ucl="ucl")
summary(fm0)
\donttest{
# Estimation using AR(1)
fm1 = update(fm0, depStruct="ARp")

# Estimation using AR(1) and t distribution
fm2 = update(fm1, distr="t")

# Comparing fitted models
criteria(list(fm0= fm0, fm1=fm1, fm2=fm2))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
