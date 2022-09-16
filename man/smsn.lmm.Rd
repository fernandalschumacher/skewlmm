\name{smsn.lmm}
\alias{smsn.lmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ML estimation of scale mixture of skew-normal linear mixed models
}
\description{
It fits a scale mixture of skew-normal linear mixed model with possible within-subject dependence structure, using an EM-type algorithm (via the DAAREM method, by default).
It provides estimates and standard errors of parameters.
}
\usage{
smsn.lmm(data, formFixed, groupVar, formRandom = ~1, depStruct = "UNC",
         timeVar = NULL, distr = "sn", covRandom = "pdSymm",
         skewind, pAR = 1, control = lmmControl())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data frame containing the variables named in \code{formFixed}, \code{formRandom}, \code{groupVar}, and \code{timeVar}.}
  \item{formFixed}{A two-sided linear formula object describing the fixed effects part of the model, with the response on the left of a \code{~} operator and the covariates, separated by \code{+} operators, on the right.}
  \item{groupVar}{A character containing the name of the variable which represents the subjects or groups in \code{data}.}
  \item{formRandom}{A one-sided linear formula object describing the random effects part of the model, with the covariates, separated by \code{+} operators, on the right of a \code{~} operator. By default, a model with random intercept is considered.}
  \item{depStruct}{A character indicating which dependence structure should be used. \code{"UNC"} for conditionally uncorrelated (\code{"CI"} is also accepted), \code{"ARp"} for AR(p) -- p is \code{length(phi)}--, \code{"CS"} for compound symmetry, \code{"DEC"} for DEC, and \code{"CAR1"} for  continuous-time AR(1).}
  \item{timeVar}{A character containing the name of the variable which represents the time in \code{data}. Meaningless if \code{depStruct="UNC"} or \code{depStruct="CS"}. For other structures, if \code{is.null(timeVar)} the observations are considered equally spaced and ordered. If \code{depStruct="ARp"}, timeVar must be an index, preferably starting at 1.}
  \item{distr}{A character indicating which distribution should be used. \code{"sn"} for skew-normal, \code{"st"} for skew-t, \code{"ssl"} for skew-slash, and \code{"scn"} for skew-contaminated normal.}
  \item{covRandom}{A character indicating which structure should be used for the random
effects scale matrix (either \code{"pdSymm"} (default), for a general positive-definite matrix, or \code{"pdDiag"}, for a diagonal matrix).}
  \item{skewind}{A vector of length equal to the number of random effects, containing
0's and 1's, indicating which elements of the skewness parameter vector should be estimated
(optional, default is an all-ones vector).}
  \item{pAR}{If \code{depStruct="ARp"}, pAR indicates the order of the autoregressive process that should be used. Otherwise, it is meaningless.}
  \item{control}{An object resulting from the function \code{lmmControl()}, containing additional options for the estimation algorithm.}
  }
\details{
It fits the model \eqn{Y_i = X_i \beta + Z_i b_i + \epsilon_i}, for \eqn{i=1,\dots,n}, where \eqn{Y_i} is a vector with \eqn{n_i} observed continuous responses, \eqn{b_i ~ SMSN(c\Delta, D,\lambda;H)} and \eqn{\epsilon_i ~ SMN(0, \Sigma_i; H)}, indexed by the same mixing distribution.

For efficiency, the DAAREM method is used for parameter estimation. In case of numerical errors, please try passing to  \code{\link{lmmControl}(algorithm = "DAAREM")} to the \code{control} argument.

For details see Schumacher, Lachos & Matos (2021).
}
\value{
An object of class "SMSN" representing the SMSN-LMM fit. Generic functions such as print and summary have methods to show the results of the fit. The functions fitted and ranef can be used to extract some of its components.

Specifically, the following components are returned:
  \item{theta}{Named vector with parameter estimates.}
  \item{iter}{Number of iterations runned.}
  \item{estimates}{A named list containing parameter estimates.}
  \item{uhat}{Estimated weights.}
  \item{loglik.track}{Vector containing the log-likelihood value from each iteration of the estimation procedure.}
  \item{random.effects}{Estimated random effects.}
  \item{std.error}{A vector with standard errors.}
  \item{loglik}{Value of the log-likelihood at last iteration.}
  \item{elapsedTime}{Time elapsed in processing, in seconds.}
  \item{error}{Convergence criterion at last iteration.}
  \item{call}{The \code{smsn.lmm} call that produced the object.}
  \item{criteria}{A list with AIC and BIC criterion.}
  \item{data}{The data frame used on \code{smsn.lmm} call.}
  \item{formula}{A list containing the formulas used on \code{smsn.lmm} call.}
  \item{depStruct}{A character indicating which dependence structure was used.}
  \item{covRandom}{A character indicating which structure was used for the random effects scale matrix.}
  \item{distr}{A character indicating which distribution was used.}
  \item{N}{The number of observations used.}
  \item{n}{The number of individuals/groups used.}
  \item{groupVar}{A character indicating the name of the grouping variable.}
  \item{control}{The control list used in the function's call.}
  \item{timeVar}{A character indicating the name of the time variable, if any.}
  \item{fitted}{A vector of fitted values, if \code{calc.bi=TRUE}.}
}
\references{
Henderson, N. C. and R. Varadhan (2019). Damped anderson acceleration with restarts
and monotonicity control for accelerating EM and EM-like algorithms. Journal of
Computational and Graphical Statistics 28(4), 834-846.

Lachos, V. H., P. Ghosh, and R. B. Arellano-Valle (2010). Likelihood based inference for skew-normal independent linear mixed models. Statistica Sinica 20, 303-322.

Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Statistics in Medicine 40(7), 1790-1810.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{lmmControl}}, \code{\link{update}}, \code{\link{predict.SMSN}}, \code{\link{residuals.SMSN}}, \code{\link{plot.SMSN}}, \code{\link{smn.lmm}}, \code{\link{smn.clmm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#simple example
dat1 = as.data.frame(nlme::Orthodont)
fm1 = smsn.lmm(dat1, formFixed=distance ~ age, groupVar="Subject",
                control=lmmControl(max.iter=30))
fm1
\donttest{
#fitting for several distributions / dependence structures
fm1 = smsn.lmm(dat1, formFixed=distance ~ age+Sex, groupVar="Subject")
fm2 = smsn.lmm(dat1, formFixed=distance ~ age+Sex, groupVar="Subject",
               distr="st")
fm3 = smsn.lmm(dat1, formFixed=distance ~ age+Sex, groupVar="Subject",
               distr="ssl")
fm4 = smsn.lmm(dat1, formFixed=distance ~ age+Sex, groupVar="Subject",
               depStruct="ARp", pAR=1)
criteria(list(fm1=fm1, fm2=fm2, fm3=fm3, fm4=fm4))
summary(fm3)

#some diagnostic tools
plot(fm3)
acf3 = acfresid(fm3, calcCI=TRUE, MCiter=100)
plot(acf3)
plot(mahalDist(fm3), nlabels=2)
healy.plot(fm3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
