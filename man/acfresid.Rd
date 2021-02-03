\name{acfresid}
\alias{acfresid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Autocorrelation function for smn.lmm or smsn.lmm residuals
%%  ~~function to do ... ~~
}
\description{
This function calculates the empirical autocorrelation function for the within-subject residuals from a smn.lmm or smsn.lmm  fit. The autocorrelation values are calculated using pairs of residuals within-subjects. The autocorrelation function is useful for investigating serial correlation models for discrete-time data, preferably equally spaced.}

\usage{
acfresid(object,maxLag,resLevel="marginal",resType="normalized",
         calcCI=FALSE,IClevel,MCiter,seed)
}

\arguments{
 \item{object}{An object inheriting from class \code{SMN} or \code{SMSN}, representing
   a fitted scale mixture of (skew) normal linear mixed model.}
 \item{maxLag}{An optional integer giving the maximum lag for which the autocorrelation should be calculated. Defaults to maximum lag in the within-subject residuals.}
 \item{resLevel}{\code{"marginal"} (default) or \code{"conditional"}. An optional character string specifying which residual should be used. For details see \code{\link{residuals.SMN}}.}
 \item{resType}{\code{"response"}, \code{"normalized"} (default), or \code{"modified"}. An optional character string specifying which type of residual should be used. For details see \code{\link{residuals.SMN}}.}
 \item{calcCI}{\code{TRUE} or \code{FALSE} (default). A logical value indicating if Monte Carlo confidence intervals should be computed for the conditionally independent model, which can be used for testing if the autocorrelations are zero.}
 \item{IClevel}{An optional numeric value in \eqn{(0,1)} indicating the confidence level that should be used in the Monte Carlo confidence intervals. Default is \code{0.95}.}
 \item{MCiter}{An optional discrete value indicating the number of Monte Carlo samples that should be used to compute the confidence intervals. Default is \code{300}.}
 \item{seed}{An optional value used to specify seeds inside the function. Default is to use a random seed.}
   }

\value{
A data frame with columns lag, ACF, and n.used representing, respectively, the lag between residuals within a pair, the corresponding empirical autocorrelation, and the number of pairs used. If \code{calcCI=TRUE}, the data frame has two extra columns containing the confidence intervals for the conditionally independent model. The returned value inherits from class acfresid.
}

\references{
Pinheiro, J. C. & Bates, D. M. (2000). Mixed-Effects Models in S and S-PLUS. Springer, New York, NY.

Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Statistics in Medicine. \url{https://doi.org/10.1002/sim.8870}.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smn.lmm}}, \code{\link{smsn.lmm}}, \code{\link{plot.acfresid}}
}
\examples{
fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
acf1<- acfresid(fm1)
acf1
plot(acf1)

\donttest{acfCI<- acfresid(fm1,calcCI=TRUE)
plot(acfCI)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
