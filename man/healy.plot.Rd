\name{healy.plot}
\alias{healy.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Healy-type plot from a smn.lmm or smsn.lmm object
%%  ~~function to do ... ~~
}
\description{
It creates a Healy-type plot from a smn.lmm or smsn.lmm object, for goodness-of-fit assessment.

}

\usage{
healy.plot(object, dataPlus = NULL, dotsize = 0.4, calcCI = FALSE,
           levelCI, MCiter, seed, ...)
}

\arguments{
 \item{object}{An object inheriting from class \code{SMN} or \code{SMSN}, representing
   a fitted scale mixture of (skew) normal linear mixed model.}
  \item{dataPlus}{Optional. Expanded dataset that should be used instead the one used for fitting. This is necessary for unbalanced datasets, since Haley's plot requires all subject to have the same number of observations.}
  \item{dotsize}{Optional. Dotsize used in ggplot.}
  \item{calcCI}{\code{TRUE} or \code{FALSE} (default). A logical value indicating if Monte Carlo confidence intervals should be computed for the conditionally independent model, which can be used for testing if the autocorrelations are zero.}
 \item{levelCI}{An optional numeric value in \eqn{(0,1)} indicating the confidence level that should be used in the Monte Carlo confidence intervals. Default is \code{0.95}.}
 \item{MCiter}{An optional discrete value indicating the number of Monte Carlo samples that should be used to compute the confidence intervals. Default is \code{300}.}
 \item{seed}{An optional value used to specify seeds inside the function. Default is to use a random seed.}
 \item{\dots}{Additional arguments.}
   }

\value{
  A \link[ggplot2:ggplot]{ggplot} object.
}

\details{It constructs a Healy-type plot (Healy, 1968) by plotting the nominal probability values \eqn{1/n,2/n,...,n/n} against the theoretical cumulative probabilities of the ordered observed Mahalanobis distances. If the fitted model is appropriate, the plot
should resemble a straight line through the origin with unit slope.
 If \code{calcCI=TRUE}, the plot presents two dashed blue lines containing approximated confidence intervals by considering that the fitted model is correct.
}

\references{
Healy, M. J. R. (1968). Multivariate normal plotting. Journal of the Royal Statistical Society: Series C (Applied Statistics), 17(2), 157-161.

Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Statistics in Medicine 40(7), 1790-1810.}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}

\seealso{
\link[ggplot2:ggplot]{ggplot}, \link{smn.lmm}, \link{smsn.lmm},
\link{mahalDist}, \link{acfresid}
}
\examples{
fm1 = smn.lmm(distance ~ age+Sex, data=nlme::Orthodont, groupVar="Subject")
healy.plot(fm1)

## computing simulated bands
\donttest{healy.plot(fm1, calcCI=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
