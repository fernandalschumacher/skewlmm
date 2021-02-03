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
healy.plot(object,dataPlus=NULL,dotsize=0.4,...)
}

\arguments{
 \item{object}{An object inheriting from class \code{SMN} or \code{SMSN}, representing
   a fitted scale mixture of (skew) normal linear mixed model.}
  \item{dataPlus}{Optional. Expanded dataset that should be used instead the one used for fitting. This is necessary for unbalanced datasets, since Haley's plot requires all subject to have the same number of observations.}
  \item{dotsize}{Optional. Dotsize used in ggplot.}
  \item{\dots}{Additional arguments.}
   }

\value{
  A ggplot object.
}

\details{It constructs a Healy-type plot (Healy, 1968) by plotting the nominal probability values \eqn{1/n,2/n,...,n/n} against the theoretical cumulative probabilities of the ordered observed Mahalanobis distances. If the fitted model is appropriate, the plot
should resemble a straight line through the origin with unit slope.
}

\references{
Healy, M. J. R. (1968). Multivariate normal plotting. Journal of the Royal Statistical Society: Series C (Applied Statistics), 17(2), 157-161.

Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Statistics in Medicine. \url{https://doi.org/10.1002/sim.8870}.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}

\seealso{
\code{\link{ggplot}}, \code{\link{smn.lmm}}, \code{\link{smsn.lmm}},
\code{\link{mahalDist}}, \code{\link{acfresid}}
}
\examples{
fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
healy.plot(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
