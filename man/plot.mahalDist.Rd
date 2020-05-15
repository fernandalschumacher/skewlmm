\name{plot.mahalDist}
\alias{plot.mahalDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Mahalanobis distance for a fitted smn.lmm or smsn.lmm
%%  ~~function to do ... ~~
}
\description{
Plot method for objects of class "mahalDist". For the total Mahalanobis distance, it gives a quantile for outlier detection, based on the Mahalanobis distance theoretical distribution.
}

\usage{
\method{plot}{mahalDist}(x,fitobject,type,level=.99,nlabels=3,\dots)
}

\arguments{
 \item{x}{An object inheriting from class \code{mahalDist}, representing
   the Mahalanobis distance from a fitted scale mixture of (skew) normal linear mixed model.}
  \item{fitobject}{An object inheriting from class \code{SMN} or \code{SMSN}, representing
   the fitted scale mixture of (skew) normal linear mixed model that was used for calculating the Mahalanobis distance.}
  \item{type}{Optional. Either \code{"total"} (default), for the standard Mahalanobis distance, \code{"error"}, for the error term of the decomposition, or \code{"b"} for the random effect term of the decomposition. For details see \code{\link{mahalDist}}.}
  \item{level}{An optional numeric value in \eqn{(0,1)} indicating the level of the quantile. It only has utility if \code{type="total"}. Default is \code{0.99}.}
  \item{nlabels}{Number of observations that should be labeled. Default is 3.}
  \item{\dots}{Additional arguments.}
   }

\value{
  A ggplot object, plotting the index versus the Mahalanobis distance, if all subject have the same number of observations; or plotting the number of observations per subject versus the Mahalanobis, otherwise.
}

\references{
Schumacher, F.L., V.H. Lachos, and L.A. Matos (2020+). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Submitted.
Preprint available at \url{https://arxiv.org/abs/2002.01040}.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{ggplot}}, \code{\link{mahalDist}}, \code{\link{smn.lmm}}, \code{\link{smsn.lmm}}
}
\examples{
  fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
  plot(mahalDist(fm1),fm1,nlabels=2)

  #the estimated quantile is stored at the attribute "info" of the plot object
  plotMD <- plot(mahalDist(fm1),fm1)
  attr(plotMD, "info")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
