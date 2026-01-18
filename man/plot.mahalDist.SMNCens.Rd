\name{plot.mahalDistCens}
\alias{plot.mahalDistCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Mahalanobis distance for a fitted smn.clmm
  %%  ~~function to do ... ~~
}
\description{
  Plot method for objects of class "mahalDistCens". It also gives a quantile for outlier detection, based on the Mahalanobis distance theoretical distribution.
}

\usage{
  \method{plot}{mahalDistCens}(x, fitobject, level = 0.99, nlabels = 3, \dots)
}

\arguments{
  \item{x}{An object inheriting from class \code{mahalDistCens}, representing
    the Mahalanobis distance from a fitted scale mixture of normal censored linear mixed model.}
  \item{fitobject}{Optional. An object inheriting from class \code{SMNclmm}, representing
    the fitted scale mixture of normal linear mixed model that was used for calculating the Mahalanobis distance.}
  \item{level}{An optional numeric value in \eqn{(0,1)} indicating the level of the quantile. Default is \code{0.99}.}
  \item{nlabels}{Number of observations that should be labeled. Default is 3.}
  \item{\dots}{Additional arguments.}
}

\value{
  A \link[ggplot2:ggplot]{ggplot} object, plotting the index versus the Mahalanobis distance, if all subject have the same number of observations; or plotting the number of observations per subject versus the Mahalanobis, otherwise.
}

\author{
  Fernanda L. Schumacher, Larissa A. Matos, Victor H. Lachos and Katherine L. Valeriano
}


\seealso{
\link[ggplot2:ggplot]{ggplot}, \link{mahalDistCens}, \link{smn.clmm}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
