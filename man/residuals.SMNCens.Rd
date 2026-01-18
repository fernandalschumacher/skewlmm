\name{residuals.SMNclmm}
\alias{residuals.SMNclmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract model residuals from smn.clmm objects
  %%  ~~function to do ... ~~
}
\description{
  The conditional residuals are obtained by subtracting the fitted values from the response vector, while the marginal residuals are obtained by subtracting only the fixed effects from the response vector.

  Censored values are imputed using their conditional expectation from the fitting algorithm.}

\usage{
  \method{residuals}{SMNclmm}(object, level = "conditional", \dots)
}

\arguments{
  \item{object}{An object inheriting from class \code{SMNclmm}, representing
    a fitted scale mixture of normal censored linear mixed model.}
  \item{level}{Either \code{"conditional"}, for obtaining conditional residuals, or \code{"marginal"}, for marginal residuals.}
  \item{\dots}{Additional arguments.}
}

\value{
  Vector with the residuals of length equal to \code{nrow(data)}.
}
\note{
The residuals are computed after imputing the censored observations.}
\author{
  Fernanda L. Schumacher, Larissa A. Matos, Victor H. Lachos and Katherine L. Valeriano
}

\seealso{
  \code{\link{smn.clmm}}, \code{\link{mahalDistCens}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
