\name{plot.SMNclmm}
\alias{plot.SMNclmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a smn.clmm object
  %%  ~~function to do ... ~~
}
\description{
  Fitted values versus residuals plot. Censored values are imputed using their conditional expectation from the fitting algorithm.
}

\usage{
  \method{plot}{SMNclmm}(x, level = "conditional", useweight = TRUE,
                         alpha = 0.3, \dots)
}

\arguments{
  \item{x}{An object inheriting from class \code{SMNclmm}, representing
    a fitted scale mixture of normal censored linear mixed model.}
  \item{level}{Level of residual that should be used. For details see \code{\link{residuals.SMNclmm}}. Default is "conditional".}
  \item{useweight}{A logical value indicating if the estimated weights should be used as color in the plot.}
  \item{alpha}{Transparency parameter to be used (\code{0<alpha<1}). Meaningless if \code{useweight = TRUE}.}
  \item{\dots}{Additional arguments.}
}

\value{
  A \link[ggplot2:ggplot]{ggplot} object.
}

\author{
  Fernanda L. Schumacher, Larissa A. Matos, Victor H. Lachos and Katherine L. Valeriano
}


\seealso{
\link[ggplot2:ggplot]{ggplot}, \link{smn.clmm}, \link{fitted.SMNclmm}, \link{residuals.SMNclmm}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
