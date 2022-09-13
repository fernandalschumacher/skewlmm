\name{summary.SMNCens}
\alias{summary.SMNCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of a smnCens.lmm object
  %%  ~~function to do ... ~~
}
\description{
  \code{summary} method for class "SMNCens".
}

\usage{
  \method{summary}{SMNCens}(object, confint.level = 0.95, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{SMNCens},
  representing a fitted scale mixture normal linear mixed model with censored responses.}
  \item{confint.level}{Level of the approximate confidence intervals presented.}
  \item{\dots}{Additional arguments.}
}
\author{
  Fernanda L. Schumacher, Larissa A. Matos, Katherine L. Valeriano, and Victor H. Lachos
}
\seealso{
  \code{\link{smnCens.lmm}}, \code{\link{plot.SMN}}, \code{\link{residuals.SMNCens}}
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
