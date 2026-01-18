\name{predict.SMNclmm}
\alias{predict.SMNclmm}

\title{Prediction of future observations from an smn.clmm object
}
\description{
  Predicted values are obtained through conditional expectation. For details, see Schumacher, Lachos, and Matos (2021).
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  \method{predict}{SMNclmm}(object, newData, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object inheriting from class \code{SMNclmm}, representing
    a fitted scale mixture normal linear mixed model.}
  \item{newData}{A data frame for which response variable should be predicted, including covariates, groupVar and possibly timeVar. If missing or \code{NULL}, fitted values are returned.}
  \item{\dots}{Additional arguments.}
}

\value{
  A data frame with covariates, groupVar and ypred, where ypred contains the predicted values from the response variable.
}

\references{
  Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Statistics in Medicine 40(7), 1790-1810.
}
\author{
  Fernanda L. Schumacher, Larissa A. Matos, Victor H. Lachos and Katherine L. Valeriano
}


\seealso{
  \code{\link{smn.clmm}}, \code{\link{fitted.SMNclmm}}
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
