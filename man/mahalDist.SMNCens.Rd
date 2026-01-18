\name{mahalDistCens}
\alias{mahalDistCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mahalanobis distance from a smn.clmm object
  %%  ~~function to do ... ~~
}
\description{
Returns the squared Mahalanobis distance from a fitted SMN-CLMM. Censored values are imputed using their conditional expectation from the fitting algorithm.}

\usage{
  mahalDistCens(object)
}

\arguments{
  \item{object}{An object inheriting from class \code{SMNclmm}, representing
    a fitted scale mixture of normal censored linear mixed model.}
}

\value{
  An object of class \code{mahalDistCens} containing the Mahalanobis distance.
}

\references{
Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Statistics in Medicine 40(7), 1790-1810.

Zeller, C. B., Labra, F. V., Lachos, V. H. & Balakrishnan, N. (2010). Influence analyses of skew-normal/independent linear mixed models. Computational Statistics & Data Analysis, 54(5).
}
\author{
Fernanda L. Schumacher, Larissa A. Matos, Victor H. Lachos and Katherine L. Valeriano
}


\seealso{
\code{\link{smn.clmm}}, \code{\link{plot.mahalDistCens}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
