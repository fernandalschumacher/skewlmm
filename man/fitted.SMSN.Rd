\name{fitted.SMSN}
\alias{fitted.SMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract smsn.lmm fitted values
%%  ~~function to do ... ~~
}
\description{
 The fitted values are obtained by adding together the
  population fitted values (based only on the fixed effects estimates)
  and the estimated contributions of the random effects to the fitted
  values at grouping levels. }

\usage{
\method{fitted}{SMSN}(object, \dots)
}

\arguments{
 \item{object}{an object inheriting from class \code{SMSN}, representing
   a fitted scale mixture skew-normal linear mixed model.}
  \item{\dots}{Additional arguments.}
   }

\value{
  Vector of fitted values with length equal to \code{nrow(data)}.
}

\references{
Lachos, V. H., P. Ghosh, and R. B. Arellano-Valle (2010). Likelihood based inference for skew-normal independent linear mixed models. Statistica Sinica 20, 303-322.

Schumacher, F.L., L.A. Matos, and V.H. Lachos (2020+). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Submitted.
Preprint available at \url{https://arxiv.org/abs/2002.01040}.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smsn.lmm}}, \code{\link{predict.SMSN}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
    fm1 <- smsn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
    fitted(fm1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
