\name{ranef}
\alias{ranef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract random effects from smsn.lmm and smn.lmm objects
%%  ~~function to do ... ~~
}
\description{
  It extracts random effects from smsn.lmm and smn.lmm objects.}

\usage{
ranef(object)
}

\arguments{
 \item{object}{an object inheriting from class \code{SMSN}, representing
   a fitted scale mixture skew-normal linear mixed model.}
   }

\value{
  Matrix of estimated random effects.
}

\references{
Lachos, V. H., P. Ghosh, and R. B. Arellano-Valle (2010). Likelihood based inference for skew-normal independent linear mixed models. Statistica Sinica 20, 303-322.

Schumacher, F.L., Lachos, V.H., and Matos, L.A. (2020+). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Submitted.
Preprint available at \url{https://arxiv.org/abs/2002.01040}.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smsn.lmm}}, \code{\link{smn.lmm}}, \code{\link{fitted.SMSN}}, \code{\link{fitted.SMN}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  fm1 <- smsn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject",tol = .0001)
  ranef(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
