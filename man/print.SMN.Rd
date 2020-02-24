\name{print.SMN}
\alias{print.SMN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print a smn.lmm object
%%  ~~function to do ... ~~
}
\description{
Print a smn.lmm object.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{print}{SMN}(x, \dots)
}

\arguments{
 \item{x}{an object inheriting from class \code{SMN}, representing
   a fitted scale mixture normal linear mixed model.}
  \item{\dots}{Additional print arguments.}
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
\code{\link{smn.lmm}}, \code{\link{summary.SMN}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
  fm1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
