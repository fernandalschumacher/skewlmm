\name{fitted.SMN}
\alias{fitted.SMN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract smn.lmm fitted values
%%  ~~function to do ... ~~
}
\description{
 The fitted values are obtained by adding together the
  population fitted values (based only on the fixed effects estimates)
  and the estimated contributions of the random effects to the fitted
  values at grouping levels. }


\usage{
\method{fitted}{SMN}(object, \dots)
}

\arguments{
 \item{object}{an object inheriting from class \code{SMN}, representing
   a fitted scale mixture normal linear mixed model.}
  \item{\dots}{Additional arguments.}
   }

\value{
  Vector of fitted values with length equal to \code{nrow(data)}.
}

\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smn.lmm}}, \code{\link{predict.SMN}}
}
\examples{
fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
fitted(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
