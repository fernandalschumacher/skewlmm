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

\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smsn.lmm}}, \code{\link{predict.SMSN}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
fm1 <- smsn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar = "Subject",
       control = lmmControl(tol = .0001))
fitted(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
