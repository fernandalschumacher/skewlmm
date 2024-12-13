\name{plot.acfresid}
\alias{plot.acfresid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot ACF for smn.lmm or smsn.lmm residuals
%%  ~~function to do ... ~~
}
\description{
Plot method for objects of class "acfresid".
}

\usage{
\method{plot}{acfresid}(x,\dots)
}

\arguments{
 \item{x}{An object inheriting from class \code{acfresid}, representing the empirical autocorrelation function from the residuals of a scale mixture of (skew) normal linear mixed model.}
  \item{\dots}{Additional arguments.}
   }

\value{
  A \link[ggplot2:ggplot]{ggplot} object.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\link[ggplot2:ggplot]{ggplot}, \link{acfresid}, \link{smn.lmm}, \link{smsn.lmm}, \link{residuals.SMN},
\link{residuals.SMSN}
}
\examples{
fm1 = smn.lmm(distance ~ age+Sex, data=nlme::Orthodont, groupVar="Subject")
plot(acfresid(fm1))

\donttest{acfCI = acfresid(fm1, calcCI=TRUE)
plot(acfCI)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
