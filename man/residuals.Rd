\name{residuals}
\alias{residuals.SMN}\alias{residuals.SMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract model residuals from smn.lmm and smsn.lmm objects
%%  ~~function to do ... ~~
}
\description{
 The conditional residuals are obtained by subtracting the fitted values from the response vector, while the marginal residuals are obtained by subtracting only the fixed effects from the response vector.}

\usage{
\method{residuals}{SMN}(object,level="conditional",type="response",\dots)

\method{residuals}{SMSN}(object,level="conditional",type="response",\dots)
}

\arguments{
 \item{object}{An object inheriting from class \code{SMN} or \code{SMSN}, representing
   a fitted scale mixture of (skew) normal linear mixed model.}
  \item{level}{Either \code{"conditional"}, for obtaining conditional residuals, or \code{"marginal"}, for marginal residuals.}
  \item{type}{An optional character string specifying the type of residuals to be used. If \code{"response"}, as by default, the "raw" residuals (observed - fitted) are used; if \code{"normalized"}, the standardized residuals (residuals pre-multiplied by the inverse square-root of the estimated variance matrix) are used; else, if \code{"modified"}, modified residuals (residuals pre-multiplied by the inverse square-root of the estimated scale matrix) are used.}
  \item{\dots}{Additional arguments.}
   }

\value{
  Vector with the residuals of length equal to \code{nrow(data)}.
}

\details{
Modified residuals are useful when the variance is not finite, such as when \eqn{\nu \le 2} for t or ST distributions, or when \eqn{\nu \le 1} for SL or SSL distributions.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smn.lmm}}, \code{\link{smsn.lmm}}, \code{\link{acfresid}},
\code{\link{mahalDist}}, \code{\link{healy.plot}}
}
\examples{fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
residuals(fm1)
plot(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
