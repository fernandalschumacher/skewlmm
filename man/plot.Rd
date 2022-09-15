\name{plot}
\alias{plot.SMN}\alias{plot.SMSN}\alias{plot.SMNCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a smn.lmm, smsn.lmm or smnCens.lmm object
%%  ~~function to do ... ~~
}
\description{
Fitted values versus residuals plot.
}

\usage{
\method{plot}{SMN}(x, type = "response", level = "conditional",
     useweight = TRUE, alpha = 0.3, \dots)

\method{plot}{SMSN}(x, type = "response", level = "conditional",
     useweight = TRUE, alpha = 0.3, \dots)

\method{plot}{SMNCens}(x, type = "response", level = "conditional",
     useweight = TRUE, alpha = 0.3, \dots)
}

\arguments{
 \item{x}{An object inheriting from class \code{SMN}, \code{SMSN} or \code{SMNCens}, representing
   a fitted scale mixture of (skew) normal linear mixed model.}
  \item{type}{Type of residual that should be used. For details see \code{\link{residuals.SMN}}. Default is "response", indicating raw residuals.}
  \item{level}{Level of residual that should be used. For details see \code{\link{residuals.SMN}}. Default is "conditional".}
  \item{useweight}{A logical value indicating if the estimated weights should be used as color in the plot.}
  \item{alpha}{Transparency parameter to be used (\code{0<alpha<1}). Meaningless if \code{useweight = TRUE}.}
  \item{\dots}{Additional arguments.}
   }

\value{
  A ggplot object.
}

\author{
Fernanda L. Schumacher, Larissa A. Matos, Victor H. Lachos and Katherine L. Valeriano
}


\seealso{
\code{\link{ggplot}}, \code{\link{smn.lmm}}, \code{\link{smsn.lmm}}, \code{\link{smnCens.lmm}}, \code{\link{fitted.SMN}},
\code{\link{fitted.SMSN}}, \code{\link{fitted.SMNCens}}, \code{\link{residuals.SMN}}, \code{\link{residuals.SMSN}}, \code{\link{residuals.SMNCens}}
}
\examples{
fm1 = smn.lmm(distance ~ age+Sex, data=nlme::Orthodont,
              groupVar="Subject", distr="t")
plot(fm1)
plot(fm1, useweight=FALSE)

\donttest{library(ggplot2)
plot(fm1) + ggtitle("t-LMM for orthodont data")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
