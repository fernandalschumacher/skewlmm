\name{plot}
\alias{plot.SMN}\alias{plot.SMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a smn.lmm or smsn.lmm object
%%  ~~function to do ... ~~
}
\description{
Fitted versus conditional standardized residuals plot.
}

\usage{
\method{plot}{SMN}(x, alpha=.3,\dots)

\method{plot}{SMSN}(x, alpha=.3,\dots)
}

\arguments{
 \item{x}{An object inheriting from class \code{SMN} or \code{SMSN}, representing
   a fitted scale mixture of (skew) normal linear mixed model.}
  \item{alpha}{Transparency parameter to be used (\code{0<alpha<1}).}
  \item{\dots}{Additional arguments.}
   }

\value{
  A ggplot object.
}

\references{
Schumacher, F.L., L.A. Matos, and V.H. Lachos (2020+). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Submitted.
Preprint available at \url{https://arxiv.org/abs/2002.01040}.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{ggplot}}, \code{\link{smn.lmm}}, \code{\link{smsn.lmm}}, \code{\link{fitted.SMN}},
\code{\link{fitted.SMSN}},\code{\link{residuals.SMN}}, \code{\link{residuals.SMSN}}
}
\examples{
  fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
  plot(fm1)

  \donttest{library(ggplot2)
  plot(fm1) + ggtitle("Model fm1")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
