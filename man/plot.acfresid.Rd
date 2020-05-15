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
  A ggplot object.
}

\references{
Schumacher, F.L., V.H. Lachos, and L.A. Matos (2020+). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Submitted.
Preprint available at \url{https://arxiv.org/abs/2002.01040}.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{ggplot}}, \code{\link{acfresid}}, \code{\link{smn.lmm}}, \code{\link{smsn.lmm}}, \code{\link{residuals.SMN}}, \code{\link{residuals.SMSN}}
}
\examples{
  fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
  plot(acfresid(fm1))

  \donttest{acfCI <- acfresid(fm1,calcCI=TRUE)
  plot(acfCI)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
