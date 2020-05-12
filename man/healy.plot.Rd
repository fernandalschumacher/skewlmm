\name{healy.plot}
\alias{healy.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates a Healy-type plot from a smn.lmm or smsn.lmm object
%%  ~~function to do ... ~~
}
\description{
xxx.
}

\usage{
healy.plot(object,dataPlus=NULL,dotsize=0.4,...)
}

\arguments{
 \item{object}{An object inheriting from class \code{SMN} or \code{SMSN}, representing
   a fitted scale mixture of (skew) normal linear mixed model.}
  \item{dataPlus}{xxx.}
  \item{dotsize}{xxx.}
  \item{\dots}{Additional arguments.}
   }

\value{
  A ggplot object.
}

\references{
Healy, M. J. R. (1968). Multivariate normal plotting. Journal of the Royal Statistical Society: Series C (Applied Statistics), 17(2), 157-161.

Schumacher, F.L., L.A. Matos, and V.H. Lachos (2020+). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Submitted.
Preprint available at \url{https://arxiv.org/abs/2002.01040}.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}

\seealso{
\code{\link{ggplot}}, \code{\link{smn.lmm}}, \code{\link{smsn.lmm}},
\code{\link{mahalDist}}, \code{\link{acfresid}}
}
\examples{
  fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
  healy.plot(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
