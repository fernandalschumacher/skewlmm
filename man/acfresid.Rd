\name{acfresid}
\alias{acfresid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Autocorrelation function for smn.lmm or smsn.lmm residuals
%%  ~~function to do ... ~~
}
\description{
    xx. }

\usage{
acfresid(object,maxLag,resLevel="marginal",resType="normalized",
                     calcCI=FALSE,IClevel,MCiter,seed)
}

\arguments{
 \item{object}{an object inheriting from class \code{SMN} or \code{SMSN}, representing
   a fitted scale mixture of (skew) normal linear mixed model.}
 \item{maxLag}{xxx.}
 \item{resLevel}{xxx.}
 \item{resType}{xxx.}
 \item{calcCI}{xxx.}
 \item{IClevel}{xxx.}
 \item{MCiter}{xxx.}
 \item{seed}{xxx.}
   }

\value{
  A data frame...
}

\references{
Pinheiro, J. C. & Bates, D. M. (2000). Mixed-Effects Models in S and S-PLUS. Springer, New York, NY.

Schumacher, F.L., L.A. Matos, and V.H. Lachos (2020+). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Submitted.
Preprint available at \url{https://arxiv.org/abs/2002.01040}.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smn.lmm}}, \code{\link{smsn.lmm}}, \code{\link{plot.ACF}}
}
\examples{
  fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
  acf1<- acfresid(fm1)
  acf1
  plot(acf1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
