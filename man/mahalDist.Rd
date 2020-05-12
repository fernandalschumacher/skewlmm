\name{mahalDist}
\alias{mahalDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mahalanobis distance from a smn.lmm or smsn.lmm object
%%  ~~function to do ... ~~
}
\description{
    xx. }

\usage{
mahalDist(object,decomposed=FALSE,dataPlus=NULL)
}

\arguments{
 \item{object}{an object inheriting from class \code{SMN} or \code{SMSN}, representing
   a fitted scale mixture of (skew) normal linear mixed model.}
 \item{decomposed}{xxx.}
 \item{dataPlus}{xxx.}
   }

\value{
  A vector, if \code{decomposed=FALSE}, or a data frame, is \code{decomposed=TRUE}.
}

\references{
Schumacher, F.L., L.A. Matos, and V.H. Lachos (2020+). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Submitted.
Preprint available at \url{https://arxiv.org/abs/2002.01040}.

Zeller, C. B., Labra, F. V., Lachos, V. H. & Balakrishnan, N. (2010). Influence analyses of skew-normal/independent linear mixed models. Computational Statistics & Data Analysis, 54(5).
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smn.lmm}}, \code{\link{smsn.lmm}}, \code{\link{plot.mahalDist}}
}
\examples{
  fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
  mahalDist(fm1)
  plot(mahalDist(fm1),fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
