\name{mahalDist}
\alias{mahalDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mahalanobis distance from a smn.lmm or smsn.lmm object
%%  ~~function to do ... ~~
}
\description{
    Returns the squared Mahalanobis distance from a fitted SMN-LMM or SMSN-LMM.}

\usage{
mahalDist(object,decomposed=FALSE,dataPlus=NULL)
}

\arguments{
 \item{object}{An object inheriting from class \code{SMN} or \code{SMSN}, representing
   a fitted scale mixture of (skew) normal linear mixed model.}
 \item{decomposed}{Logical. If \code{TRUE}, the Mahalanobis distance is decomposed in an error term and a random effect term. Default is \code{FALSE}.}
 \item{dataPlus}{Optional. Expanded dataset that should be used instead the one used for fitting, useful for using Healy's plot with missing data.}
   }

\value{
  A vector containing the Mahalanobis distance, if \code{decomposed=FALSE}, or a data frame containing the Mahalanobis distance and its decomposition in error term and random effect (b) term, if \code{decomposed=TRUE}.
}

\references{
Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Statistics in Medicine. \doi{10.1002/sim.8870}.

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
