\name{errorVar}
\alias{errorVar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Return error variance associated with \code{times}
%%  ~~function to do ... ~~
}
\description{
  It returns a variance matrix associated with the error term at time \code{times}.
  Can be applied to a smsn.lmm object or to a specific dependence structure with chosen parameter values.
}
\usage{
errorVar(times, object = NULL, sigma2 = NULL, depStruct = NULL,
    phi = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{A vector containing the times for which the matrix should be calculated.}
  \item{object}{A smsn.lmm object for which the variance should be extracted.}
  \item{sigma2}{Common variance parameter, such that \eqn{\Sigma = sigma2*R}. Only evaluated if object \code{is.null(object)}.}
  \item{depStruct}{Dependence structure. \code{"CI"} for condicional independence, \code{"ARp"} for AR(p) -- p is \code{length(phi)}--, \code{"CS"} for compound symmetry, \code{"DEC"} for DEC, and \code{"CAR1"} for  continuous-time AR(1). Only evaluated if object \code{is.null(object)}.}
  \item{phi}{Parameter vector indexing the dependence structure. Only evaluated if object \code{is.null(object)}.}
}
\value{
  Matrix of dimension \code{length(times)}.
}
\references{
Schumacher, F.L., Matos, L.A., and Lachos, V.H. (2020+) "Scale mixture of skew-normal linear mixed models with within-subject serial dependence". Submitted.
Preprint available at \url{https://arxiv.org/abs/2002.01040}.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smsn.lmm}}, \code{\link{summary.SMSN}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\donttest{
  fm1 <- smsn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
  errorVar(times=1:4,fm1)
  #}
  errorVar(times=1:5,sigma2 = 1, depStruct = "ARp",phi = .5)
  errorVar(times=1:5,sigma2 = 1, depStruct = "DEC",phi = c(.5,.8))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
