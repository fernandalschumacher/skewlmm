\name{print.SMSN}
\alias{print.SMSN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print a smsn.lmm object
%%  ~~function to do ... ~~
}
\description{
Print a smsn.lmm object.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{print}{SMSN}(x, \dots)
}

\arguments{
 \item{x}{an object inheriting from class \code{SMSN}, representing
   a fitted scale mixture skew-normal linear mixed model.}
  \item{\dots}{Additional print arguments.}
   }

\references{
Schumacher, F.L., Matos, L.A., and Lachos, V.H. (2020+) "Scale mixture of skew-normal linear mixed models with within-subject serial dependence". Submitted.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smsn.lmm}}, \code{\link{summary.SMSN}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
    fm1 <- smsn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
    fm1
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
