\name{weight_plot}
\alias{weight_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weight plot for smn.lmm or smsn.lmm object
%%  ~~function to do ... ~~
}
\description{
Estimated weights versus Mahalanobis distance plot
}

\usage{
weight_plot(object)
}

\arguments{
 \item{object}{An object inheriting from class \code{SMN} or \code{SMSN}, representing
   a fitted scale mixture of (skew) normal linear mixed model.}
   }

\value{
  A ggplot object.
}

\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\link{ggplot}, \link{smn.lmm}, \link{smsn.lmm}, \link{fitted.SMN}, \link{fitted.SMSN},
\link{residuals.SMN}, \link{residuals.SMSN}
}
\examples{
fm1 = smn.lmm(distance ~ age+Sex, data=nlme::Orthodont,
              groupVar="Subject", distr="t")
weight_plot(fm1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
