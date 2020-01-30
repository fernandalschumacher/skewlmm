\name{predict.SMSN}
\alias{predict.SMSN}

\title{Prediction of future observations from an smsn.lmm object
}
\description{
  Predicted values are obtained through conditional expectation. For details, see Schumacher, Matos and Lachos (2020).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{predict}{SMSN}(object, newData,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{an object inheriting from class \code{SMSN}, representing
   a fitted scale mixture skew-normal linear mixed model.}
 \item{newData}{a data frame for which response variable should be predicted, including covariates, groupVar and possibly timeVar.}
  \item{\dots}{Additional arguments.}
   }

\value{
  A data frame with covariates, groupVar and ypred, where ypred contains the predicted values from the response variable.
}

\references{
Schumacher, F.L., Matos, L.A., and Lachos, V.H. (2020+) "Scale mixture of skew-normal linear mixed models with within-subject serial dependence". Submitted.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smsn.lmm}}, \code{\link{fitted.SMSN}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
    dat1 <- nlme::Orthodont
    fm1 <- smsn.lmm(distance ~ age + Sex, data = subset(dat1,age<14), groupVar="Subject")
    predict(fm1,subset(dat1,age==14))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
