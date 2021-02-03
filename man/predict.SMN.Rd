\name{predict.SMN}
\alias{predict.SMN}

\title{Prediction of future observations from an smn.lmm object
}
\description{
  Predicted values are obtained through conditional expectation. For details, see Schumacher, Lachos, and Matos (2020).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{predict}{SMN}(object, newData,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{an object inheriting from class \code{SMN}, representing
   a fitted scale mixture normal linear mixed model.}
 \item{newData}{a data frame for which response variable should be predicted, including covariates, groupVar and possibly timeVar.}
  \item{\dots}{Additional arguments.}
   }

\value{
  A data frame with covariates, groupVar and ypred, where ypred contains the predicted values from the response variable.
}

\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smn.lmm}}, \code{\link{fitted.SMN}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dat1 <- nlme::Orthodont
fm1 <- smn.lmm(distance ~ age + Sex, data = subset(dat1,age<14), groupVar="Subject")
predict(fm1,subset(dat1,age==14))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
