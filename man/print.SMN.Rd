\name{print.SMN}
\alias{print.SMN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print a smn.lmm object
%%  ~~function to do ... ~~
}
\description{
Print a smn.lmm object.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{print}{SMN}(x, \dots)
}

\arguments{
 \item{x}{an object inheriting from class \code{SMN}, representing
   a fitted scale mixture normal linear mixed model.}
  \item{\dots}{Additional print arguments.}
   }

\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smn.lmm}}, \code{\link{summary.SMN}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar="Subject")
fm1
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
