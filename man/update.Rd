\name{update}
\alias{update.SMN}\alias{update.SMSN}\alias{update.SMNCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Update for SMSN/SMN/SMNCens objects
%%  ~~function to do ... ~~
}
\description{
  It fits a SM(S)N-(C)LMM by updating a fitted object.
}
\usage{
\method{update}{SMN}(object, \dots, evaluate = TRUE)

\method{update}{SMSN}(object, \dots, evaluate = TRUE)

\method{update}{SMNCens}(object, \dots, evaluate = TRUE)
}
\arguments{
  \item{object}{A smsn.lmm, smn.lmm or smn.clmm object containing the fitted model to be updated.}
  \item{\dots}{Arguments to be changed.}
  \item{evaluate}{A logical value indicating if the new class should be evaluated. If \code{FALSE}, the call is returned.}
}
\value{
An object resulting from the \code{smsm.lmm()}, \code{smn.lmm()} or \code{smn.clmm()} function.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos, Victor H. Lachos and Katherine L. Valeriano
}

\seealso{
\code{\link{smsn.lmm}}, \code{\link{smn.lmm}}, \code{\link{smn.clmm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
fm1 = smn.lmm(nlme::Orthodont, formFixed=distance ~ age+Sex, groupVar="Subject")
fm2 = update(fm1, distr="t")
fm2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{robust}% __ONLY ONE__ keyword per line
