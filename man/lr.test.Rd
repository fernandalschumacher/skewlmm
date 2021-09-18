\name{lr.test}
\alias{lr.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Likelihood-ratio test for SMSN/SMN objects
%%  ~~function to do ... ~~
}
\description{
  It performs a likelihood-ratio test for two nested SMSN-LMM or SMN-LMM.
}
\usage{
lr.test(obj1, obj2, level = 0.05)
}
\arguments{
  \item{obj1, obj2}{smsn.lmm or smn.lmm objects containing the fitted models to be tested.}
  \item{level}{The significance level that should be used. If \code{quiet = TRUE}, this is ignored.}
}
\value{
  \item{statistic}{The test statistic value.}
  \item{p.value}{The p-value from the test.}
  \item{df}{The degrees of freedom used on the test.}
}
\references{
Schumacher, F. L., Lachos, V. H., and Matos, L. A. (2021). Scale mixture of skew-normal linear mixed models with within-subject serial dependence. Statistics in Medicine 40(7), 1790-1810.
}
\author{
Fernanda L. Schumacher, Larissa A. Matos and Victor H. Lachos
}


\seealso{
\code{\link{smsn.lmm}}, \code{\link{smn.lmm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
fm1<- smn.lmm(nlme::Orthodont, formFixed = distance ~ age+Sex, groupVar = "Subject")
fm2<- smsn.lmm(nlme::Orthodont, formFixed = distance ~ age+Sex, groupVar = "Subject",
               control = lmmControl(tol = .0001))
lr.test(fm1,fm2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}% use one of  RShowDoc("KEYWORDS")
