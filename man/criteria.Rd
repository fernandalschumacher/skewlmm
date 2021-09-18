\name{criteria}
\alias{criteria}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extracts criteria for model comparison of SMSN/SMN objects
%%  ~~function to do ... ~~
}
\description{
  It extracts criteria for model comparison of sereval SMSN-LMM and/or SMN-LMM.
}
\usage{
criteria(lobjects)
}
\arguments{
  \item{lobjects}{A list containg the smsn.lmm or smn.lmm objects to be compared.}
}
\value{
A \code{data.frame} containing for each model the maximum log-likelihood value, the number of parameters, the AIC, and the BIC.
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
fm_norm <- smn.lmm(nlme::Orthodont, formFixed = distance ~ age+Sex, groupVar = "Subject")
fm_t <- update(fm_norm, distr="t")
fm_sn<- smsn.lmm(nlme::Orthodont, formFixed = distance ~ age+Sex, groupVar = "Subject")
criteria(list(fm_norm = fm_norm,fm_t = fm_t,fm_sn = fm_sn))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
