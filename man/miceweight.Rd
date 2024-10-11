\docType{data}
\name{miceweight}
\alias{miceweight}
\title{
Data set for clinical trial measuring mice weight
}
\description{
Sintetic longitudinal data set based on a clinical trial designed to test two diet treatments in comparison to a control group.
The weight of 52 mice (reported in grams) is measured weekly from baseline until week 10 of treatment.
}
\usage{
data(miceweight)
}
\format{
This data frame contains the following columns:
\describe{
    \item{\code{treat}}{Treatment received (treatment 1 (T1), treatment 2 (T2), or control (C))}
    \item{\code{mouseID}}{Mouse ID}
    \item{\code{week}}{Week since treatment started}
    \item{\code{weight}}{Mouse weight in grams}
}
}
\details{
Dietary supplementations were administered for groups T1 and T2 starting at week 5.

This dataset was created using minor perturbations from a real study for which the original data is not publicly available. Data features such as covariances and skewness were preserved.
}
\examples{
library(ggplot2)

data(miceweight)
ggplot(miceweight) + geom_line(aes(x=week, y=weight, group=mouseID)) +
  facet_wrap(~treat) + theme_bw()
}
\seealso{
\code{\link{smsn.lmm}}
}
\keyword{datasets}
