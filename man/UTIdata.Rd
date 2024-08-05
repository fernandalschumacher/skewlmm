\docType{data}
\name{UTIdata}
\alias{UTIdata}
\title{
Data set for Unstructured Treatment Interruption Study
}
\description{
Data set from a study of Unstructured Treatment Interruption in HIV-infected adolescents in four institutions in the USA. The main outcome is the HIV-1 RNA viral load, which is subject to censoring below the lower limit of detection of the assay (50 copies/mL). The censored observations are indicated by the variable RNAcens.
}
\usage{
data(UTIdata)
}
\format{
This data frame contains the following columns:
\describe{
    \item{\code{Patid}}{patient ID}
    \item{\code{Days.after.TI}}{days after treatment interruption}
    \item{\code{Fup}}{follow-up months}
    \item{\code{RNA}}{viral load RNA}
    \item{\code{RNAcens}}{censoring indicator for viral load}
}
}
\details{
This dataset was copied from the \code{lmec} package, which was discontinued from CRAN in May 2022.
}
\source{
Saitoh, A., Foca, M, et al. (2008), Clinical outcome in perinatally acquired HIV-infected children and adolescents after unstructured treatment interruption, Pediatrics,121, e513-e521.
}
\examples{
library(ggplot2)

data(UTIdata)
ggplot(UTIdata) + geom_line(aes(x=Fup, y=log10(RNA), group=Patid)) +
  geom_hline(yintercept=log10(50), color="red", linetype="dashed") +
  geom_hline(yintercept=log10(400), color="red", linetype="dashed") +
  labs(x="Time", y=bquote(log["10"]("RNA"))) + theme_bw()

# Proportion of censoring
prop.table(table(UTIdata$RNAcens))
}
\seealso{
\code{\link{smn.clmm}}
}
\keyword{datasets}
