\name{DrugsAnnot}
\alias{DrugsAnnot}
\docType{data}
\title{
Annotation of a set of drugs
}
\description{
  This data set is a matrix represents the drugs annotation that
  contains drugs' names, targets, mechanism of action , clinical phase,
  disease area, and indication.
}
\usage{data("DrugsAnnot")}
\format{
  A data frame with 6125 observations on the following 6 variables.
  \describe{
    \item{\code{pert_iname}}{a character vector}
    \item{\code{clinical_phase}}{a character vector}
    \item{\code{moa}}{a character vector}
    \item{\code{target}}{a character vector}
    \item{\code{disease_area}}{a character vector}
    \item{\code{indication}}{a character vector}
  }
}
\details{ 
L1000 repurposing drugs . For research use only. Do not use the
Repurposing Hub to make clinical treatment decisions. BROAD DOES NOT
GUARANTEE OR WARRANT THE ACCURACY OF THE DATA WITHIN THE REPURPOSING
HUB.
}
\source{
"The Drug Repurposing Hub, Broad Institute"	http://www.broadinstitute.org/repurposing
}
\references{
Corsello SM, et al. Nature Medicine. 2017 Apr 7;23(4):405-408. doi: 10.1038/nm.4306
}
\examples{
head(metaLINCS::DrugsAnnot)
## str(DrugsAnnot) 
}
\keyword{datasets}
