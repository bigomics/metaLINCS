\name{mFC}
\alias{mFC}
\docType{data}
\title{
Fold change matrix of differential gene expresion from RNAseq data of Multiple myeloma
}
\description{There are 6 contrasts :
1- glucocorticoid sensitive vs glucocorticoid resistant ("Resistant.vs.Sensitive").2- Treated with Withaferin-A vs untreated ("WithaferinA.vs.Untreated"). 3- untreated glucocorticoid resistant cells vs untreated glucocorticoid sensitive cells ("UTR.vs.UTS"). 4-  Treated with Withaferin-A glucocorticoid resistant cells vs Treated with Withaferin-A glucocorticoid sensitive cells ("WAR.vs.WAS"). 5- Treated Resistant vs Untreated Resistant ("WAR.vs.UTR"). 6- Treated sensitive vs Untreated sensitive ("WAS.vs.UTS") }
\usage{data("mFC")}
\format{
  A data frame with 5332 observations on the following 6 variables.
  \describe{
    \item{\code{Resistant.vs.Sensitive}}{a numeric vector}
    \item{\code{WithaferinA.vs.Untreated}}{a numeric vector}
    \item{\code{UTR.vs.UTS}}{a numeric vector}
    \item{\code{WAR.vs.WAS}}{a numeric vector}
    \item{\code{WAR.vs.UTR}}{a numeric vector}
    \item{\code{WAS.vs.UTS}}{a numeric vector}
  }
}
\details{
The RNAseq data from (Logie et al 2021) study. They compared the therapeutic efficacy of the phytochemical kinase inhibitor withaferin A with the clinically approved BTK inhibitor ibrutinib to target hyperactivated tyrosine kinase signaling in glucocorticoid-resistant multiple myeloma cells. The results demonstrate that withaferin-A induced cell death of glucocorticoid-resistant MM1R cells involves covalent cysteine targeting of multiple Hinge-6 domain type tyrosine kinases of the kinase cysteinome classification, including BTK.
}
\source{
GEO accession number GSE162475
}
\references{
Logie E, Chirumamilla CS, Perez-Novo C, Shaw P et al. Covalent Cysteine Targeting of Bruton's Tyrosine Kinase (BTK) Family by Withaferin-A Reduces Survival of Glucocorticoid-Resistant Multiple Myeloma MM1 Cell. Cancers (Basel) 2021 Mar 31;13(7). PMID: 33807411
}
\examples{
data(mFC)
## maybe str(mFC) ; plot(mFC) ...
}
\keyword{datasets}
