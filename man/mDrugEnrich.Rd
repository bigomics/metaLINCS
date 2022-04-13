\name{mDrugEnrich}
\alias{mDrugEnrich}
\docType{data}
\title{
The drugs activity on 1001 genes
}
\description{
a large matrix represents the gene expression fold change of 1001 genes in the presence of 20220 different drugs concentrations and working time.
}
\usage{data("mDrugEnrich")}
\format{
  The format is:
 num [1:1001, 1:20220] -12.7 0.82 -1.48 -1.89 4.02 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:1001] "PSME1" "CISD1" "SPDEF" "ATF1" ...
  ..$ : chr [1:20220] "crizotinib_HEPG2_10um_24h" "avicin-d_JURKAT_10um_24h" "BRD-K35716340_WSUDLCL2_10um_6h" "trequinsin_HCC515_0.4um_24h" ...
}
\details{
The LINCS L1000 project has collected gene expression profiles for thousands of perturbagens at a variety of time points, doses, and cell lines. A full list of the chemical and genetic perturbations used can be found on the CLUE website along with their descriptions.
}
\source{
https://lincsproject.org/LINCS/tools/workflows/find-the-best-place-to-obtain-the-lincs-l1000-data
}
\references{
Duan Q et al. LINCS Canvas Browser: interactive web app to query, browse and interrogate LINCS L1000 gene expression signatures. Nucleic Acids Res. 2014 Jul;42 (Web Server issue):W449-60.
}
\examples{
data(mDrugEnrich)
## maybe str(mDrugEnrich) ; plot(mDrugEnrich) ...
}
\keyword{datasets}
