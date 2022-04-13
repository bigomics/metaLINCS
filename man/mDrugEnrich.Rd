\name{mDrugEnrich}
\alias{mDrugEnrich}
\docType{data}
\title{
The drugs activity on 1001 genes
}
\description{

  Large matrix representing the gene expression fold change of 1001
  genes in the presence of 20220 different drugs concentrations and
  working time.

}
\usage{
   SpaceLINCS::mDrugEnrich
}
\format{
  
  A gene expression matrix (log2FC) with 1001 rows and 20220 columns
  corresponding to 1001 genes and 20220 drug perturbation experiments.

}
\details{

  The LINCS L1000 project has collected gene expression profiles for
  thousands of perturbagens at a variety of time points, doses, and cell
  lines. A full list of the chemical and genetic perturbations used can
  be found on the CLUE website along with their descriptions.

}

\source{

  https://lincsproject.org/LINCS/tools/workflows/find-the-best-place-to-obtain-the-lincs-l1000-data

}
\references{

  Duan Q et al. LINCS Canvas Browser: interactive web app to query,
  browse and interrogate LINCS L1000 gene expression signatures. Nucleic
  Acids Res. 2014 Jul;42 (Web Server issue):W449-60.

}
\examples{
head(SpaceLINCS::mDrugEnrich)
## maybe str(mDrugEnrich) ; plot(mDrugEnrich) ...
}
\keyword{datasets}
