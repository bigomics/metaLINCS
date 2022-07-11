## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  install.packages("remotes")
#  remotes::install_github("bigomics/metaLINCS")

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  ## source code within computeConnectivityEnrichment function, do not
#  ## run it, it is for method illustration not part of the workflow
#  
#  ##  first calculation of rank correlation
#  gg   <- intersect(rownames(mDrugEnrich), rownames(mFC))
#  rnk1 <- apply(mDrugEnrich[gg, , drop = FALSE], 2, rank, na.last = "keep")
#  rnk2 <- apply(mFC[gg, , drop = FALSE], 2, rank, na.last = "keep")
#  R1   <- stats::cor(rnk1, rnk2, use = "pairwise")
#  
#  # Then, we calculate the perturbation enrichment using GSEA
#  xdrugs <- gsub("[@_].*$", "", colnames(mDrugEnrich))
#  meta.gmt <- tapply(colnames(mDrugEnrich), xdrugs, list)
#  res <- list()
#  for (i in 1:ncol(R1)) {
#     suppressWarnings(res[[i]] <- fgsea::fgseaSimple(meta.gmt, stats = R1[, i], nperm = 1000))
#  }

## -----------------------------------------------------------------------------
library(metaLINCS)

## -----------------------------------------------------------------------------
head(mFC)

## -----------------------------------------------------------------------------
head(DrugsAnnot)

## -----------------------------------------------------------------------------
dim(mDrugEnrich)
head(mDrugEnrich)[,1:2]

## ----message=FALSE, warning=FALSE, paged.print=FALSE--------------------------
res <- computeConnectivityEnrichment(mFC, nprune=0)
names(res)

## -----------------------------------------------------------------------------
res1 <- selectResult(res,1)
head(res1)

## -----------------------------------------------------------------------------
moa <- computeMoaEnrichment(res) 
names(moa)

## -----------------------------------------------------------------------------
## Get the mechanism of action results for the first contrast
head(moa[[1]]$drugClass)
head(moa[[1]]$targetGene)

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  ## Plot the drugs connectivity using plotDrugConnectivity()
#  plotDrugConnectivity(res, contr="Resistant.vs.Sensitive")

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  ## Plot the mechanism of action using plotMOA()
#  plotMOA(moa, contr="WithaferinA.vs.Untreated", type="drugClass", ntop=20)

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  ## Plot the mechanism of action using plotMOA()
#  plotMOA(moa, contr="WithaferinA.vs.Untreated", type="targetGene", ntop=20)

## ---- eval=FALSE, echo=TRUE---------------------------------------------------
#  ## Plot the drugs activity map using dseaPlotActmap()
#  plotActivationMap(res, nterms = 60, nfc=20)

