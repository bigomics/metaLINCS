# SpaceLINCS: Stratified Perturbation Analysis by Connectivity Enrichment of LINCS L1000 signatures

SpaceLINCS calculates and visualizes the correlation between your
experimental gene expression profile with perturbations signatures
from the [LINCS L1000](https://lincsproject.org/LINCS/) drug
connectivity map. Summarizing the analysis with these perturbation
databases is difficult because they consist of more than a million of
profiles, corresponding to different cell lines and varying treatment
concentrations. SpaceLINCS attempts to efficiently calculate and
easily visualize the results by performing stratified enrichment tests
on the connectivity scores. In this way, mechanism-of-action or gene
targets are easily evident from the analysis.

## Installation

You can install the development version of SpaceLINCS from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("bigomics/SpaceLINCS")
```

## Example

This is a basic example which shows you how to use SpaceLINCS:

```{r}
load(system.file("data","mFC.rda", package = "SpaceLINCS"))
load(system.file("data","DrugsAnnot.rda", package = "SpaceLINCS"))
load(system.file("data","mDrugEnrich.rda", package = "SpaceLINCS"))

drugs <- computePeturbEnrichment(mFC, mDrugEnrich, DrugsAnnot, methods = c("GSEA", "cor"))
dsea <- getActiveDSEA(mDrugEnrich, DrugsAnnot, mFC)

### Get the mechanism of action data for each drug getMOA()
Moa <- getMOA(dsea)
DMoa <- Moa$drugClass
Dtarget <- Moa$geneTargets

### Plot the drugs connectivity using plotDrugConnectivity()
plotDrugConnectivity(dsea = dsea)

### Plot the mechanism of action using plotMOA()
plotMOA(dsea)

### Plot the drugs activity map using dseaPlotActmap()
dseaPlotActmap (dsea, drugs, method = "GSEA", contr=colnames(mFC)[1], nterms = 50, nfc=20)
```
