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
library(SpaceLINCS)

## First we compute the connectivity enrichment    
res <- computeConnectivityEnrichment(mFC, nprune=0)
names(res)

## Now compute the MOA enrichment
moa <- computeMoaEnrichment(res) 
names(moa)

## Plot the drugs connectivity using plotDrugConnectivity()
plotDrugConnectivity(res, contr=1)

## Plot the mechanism of action using plotMOA()
plotMOA(moa, contr=k, type="drugClass", ntop=20)
plotMOA(moa, contr=k, type="targetGene", ntop=20)
    
## Plot the drugs activity map using plotActivationMap()
plotActivationMap(res, nterms = 60, nfc=20, rot=FALSE)
```
