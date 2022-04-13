---
output:
  pdf_document: 
    toc: yes
    number_sections: yes
  html_document: default
---

\title{SpaceLINCS Reference Manual}
\author{Ivo Kwee and Layal Abo Khayal  \thanks{funded by the BigOmics}}

\maketitle

<!-- toc -->

Januar 03, 2022

## DESCRIPTION

    Package: SpaceLINCS
    Title: Stratified Perturbation Analysis by Correlation Enrichment of
        connectivity
    Version: 0.0.0.9000
    Authors@R: c(
        person("IVO", "KWEE", , "kwee@bigomics.ch", role = "aut", comment = c(ORCID = "0000-0002-2751-4218")),
        pperson("Layal", "Abo Khayal", , "Layal.abo.khayal@gmail.com", role = "cre")
        )
    Description: 'SpaceLINCS' calculate and visualise the correlations
        between expermental gene expresion profiles and the signatures from
        databases that have focused on detailing gene expression profiles
        based on the exposure of different cell lines to various perturbagen,
        such as varying concentration of various compounds, or identifying
        signatures that correlate with drug responses, or even the alteration
        in gene expression profile caused by silencing or activating a gene or
        a set of genes.
    License: GPL (>= 3)
    URL: https://public.bigomics.ch/app/omicsplayground_contrib
    Imports: 
        corrplot,
        fgsea,
        gplots,
        graphics,
        Matrix,
        qlcMatrix,
        stats,
        utils,
        uwot,
        WGCNA
    Encoding: UTF-8
    LazyData: true
    Depends: R (>= 2.10)
    LazyDataCompression: bzip2
    Roxygen: list(markdown = TRUE)
    RoxygenNote: 7.1.2

# computePeturbEnrichment

Compute drug Enrichment, drugs activity or drugs sensitivity based on database and annotation provided as input parameter

## Description

Compute drug Enrichment, drugs activity or drugs sensitivity based on database and annotation provided as input parameter

## Usage

``` r
computePeturbEnrichment(
  mFC,
  mDrugEnrich,
  DrugsAnnot,
  methods = c("GSEA", "cor"),
  nmin = 15,
  nprune = 0,
  contrast = NULL
)
```

## Arguments

+---------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Argument      | Description                                                                                                                                                                                                |
+===============+============================================================================================================================================================================================================+
| `mFC`         | = A matrix of differential gene expression fold-change of own experiment, the rows name of the matrix must be the genes names                                                                              |
+---------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `mDrugEnrich` | = a large matrix represents the gene expression fold change in the presence of different perturbations, either various drugs concentrations or other genetic engineering techniques as gene-knock-in etc.. |
+---------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `DrugsAnnot`  | =is matrix represents the drugs annotation that contains the targets, mechanism of action , clinical phase, disease area, and indication                                                                   |
+---------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `methods`     | the computation methods of Enrichment, either using GSEA algorithm or the rank correlation.                                                                                                                |
+---------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `nmin`        | the minimum number of drugs                                                                                                                                                                                |
+---------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `nprune`      | takes only the (nprune) top matching drugs for each comparison to reduce the size of the matrix. default is 0 to take the full matrix                                                                      |
+---------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `contrast`    | is a character string represent the two compared conditions(contrast) as it is provided from the fold-change matrix                                                                                        |
+---------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

## Value

drugs a list of drugs enrichment stats and annotations

## Examples

``` r
# from the data-sets provided as examples within the package load the .rda files
# DrugsAnnot, mDrugEnrich, mFC
drugs <- computePeturbEnrichment(
mFC = mFC, mDrugEnrich = mDrugEnrich, DrugsAnnot = DrugsAnnot,
methods = "GSEA", nprune = 250, contrast = NULL
)
```

# DrugsAnnot

Annotation of a set of drugs

## Description

This data set is a matrix represents the drugs annotation that contains drugs' names, targets, mechanism of action , clinical phase, disease area, and indication.

## Format

A data frame with 6125 observations on the following 6 variables. "pert_iname","clinical_phase", "moa", "target", "disease_area" and indication

## Usage

``` r
data("DrugsAnnot")
```

## Details

L1000 repurposing drugs . For research use only. Do not use the Repurposing Hub to make clinical treatment decisions. BROAD DOES NOT GUARANTEE OR WARRANT THE ACCURACY OF THE DATA WITHIN THE REPURPOSING HUB.

## References

Corsello SM, et al. Nature Medicine. 2017 Apr 7;23(4):405-408. doi: 10.1038/nm.4306

## Examples

``` r
data(DrugsAnnot)
## maybe str(DrugsAnnot) ; plot(DrugsAnnot) ...
```

# dseaPlotActmap

Plot drug Activity Map

## Description

Plot drug Activity Map

## Usage

``` r
dseaPlotActmap(dsea, drugs, method = "GSEA", contr, nterms = 50, nfc = 20)
```

## Arguments

+------------+---------------------------------------------------------------------------------------------------------------------+
| Argument   | Description                                                                                                         |
+============+=====================================================================================================================+
| `dsea`     | drugs set enrichment analysis object as output of getActiveDSEA                                                     |
+------------+---------------------------------------------------------------------------------------------------------------------+
| `drugs`    | a list of drugs enrichment stats and annotations the output of computePeturbEnrichment()                            |
+------------+---------------------------------------------------------------------------------------------------------------------+
| `method`   | the computation methods of Enrichment, either using GSEA algorithm or the rank correlation.                         |
+------------+---------------------------------------------------------------------------------------------------------------------+
| `contr`    | is a character string represent the two compared conditions(contrast) as it is provided from the fold-change matrix |
+------------+---------------------------------------------------------------------------------------------------------------------+
| `nterms`   | integer                                                                                                             |
+------------+---------------------------------------------------------------------------------------------------------------------+
| `nfc`      | integer                                                                                                             |
+------------+---------------------------------------------------------------------------------------------------------------------+

## Value

plot of drug activity map

## Examples

``` r
# from the data-sets provided as examples within the package load the .rda files
#DrugsAnnot, mDrugEnrich, mFC
dsea = getActiveDSEA(mDrugEnrich, DrugsAnnot, mFC, contr = colnames(mFC)[1])
drugs= computePeturbEnrichment(mFC = mFC, mDrugEnrich = mDrugEnrich, DrugsAnnot = DrugsAnnot,
methods = "GSEA", nprune = 250, contrast = NULL)
dseaPlotActmap (dsea, drugs, method = "GSEA", contr=colnames(mFC)[1], nterms = 50, nfc=20)
```

# getActiveDSEA

Create the drug set enrichment object that is used in the visualization

## Description

Create the drug set enrichment object that is used in the visualization

## Usage

``` r
getActiveDSEA(mDrugEnrich, DrugsAnnot, mFC, contr = NULL)
```

## Arguments

+---------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Argument      | Description                                                                                                                                                            |
+===============+========================================================================================================================================================================+
| `mDrugEnrich` | a large matrix represents the gene expression fold change in the presence of different perturbations, either various drugs concentrations or other genetic engineering |
+---------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `DrugsAnnot`  | is matrix represents the drugs annotation that contains the targets, mechanism of action , clinical phase, disease area, and indication                                |
+---------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `mFC`         | A matrix of differential gene expression fold-change of own experiment, the rows name of the matrix must be the genes names                                            |
+---------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `contr`       | is a character string represent the two compared conditions(contrast) as it is provided from the fold-change matrix                                                    |
+---------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

## Value

drug set enrichment object

## Examples

``` r
# from the data-sets provided as examples within the package load the .rda files
# DrugsAnnot, mDrugEnrich, mFC
dsea <- getActiveDSEA(mDrugEnrich, DrugsAnnot, mFC, contr = colnames(mFC)[1])
```

# getMOA

Get the mechanism of action object

## Description

Get the mechanism of action object

## Usage

``` r
getMOA(dsea)
```

## Arguments

+---------------+------------------------------------------------------------------------------+
| Argument      | Description                                                                  |
+===============+==============================================================================+
| `dsea`        | is a drug set enrichment object as an output of the function getActiveDSEA() |
+---------------+------------------------------------------------------------------------------+

## Value

MoA Mechanism of Action object

## Examples

``` r
# from the data-sets provided as examples within the package load the .rda files
# DrugsAnnot, mDrugEnrich, mFC
dsea <- getActiveDSEA(mDrugEnrich, DrugsAnnot, mFC, contr = colnames(mFC)[1])
Moa <- getMOA(dsea)
```

# mDrugEnrich

The drugs activity on 1001 genes

## Description

a large matrix represents the gene expression fold change of 1001 genes in the presence of 20220 different drugs concentrations and working time.

## Format

The format is: num \[1:1001, 1:20220\] -12.7 0.82 -1.48 -1.89 4.02 ... - attr(\*, "dimnames")=List of 2 ..\$ : chr \[1:1001\] "PSME1" "CISD1" "SPDEF" "ATF1" ... ..\$ : chr \[1:20220\] "crizotinib_HEPG2_10um_24h" "avicin-d_JURKAT_10um_24h" "BRD-K35716340_WSUDLCL2_10um_6h" "trequinsin_HCC515_0.4um_24h" ...

## Usage

``` r
data("mDrugEnrich")
```

## Details

The LINCS L1000 project has collected gene expression profiles for thousands of perturbagens at a variety of time points, doses, and cell lines. A full list of the chemical and genetic perturbations used can be found on the CLUE website along with their descriptions.

## References

Duan Q et al. LINCS Canvas Browser: interactive web app to query, browse and interrogate LINCS L1000 gene expression signatures. Nucleic Acids Res. 2014 Jul;42 (Web Server issue):W449-60.

## Examples

``` r
data(mDrugEnrich)
## maybe str(mDrugEnrich) ; plot(mDrugEnrich) ...
```

# mFC

Fold change matrix of differential gene expresion from RNAseq data of Multiple myeloma

## Description

There are 6 contrasts : 1- glucocorticoid sensitive vs glucocorticoid resistant ("Resistant.vs.Sensitive").2- Treated with Withaferin-A vs untreated ("WithaferinA.vs.Untreated"). 3- untreated glucocorticoid resistant cells vs untreated glucocorticoid sensitive cells ("UTR.vs.UTS"). 4- Treated with Withaferin-A glucocorticoid resistant cells vs Treated with Withaferin-A glucocorticoid sensitive cells ("WAR.vs.WAS"). 5- Treated Resistant vs Untreated Resistant ("WAR.vs.UTR"). 6- Treated sensitive vs Untreated sensitive ("WAS.vs.UTS")

## Format

A data frame with 5332 observations on the following 6 variables. Resistant.vs.Sensitive, "WithaferinA.vs.Untreated, UTR.vs.UTS, WAR.vs.WAS , WAR.vs.UTR, WAS.vs.UTS.

## Usage

``` r
data("mFC")
```

## Details

The RNAseq data from (Logie et al 2021) study. They compared the therapeutic efficacy of the phytochemical kinase inhibitor withaferin A with the clinically approved BTK inhibitor ibrutinib to target hyperactivated tyrosine kinase signaling in glucocorticoid-resistant multiple myeloma cells. theur results demonstrate that withaferin-A induced cell death of glucocorticoid-resistant MM1R cells involves covalent cysteine targeting of multiple Hinge-6 domain type tyrosine kinases of the kinase cysteinome classification, including BTK.

## References

Logie E, Chirumamilla CS, Perez-Novo C, Shaw P et al. Covalent Cysteine Targeting of Bruton's Tyrosine Kinase (BTK) Family by Withaferin-A Reduces Survival of Glucocorticoid-Resistant Multiple Myeloma MM1 Cell. Cancers (Basel) 2021 Mar 31;13(7). PMID: 33807411

## Examples

``` r
data(mFC)
## maybe str(mFC) ; plot(mFC) ...
```

# `plotDrugConnectivity`

Plot Drug Connectivity

## Description

Plot Drug Connectivity

## Usage

``` r
plotDrugConnectivity(dsea)
```

## Arguments

| Argument | Description                                                    |
|----------|----------------------------------------------------------------|
| `dsea`   | is a drug set Enrichment object is the output of getActiveDSEA |

## Value

a plot of drug connectivity

## Examples

``` r
# from the data-sets provided as examples within the package load the .rda files
# DrugsAnnot, mDrugEnrich, mFC
dsea <- getActiveDSEA(mDrugEnrich, DrugsAnnot, mFC, contr = colnames(mFC)[1])
plotDrugConnectivity(dsea)
```

# plotMOA

Plot mechanism of action

## Description

Plot mechanism of action

## Usage

``` r
plotMOA(dsea, ntop = 16)
```

## Arguments

+---------------+------------------------------------------------------------------------------+
| Argument      | Description                                                                  |
+===============+==============================================================================+
| `dsea`        | is dsea object which is the output of the function getActiveDSEA()           |
+---------------+------------------------------------------------------------------------------+
| `ntop`        | the number of the top enteries (genes or drug), ntop = 16 as a default value |
+---------------+------------------------------------------------------------------------------+

## Value

plot of mechanism of action

## Examples

``` r
# from the data-sets provided as examples within the package load the .rda files
# DrugsAnnot, mDrugEnrich, mFC
dsea <- getActiveDSEA(mDrugEnrich, DrugsAnnot, mFC, contr = colnames(mFC)[1])
plotMOA(dsea)
```
