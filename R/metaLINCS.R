#' metaLINCS: an R package for enrichment analysis of LINCS L1000 drug signatures using meta-level connectivity mapping
#'
#' The metaLINCS package provides three categories of important functions:
#' computePeturbEnrichment, computeComboEnrichment and several plotting functions
#'
#' @section metaLINCS functions:
#'
#'
#' computeConnectivityEnrichment: Compute Connectivity Enrichment  based on database and annotation provided as input parameter
#' USAGE:
#' computeConnectivityEnrichment(mFC = mFC, mDrugEnrich = mDrugEnrich, nprune = 250, contrast = NULL)
#'
#' plotActivationMap: Plot drug Activity Map.
#' USAGE:
#' plotActivationMap(res, nterms = 50, nfc = 20)
#'
#' computeMoaEnrichment: Compute mechanism-of-action enrichment
#' USAGE:
#' computeMoaEnrichment(dsea)
#'
#' plotDrugConnectivity: Plot Drug Connectivity
#' USAGE:
#' plotDrugConnectivity(dsea)
#'
#' plotMOA: Plot mechanism of action bargraph
#' USAGE:
#' plotMOA(dsea, ntop = 16)
#'
#'
#' @docType package
#' @name metaLINCS
NULL
