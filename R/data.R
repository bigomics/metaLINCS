#' Drug annotation table
#'
#' Mechanism-of-action (MoA) and target annotation for the perturbagens in
#' \code{mDrugEnrich}, keyed by \code{pert_iname}.
#'
#' @format A data frame with 6125 rows and 6 columns:
#' \describe{
#'   \item{pert_iname}{drug/perturbagen name}
#'   \item{clinical_phase}{clinical development phase}
#'   \item{moa}{mechanism(s) of action, "|"-separated}
#'   \item{target}{gene target(s), "|"-separated}
#'   \item{disease_area}{associated disease area}
#'   \item{indication}{clinical indication}
#' }
#' @source LINCS L1000 / Connectivity Map (clue.io) drug metadata
"DrugsAnnot"

#' LINCS L1000 drug perturbation signatures
#'
#' Gene expression fold-change matrix for thousands of drug and genetic
#' perturbation profiles from the LINCS L1000 project, used as the reference
#' connectivity database in \code{computeConnectivityEnrichment()}.
#'
#' @format A matrix with 1001 genes (rows) and 20220 perturbation profiles
#'   (columns). Column names encode the perturbagen, e.g.
#'   \code{"<drug>_<cell line>_<dose/time>"}.
#' @source LINCS L1000 project (clue.io)
"mDrugEnrich"

#' Example differential gene expression fold-changes
#'
#' Example fold-change matrix from an in-house experiment, used throughout
#' the package examples and vignette as the query signature for
#' \code{computeConnectivityEnrichment()}.
#'
#' @format A data frame with 5332 genes (rows) and 6 contrasts (columns).
"mFC"
