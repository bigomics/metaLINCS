#' Compute Connectivity Enrichment
#'
#' @param mFC A matrix of differential gene expression fold-change of own experiment, the rows name of the matrix must be the genes names
#' @param names Names in case mFC is a vector
#' @param mDrugEnrich a large matrix represents the gene expression fold change in the presence of different perturbations, either various drugs concentrations or other genetic engineering techniques as gene-knock-in etc..
#' @param nmin the minimum number of experiments in the drug set
#' @param nprune  takes only the (nprune) top matching drugs for each comparison to reduce the size of the matrix. default is 0 to take the full matrix
#' @return list of drugs enrichment statistics and annotations
#' @export
#' @import stats
#' @importFrom Matrix head
#' @importFrom fgsea  fgseaSimple
#'
#' @examples # from the data-sets provided as examples within the package load the .rda files
#' # DrugsAnnot, mDrugEnrich, mFC
#' res <- computeConnectivityEnrichment(mFC, mDrugEnrich, nmin=15, nprune = 250)
computeConnectivityEnrichment <- function(mFC, names=NULL,
                                          mDrugEnrich = SpaceLINCS::mDrugEnrich,
                                          nmin = 15, nprune = 250)
{
    
    if(NCOL(mFC)==1 && !is.null(names(mFC))) {
        mFC <- matrix(mFC, ncol=1, dimnames=list(names(mFC),'FC'))
    }
    if(NCOL(mFC)==1 && is.null(names(mFC)) && !is.null(names)) {
        mFC <- matrix(mFC, ncol=1, dimnames=list(names,'FC'))
    }
    if(NCOL(mFC)==1 && is.null(names(mFC)) && is.null(names)) {
        stop("must supply names")
    }
    
    xdrugs <- gsub("[@_].*$", "", colnames(mDrugEnrich))
    ndrugs <- length(table(xdrugs))

    message("number of profiles: ", ncol(mDrugEnrich))
    message("number of drugs: ", ndrugs)

    ## create drug meta sets
    meta.gmt <- tapply(colnames(mDrugEnrich), xdrugs, list)
    meta.gmt <- meta.gmt[which(sapply(meta.gmt, length) >= nmin)]
    length(meta.gmt)
    if (length(meta.gmt) == 0) {
        message("WARNING::: computeDrugEnrichment : no valid genesets!!")
        return(NULL)
    }

    ## first level (rank) correlation
    message("Calculating first level rank correlation ...")
    gg <- intersect(rownames(mDrugEnrich), rownames(mFC))
    length(gg)
    if (length(gg) < 20) {
        message("WARNING::: computeDrugEnrichment : not enough common genes!!")
        return(NULL)
    }
    rnk1 <- apply(mDrugEnrich[gg, , drop = FALSE], 2, rank, na.last = "keep")
    rnk2 <- apply(mFC[gg, , drop = FALSE], 2, rank, na.last = "keep")
    R1 <- stats::cor(rnk1, rnk2, use = "pairwise")
    dim(R1)

    R1 <- R1 + 1e-8 * matrix(stats::rnorm(length(R1)), nrow(R1), ncol(R1))
    colnames(R1) <- colnames(mFC)
    rownames(R1) <- colnames(mDrugEnrich)
    dim(R1)

    ## experiment to drug
    res <- list()
    message("Calculating connectivity enrichment using GSEA ...")
    res0 <- list()
    i <- 1
    for (i in 1:ncol(R1)) {
        suppressWarnings(res0[[i]] <- fgsea::fgseaSimple(meta.gmt, stats = R1[, i], nperm = 1000))
    }
    names(res0) <- colnames(R1)
    length(res0)
    
    mNES <- sapply(res0, function(x) x$NES)
    mQ <- sapply(res0, function(x) x$padj)
    mP <- sapply(res0, function(x) x$pval)
    if (length(res0) == 1) {
        mNES <- cbind(mNES)
        mP <- cbind(mP)
        mQ <- cbind(mQ)
    }
    
    pw <- res0[[1]]$pathway
    rownames(mNES) <- rownames(mQ) <- rownames(mP) <- pw
    colnames(mNES) <- colnames(mQ) <- colnames(mP) <- colnames(mFC)
    msize <- res0[[1]]$size
    dim(R1)
    res <- list(X = mNES, Q = mQ, P = mP, size = msize)
    names(res)

    ## this takes only the top matching drugs for each comparison to
    ## reduce the size of the matrices
    if (nprune > 0) {
        message(" pruning : nprune = ", nprune)
        ## reduce solution set with top-N of each comparison
        ##mtop <- apply(abs(res$X), 2, function(x) Matrix::head(order(-x), nprune)) ## absolute NES!!
        rx   <- apply(abs(res$X), 2, order)
        rownames(rx) <- rownames(res$X)
        mtop <- names(head(sort(rowMeans(rx), decreasing=TRUE),ntop))
        length(mtop)
        top.idx <- unique(unlist(meta.gmt[mtop]))
        length(top.idx)
        res$X <- res$X[mtop, , drop = FALSE]
        res$P <- res$P[mtop, , drop = FALSE]
        res$Q <- res$Q[mtop, , drop = FALSE]
        res$size <- res$size[mtop]
    }

    ## statistics (running metric)
    sel.drugs <- rownames(res$X)
    xsel <- which(xdrugs %in% sel.drugs)
    rp <- as.array(R1[xsel, ])
    res$stats <- rp
    res$drug <- xdrugs[xsel]
    names(res)
    
    return(res)
}
