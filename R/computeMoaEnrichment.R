#' Get the mechanism of action object
#'
#' @param res is the output object from `computeConnectivityEnrichment()`
#' @param annot is a drug annotation table with `moa` and `target` columns
#'
#' @return MoA Mechanism of Action object
#' @export
#' @importFrom fgsea  fgsea
#'
#' @examples # from the data-sets provided as examples within the package load the .rda files
#' # DrugsAnnot, mDrugEnrich, mFC
#' res <- computeConnectivityEnrichment(mFC)
#' moa <- computeMoaEnrichment(res)
computeMoaEnrichment <- function(res, annot = metaLINCS::DrugsAnnot ) {

    ##annot = metaLINCS::DrugsAnnot

    ## --------------- attach annotation
    annot$drug <- annot$pert_iname
    annot <- annot[match(rownames(res$X), annot$drug), ]
    rownames(annot) <- rownames(res$X)
    ##Matrix::head(annot)
    annot <- annot[, c("drug", "moa", "target")]
    
    nc <- ncol(res$X)
    nn <- colnames(res$X)
    moa <- list()
    for(i in 1:nc) {
        message("Computing MOA enrichment for ",nn[i]," ...")
        moa[[i]] <- .computeMoaEnrichment1(res, annot, i)
    }
    names(moa) <- nn
    return(moa)

}

.getActiveDSEA <- function(res, annot, contr)
{
    nes <- round(res$X[, contr], 4)
    pv <- round(res$P[, contr], 4)
    qv <- round(res$Q[, contr], 4)
    drug <- rownames(res$X)
    nes[is.na(nes)] <- 0
    qv[is.na(qv)] <- 1
    pv[is.na(pv)] <- 1
    
    ## compile results matrix
    jj <- match(toupper(drug), toupper(rownames(annot)))
    annot1 <- annot[jj, c("moa", "target")]
    dt <- data.frame(drug = drug, NES = nes, pval = pv, padj = qv, annot1)    
    ##sel <- which(dt$moa != "" | dt$target != "")
    ##dt <- dt[sel, , drop = FALSE]    
    return(dt)
}


.computeMoaEnrichment1 <- function(res, annot, select) {
    ## meta-GSEA on MOA terms
        
    dt <- .getActiveDSEA(res, annot, contr=select)
    
    moa.list <- lapply(
        as.character(dt$moa),
        function(s) trimws(strsplit(s, split = "[\\|;,]")[[1]])
    )
    
    names(moa.list) <- rownames(dt)
    moa <- setdiff(unlist(moa.list), c("", NA, " "))
    gmt <- lapply(moa, function(g) names(which(sapply(moa.list, function(t) (g %in% t)))))
    names(gmt) <- moa
    
    targets.list <- lapply(
        as.character(dt$target),
        function(s) trimws(strsplit(s, split = "[\\|;,]")[[1]])
    )
    names(targets.list) <- rownames(dt)
    targets <- setdiff(unique(unlist(targets.list)), c(NA, "", " "))
    gmtar <- lapply(targets, function(g) {
        names(which(sapply(targets.list, function(t) (g %in% t))))
    })
    names(gmtar) <- targets
    
    rnk <- dt$NES
    names(rnk) <- rownames(dt)
    
    ## get drug class
    suppressWarnings(
        moa.class <- fgsea::fgsea(gmt, rnk, nperm = 20000)
    )
    moa.class <- moa.class[order(-abs(moa.class$NES)), ]
    
    ## get gene targets
    suppressWarnings(
        moa.target <- fgsea::fgsea(gmtar, rnk, nperm = 20000)
    )
    moa.target <- moa.target[order(-abs(moa.target$NES)), ]
    
    moa <- list("drugClass" = moa.class, "targetGene" = moa.target)
    
    return(moa)
}
