#' Plot drug Activity Map
#'
#' @param dsea drugs set enrichment analysis object as output of getActiveDSEA
#' @param drugs a list of drugs enrichment stats and annotations the output of computePeturbEnrichment()
#' @param method the computation methods of Enrichment, either using GSEA algorithm or the rank correlation.
#' @param contr is a character string represent the two compared conditions(contrast) as it is provided from the fold-change matrix
#' @param nterms integer
#' @param nfc  integer
#'
#' @return plot of drug activity map
#' @export
#' @import stats
#' @import utils
#' @importFrom corrplot corrplot
#' @importFrom gplots bluered
#'
#' @examples# from the data-sets provided as examples within the package load the .rda files
#' #DrugsAnnot, mDrugEnrich, mFC
#' dsea = getActiveDSEA(mDrugEnrich, DrugsAnnot, mFC, contr = colnames(mFC)[1])
#' drugs= computePeturbEnrichment(mFC = mFC, mDrugEnrich = mDrugEnrich, DrugsAnnot = DrugsAnnot,
#' methods = "GSEA", nprune = 250, contrast = NULL)
#' plotActivationMap (dsea, drugs, method = "GSEA", contr=colnames(mFC)[1], nterms = 50, nfc=20)
#'
plotActivationMap <- function(res, nterms = 60, nfc = 20, rot=FALSE)
{
    ##nterms=50;nfc=20
    nes <-res$X
    qv <- res$Q
    score <- nes * (1 - qv)**2    
    score[is.na(score)] <- 0
    dim(score)
    
    ## reduce score matrix
    score <- score[order(-rowMeans(score**2)),,drop = FALSE] ## sort by score
    
    if (nrow(score) <= 1) {
        return(NULL)
    }
    
    score <- utils::head(score, nterms) ## max number of terms
    score <- score[, utils::head(order(-colSums(score**2)), nfc), drop = FALSE] ## max contrs/FC
    score <- score + 1e-3 * matrix(stats::rnorm(length(score)), nrow(score), ncol(score))
    
    ## Normalize columns of the activation matrix.
    score <- t(t(score) / (1e-8 + sqrt(colMeans(score**2))))
    
    score <- sign(score) * abs(score)**3 ## fudging
    score <- score / (1e-8 + max(abs(score), na.rm = TRUE))
    
    if (NCOL(score) > 1) {
        d1 <- stats::as.dist(1 - stats::cor(t(score), use = "pairwise"))
        d2 <- stats::as.dist(1 - stats::cor(score, use = "pairwise"))        
        d1[is.na(d1)] <- 1
        d2[is.na(d2)] <- 1
        jj <- 1
        ii <- 1:nrow(score)
        ii <- stats::hclust(d1)$order
        jj <- stats::hclust(d2)$order
        score <- score[ii, jj, drop = FALSE]
    } else {
        score <- score[order(-score[, 1]), , drop = FALSE]
    }

    cex2 <- 1
    colnames(score) <- substring(colnames(score), 1, 30)
    rownames(score) <- substring(rownames(score), 1, 50)
    cex2 <- 0.85
    par(mfrow = c(1, 1), mar = c(1, 1, 1, 1), oma = c(0, 1, 0, 0))

    if(rot) {
        ## vertical plot
        corrplot::corrplot(score,
                           is.corr = FALSE, cl.pos = "n", col = gplots::bluered(100),
                           tl.cex = 0.9 * cex2, tl.col = "grey20", tl.srt = 90
                           )
    } else {
        ## horizontal plot
        corrplot::corrplot(t(score),
                           is.corr = FALSE, cl.pos = "n", col = gplots::bluered(100),
                           tl.cex = 0.9 * cex2, tl.col = "grey20", tl.srt = 90
                           )

    }
}
