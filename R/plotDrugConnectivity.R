#' Plot Drug Connectivity
#'
#' @param dsea is a drug set Enrichment object is the output of getActiveDSEA
#'
#' @return a plot of drug connectivity
#' @export
#' @import graphics
#' @import utils
#'
#' @examples # from the data-sets provided as examples within the package load the .rda files
#' # DrugsAnnot, mDrugEnrich, mFC
#' dsea <- getActiveDSEA(mDrugEnrich, DrugsAnnot, mFC, contr = colnames(mFC)[1])
#' plotDrugConnectivity(dsea)
plotDrugConnectivity <- function(res, contr, drugs=NULL, nplots=16) {

    dt <- selectResult(res, contr)

    if (nrow(dt) == 0) {
        return(NULL)
    }
    
    rnk <- res$stats[,contr]
    if (length(rnk) == 0) {
        return(NULL)
    }

    dt <- dt[order(-abs(dt$NES)),]

    if(!is.null(drugs)) {
        sel <- intersect(drugs,rownames(dt))
        dt <- dt[sel,]
    }
    
    ## set layout
    dt <- utils::head(dt, nplots)
    lab.cex <- 0.75
    xlab <- ""
    ylab <- ""
    nc <- ceiling(sqrt(nrow(dt)))
    par(oma = c(0, 1.9, 0, 0))
    par(mfrow = c(nc, nc), mar = c(0.8, 1.0, 1.3, 0.8), mgp = c(1.9, 0.6, 0))

    ## generate plot
    i <- 1
    for (i in 1:nrow(dt)) {
        dx <- rownames(dt)[i]
        
        gmtdx <- grep(dx, names(rnk), fixed = TRUE, value = TRUE) ## L1000 naming
        length(gmtdx)
        
        dx1 <- substring(dx, 1, 26)
        par(cex.axis = 0.001)
        if (i %% nc == 1) par(cex.axis = 0.98)
        suppressWarnings(
            gsea.enplot(rnk, gmtdx,
                        main = dx1, cex.main = 1.2,
                        xlab = xlab, ylab = ylab
                        )
        )
        nes <- round(dt$NES[i], 2)
        qv <- round(dt$padj[i], 4)
        tt <- c(paste("NES=", nes), paste("q=", qv))
        graphics::legend("topright", legend = tt, cex = 0.8, y.intersp = 0.85, bty = "n")
        if (i %% nc == 1 && nrow(dt) > 1) {
            graphics::mtext("rank metric", side = 2, line = 1.8, cex = lab.cex)
        }
    }
    
}



#' Select result
#'
#' @param res is the result from computeConnectivityEnrichment()
#' @param contr contrast selection
#'
#' @return selected result for selected contrast
#' @export
#' 
#' @examples
#' res1 <- selectResult(res,1)
selectResult <- function(res, contr)
{
    nes <- round(res$X[, contr], 4)
    pv <- round(res$P[, contr], 4)
    qv <- round(res$Q[, contr], 4)
    drug <- rownames(res$X)
    nes[is.na(nes)] <- 0
    qv[is.na(qv)] <- 1
    pv[is.na(pv)] <- 1
    
    ## compile results matrix
    dt <- data.frame(drug = drug, NES = nes, pval = pv, padj = qv)    
    return(dt)
}

