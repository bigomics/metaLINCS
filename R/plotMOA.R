#' Plot mechanism of action
#'
#' @param dsea is dsea object which is the output of the function getActiveDSEA()
#' @param ntop the number of the top enteries (genes or drug), ntop = 16 as a default value
#'
#' @return plot of mechanism of action
#' @export
#' @import utils
#' @import graphics
#'
#' @examples # from the data-sets provided as examples within the package load the .rda files
#' # DrugsAnnot, mDrugEnrich, mFC
#' dsea <- getActiveDSEA(mDrugEnrich, DrugsAnnot, mFC, contr = colnames(mFC)[1])
#' plotMOA(dsea)
plotMOA <- function(moa, contr=NULL, type=c("drugClass","targetGene"), ntop = 20) {

    if(length(moa)>1 && is.null(contr)) {
        stop("Multiple contrasts detected. You must select one contrast.")
    }

    type <- type[1]
    df <- moa[[contr]][[type]]
    
    jj <- unique(c(utils::head(order(-df$NES), ntop), utils::tail(order(-df$NES), ntop)))
    moa.top <- df$NES[jj]
    names(moa.top) <- df$pathway[jj]
    par(mfrow = c(2, 1), mar = c(3.5, 4, 3.5, 3.5), mgp = c(2.5, 1, 0.5), cex=1)
    graphics::barplot(moa.top,
                      horiz = FALSE, las = 3,
                      ylab = "Enrichment score  (NES)",
                      main = "Mechanism of action",
                      cex.names = 1.1
                      )
}
