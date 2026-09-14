## Runs GSEA using hpgsea if installed (much faster), otherwise falls
## back to fgsea::fgseaSimple. Always returns a data.frame with columns
## set, set_size, NES, p_value, adj_p_value.
.runGSEA <- function(gene_sets, stats, nperm) {
    if (requireNamespace("hpgsea", quietly = TRUE)) {
        res <- suppressWarnings(
            hpgsea::hpgsea(stats = stats, gene_sets = gene_sets, nperm = nperm, sort = FALSE)
        )
        return(res[, c("set", "set_size", "NES", "p_value", "adj_p_value")])
    }

    if (!requireNamespace("fgsea", quietly = TRUE)) {
        stop("Either 'hpgsea' or 'fgsea' must be installed to run GSEA.")
    }
    res <- suppressWarnings(fgsea::fgseaSimple(gene_sets, stats = stats, nperm = nperm))
    data.frame(
        set = res$pathway, set_size = res$size, NES = res$NES,
        p_value = res$pval, adj_p_value = res$padj,
        stringsAsFactors = FALSE
    )
}
