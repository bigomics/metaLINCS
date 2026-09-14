## Runs GSEA using hpgsea if installed (much faster), otherwise falls
## back to fgsea::fgseaSimple. Always returns a data.frame with fgsea's
## column names: pathway, pval, padj, NES, size.
.runGSEA <- function(gene_sets, stats, nperm) {
    if (requireNamespace("hpgsea", quietly = TRUE)) {
        res <- suppressWarnings(
            hpgsea::hpgsea(stats = stats, gene_sets = gene_sets, nperm = nperm, sort = FALSE)
        )
        return(data.frame(
            pathway = res$set, pval = res$p_value, padj = res$adj_p_value,
            NES = res$NES, size = res$set_size,
            stringsAsFactors = FALSE
        ))
    }

    if (!requireNamespace("fgsea", quietly = TRUE)) {
        stop("Either 'hpgsea' or 'fgsea' must be installed to run GSEA.")
    }
    res <- suppressWarnings(fgsea::fgseaSimple(gene_sets, stats = stats, nperm = nperm))
    res[, c("pathway", "pval", "padj", "NES", "size")]
}
