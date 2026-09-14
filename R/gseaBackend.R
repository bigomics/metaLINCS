## Runs GSEA using hpgsea if installed (much faster), otherwise falls
## back to fgsea::fgseaSimple. Always returns a data.frame with fgsea's
## column names: pathway, pval, padj, NES, size.
.emptyGSEA <- function() {
    data.frame(
        pathway = character(0), pval = numeric(0), padj = numeric(0),
        NES = numeric(0), size = integer(0), stringsAsFactors = FALSE
    )
}

.runGSEA <- function(gene_sets, stats, nperm) {
    if (requireNamespace("hpgsea", quietly = TRUE)) {
        ## hpgsea requires min_size >= 2 and errors out (instead of just
        ## dropping sets, like fgsea does) when nothing meets it
        res <- tryCatch(
            suppressWarnings(
                hpgsea::hpgsea(stats = stats, gene_sets = gene_sets, nperm = nperm, sort = FALSE)
            ),
            error = function(e) {
                if (grepl("min_size", conditionMessage(e))) return(NULL)
                stop(e)
            }
        )
        if (is.null(res)) return(.emptyGSEA())
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
