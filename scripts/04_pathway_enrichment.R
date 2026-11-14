# =============================================================================
# 04_pathway_enrichment.R
# Gene Ontology (Biological Process) over-representation analysis for genes
# up- and down-regulated by dexamethasone.
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

res_df <- readRDS("results/intermediate/res_df.rds")

alpha      <- 0.05
lfc_cutoff <- 1

tested <- subset(res_df, !is.na(padj))
up   <- tested$gene_id[tested$padj < alpha & tested$log2FC_shrunk >  lfc_cutoff]
down <- tested$gene_id[tested$padj < alpha & tested$log2FC_shrunk < -lfc_cutoff]

# The "universe" is every gene we actually tested, not the whole genome.
# Using the whole genome would inflate significance.
universe <- tested$gene_id

run_go <- function(genes, label) {
  cat(sprintf("Running GO enrichment for %d %s genes...\n", length(genes), label))
  ego <- enrichGO(gene          = genes,
                  universe      = universe,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENSEMBL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)   # convert IDs to gene symbols

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat(sprintf("  No significant GO terms for %s genes.\n", label))
    return(invisible(NULL))
  }

  write.csv(as.data.frame(ego),
            sprintf("results/go_enrichment_%s.csv", label), row.names = FALSE)

  p <- dotplot(ego, showCategory = 15) +
    ggtitle(sprintf("GO Biological Process: %s-regulated genes", label))
  ggsave(sprintf("figures/09_go_%s.png", label), p, width = 8, height = 7, dpi = 150)

  cat("  Top terms:\n")
  print(head(as.data.frame(ego)[, c("Description", "Count", "p.adjust")], 10),
        row.names = FALSE)
  invisible(ego)
}

run_go(up,   "up")
run_go(down, "down")
cat("Enrichment results written to results/ and figures/\n")
