# =============================================================================
# 02_differential_expression.R
# Test for genes whose expression changes with dexamethasone using DESeq2,
# shrink fold changes, and annotate Ensembl IDs with gene symbols.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

dds <- readRDS("results/intermediate/dds_filtered.rds")

# --- 1. Fit the model ---------------------------------------------------------
# DESeq() estimates size factors (sequencing depth), gene-wise dispersions
# (biological variability), and fits a negative binomial GLM for each gene.
dds <- DESeq(dds)
print(resultsNames(dds))

png("figures/04_dispersion.png", width = 900, height = 700, res = 150)
plotDispEsts(dds, main = "Dispersion estimates")
dev.off()

# --- 2. Extract results: treated vs untreated --------------------------------
alpha <- 0.05
res <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = alpha)
summary(res)

# --- 3. Shrink log2 fold changes ---------------------------------------------
# Genes with low counts can show huge but unreliable fold changes. apeglm
# shrinks these towards zero, giving better estimates for ranking and plotting.
res_shrunk <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")

# --- 4. Annotate with gene symbols -------------------------------------------
symbols <- mapIds(org.Hs.eg.db,
                  keys      = rownames(res),
                  keytype   = "ENSEMBL",
                  column    = "SYMBOL",
                  multiVals = "first")

res_df <- data.frame(
  gene_id          = rownames(res),
  symbol           = unname(symbols),
  baseMean         = res$baseMean,
  log2FC           = res$log2FoldChange,
  log2FC_shrunk    = res_shrunk$log2FoldChange,
  lfcSE            = res$lfcSE,
  pvalue           = res$pvalue,
  padj             = res$padj,
  stringsAsFactors = FALSE
)
# Fall back to the Ensembl ID where no symbol exists
res_df$symbol[is.na(res_df$symbol)] <- res_df$gene_id[is.na(res_df$symbol)]
res_df <- res_df[order(res_df$padj), ]

# --- 5. Define significant genes ---------------------------------------------
# Significant = adjusted p < 0.05 AND at least a 2-fold change (|log2FC| > 1).
lfc_cutoff <- 1
sig <- subset(res_df, !is.na(padj) & padj < alpha & abs(log2FC_shrunk) > lfc_cutoff)

cat(sprintf("\nGenes tested: %d\n", nrow(res_df)))
cat(sprintf("Significant (padj < %.2f, |log2FC| > %d): %d  (%d up, %d down)\n",
            alpha, lfc_cutoff, nrow(sig),
            sum(sig$log2FC_shrunk > 0), sum(sig$log2FC_shrunk < 0)))
cat("\nTop 15 genes by adjusted p-value:\n")
print(head(res_df[, c("symbol", "baseMean", "log2FC_shrunk", "padj")], 15),
      row.names = FALSE)

# --- 6. Save -----------------------------------------------------------------
write.csv(res_df, "results/deseq2_all_genes.csv", row.names = FALSE)
write.csv(sig,    "results/deseq2_significant_genes.csv", row.names = FALSE)
saveRDS(dds,    "results/intermediate/dds_fitted.rds")
saveRDS(res_df, "results/intermediate/res_df.rds")
cat("\nResults written to results/\n")
