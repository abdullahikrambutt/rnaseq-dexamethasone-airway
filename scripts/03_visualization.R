# =============================================================================
# 03_visualization.R
# MA plot, volcano plot, heatmap of top genes, and per-gene count plots.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

dds    <- readRDS("results/intermediate/dds_fitted.rds")
vsd    <- readRDS("results/intermediate/vsd.rds")
res_df <- readRDS("results/intermediate/res_df.rds")

alpha      <- 0.05
lfc_cutoff <- 1

res_df$status <- "Not significant"
res_df$status[!is.na(res_df$padj) & res_df$padj < alpha &
              res_df$log2FC_shrunk >  lfc_cutoff] <- "Up in dexamethasone"
res_df$status[!is.na(res_df$padj) & res_df$padj < alpha &
              res_df$log2FC_shrunk < -lfc_cutoff] <- "Down in dexamethasone"
status_colours <- c("Up in dexamethasone"   = "#D55E00",
                    "Down in dexamethasone" = "#0072B2",
                    "Not significant"       = "grey75")

# --- 1. MA plot: fold change vs mean expression ------------------------------
p_ma <- ggplot(res_df, aes(log10(baseMean), log2FC_shrunk, colour = status)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = status_colours) +
  labs(x = "log10 mean normalised count", y = "Shrunken log2 fold change",
       colour = NULL, title = "MA plot: dexamethasone vs untreated") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("figures/05_ma_plot.png", p_ma, width = 7, height = 5, dpi = 150)

# --- 2. Volcano plot ----------------------------------------------------------
plot_df <- subset(res_df, !is.na(padj))
top_labels <- head(plot_df[plot_df$status != "Not significant", ], 15)

p_volcano <- ggplot(plot_df, aes(log2FC_shrunk, -log10(padj), colour = status)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  geom_text_repel(data = top_labels, aes(label = symbol),
                  colour = "black", size = 3, max.overlaps = 30) +
  scale_colour_manual(values = status_colours) +
  labs(x = "Shrunken log2 fold change", y = "-log10 adjusted p-value",
       colour = NULL, title = "Volcano plot: dexamethasone vs untreated") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("figures/06_volcano.png", p_volcano, width = 7, height = 6, dpi = 150)

# --- 3. Heatmap of the top 30 differentially expressed genes ---------------
top30 <- head(plot_df$gene_id[plot_df$status != "Not significant"], 30)
mat <- assay(vsd)[top30, ]
mat <- t(scale(t(mat)))   # z-score each gene across samples
rownames(mat) <- res_df$symbol[match(top30, res_df$gene_id)]

annotation <- data.frame(Treatment = vsd$dex, Donor = vsd$cell,
                         row.names = colnames(vsd))
pheatmap(mat,
         annotation_col = annotation,
         main = "Top 30 differentially expressed genes (row z-scores)",
         filename = "figures/07_top_genes_heatmap.png",
         width = 8, height = 9)

# --- 4. Normalised counts for the top 6 genes --------------------------------
top6 <- head(plot_df$gene_id[plot_df$status != "Not significant"], 6)
counts_long <- do.call(rbind, lapply(top6, function(g) {
  d <- plotCounts(dds, gene = g, intgroup = c("dex", "cell"), returnData = TRUE)
  d$symbol <- res_df$symbol[res_df$gene_id == g]
  d
}))
counts_long$symbol <- factor(counts_long$symbol,
                             levels = res_df$symbol[match(top6, res_df$gene_id)])

p_counts <- ggplot(counts_long, aes(dex, count, colour = cell, group = cell)) +
  geom_point(size = 2.5) +
  geom_line(alpha = 0.5) +       # lines connect the same donor
  scale_y_log10() +
  facet_wrap(~ symbol, scales = "free_y") +
  labs(x = NULL, y = "Normalised count (log scale)", colour = "Donor",
       title = "Top genes: paired samples per donor") +
  theme_bw()
ggsave("figures/08_top_gene_counts.png", p_counts, width = 8, height = 5.5, dpi = 150)

cat("Plots written to figures/\n")
