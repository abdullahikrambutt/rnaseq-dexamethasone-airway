# =============================================================================
# 01_load_and_qc.R
# Load the airway RNA-seq dataset, filter low-count genes, and run
# sample-level quality control (library sizes, PCA, sample distances).
# =============================================================================

suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

dir.create("results/intermediate", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# --- 1. Load data -------------------------------------------------------------
# airway: 8 samples of human airway smooth muscle cells from 4 donors (cell
# lines), each either untreated or treated with 1 uM dexamethasone for 18 h.
# Counts are stored as a SummarizedExperiment: genes (Ensembl IDs) x samples.
data("airway")
se <- airway

# Make "untreated" the reference level, so positive log2 fold changes mean
# "higher in dexamethasone-treated cells".
se$dex  <- relevel(factor(se$dex), ref = "untrt")
se$cell <- factor(se$cell)

# Readable sample labels, e.g. "N61311_trt"
colnames(se) <- paste(se$cell, se$dex, sep = "_")

cat(sprintf("Loaded %d genes x %d samples\n", nrow(se), ncol(se)))
print(as.data.frame(colData(se))[, c("cell", "dex")])

# --- 2. Build the DESeq2 object ---------------------------------------------
# Paired design: '~ cell + dex' compares treated vs untreated *within* each
# donor, removing donor-to-donor variation from the treatment effect.
dds <- DESeqDataSet(se, design = ~ cell + dex)

# --- 3. Filter low-count genes ----------------------------------------------
# Keep genes with >= 10 reads in at least 4 samples (the size of the smallest
# group). This removes noise and speeds up the analysis.
keep <- rowSums(counts(dds) >= 10) >= 4
cat(sprintf("Keeping %d of %d genes after low-count filtering\n",
            sum(keep), length(keep)))
dds <- dds[keep, ]

# --- 4. Library sizes ---------------------------------------------------------
lib <- data.frame(
  sample   = colnames(dds),
  millions = colSums(counts(dds)) / 1e6,
  dex      = dds$dex
)
p_lib <- ggplot(lib, aes(x = reorder(sample, millions), y = millions, fill = dex)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Assigned reads (millions)", fill = "Treatment",
       title = "Library size per sample") +
  theme_bw()
ggsave("figures/01_library_sizes.png", p_lib, width = 6, height = 4, dpi = 150)

# --- 5. Variance-stabilising transformation + PCA ---------------------------
# Raw counts have variance that grows with the mean; VST makes samples
# comparable for distance-based methods like PCA and clustering.
vsd <- vst(dds, blind = TRUE)

pca <- plotPCA(vsd, intgroup = c("dex", "cell"), returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))
p_pca <- ggplot(pca, aes(PC1, PC2, colour = dex, shape = cell)) +
  geom_point(size = 4) +
  labs(x = paste0("PC1 (", pct[1], "% variance)"),
       y = paste0("PC2 (", pct[2], "% variance)"),
       colour = "Treatment", shape = "Donor",
       title = "PCA of variance-stabilised counts") +
  theme_bw()
ggsave("figures/02_pca.png", p_pca, width = 6, height = 4.5, dpi = 150)

# --- 6. Sample-to-sample distance heatmap ------------------------------------
sample_dists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sample_dists),
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         main = "Euclidean distance between samples (VST)",
         filename = "figures/03_sample_distances.png",
         width = 7, height = 6)

# --- 7. Save for the next step ----------------------------------------------
saveRDS(dds, "results/intermediate/dds_filtered.rds")
saveRDS(vsd, "results/intermediate/vsd.rds")
cat("QC done. Figures written to figures/\n")
