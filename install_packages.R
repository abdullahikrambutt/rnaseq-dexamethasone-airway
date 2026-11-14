# Alternative to conda: install everything with BiocManager.
# Run once:  Rscript install_packages.R
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# GO.db is installed first because enrichplot needs it during installation

BiocManager::install("GO.db", update = FALSE, ask = FALSE)
BiocManager::install(c(
  "airway", "DESeq2", "apeglm", "AnnotationDbi", "org.Hs.eg.db",
  "clusterProfiler", "enrichplot", "ggplot2", "ggrepel", "pheatmap"
), update = FALSE, ask = FALSE)
