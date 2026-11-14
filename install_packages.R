# Alternative to conda: install everything with BiocManager.
# Run once:  Rscript install_packages.R
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(c(
  "airway", "DESeq2", "apeglm", "AnnotationDbi", "org.Hs.eg.db",
  "clusterProfiler", "ggplot2", "ggrepel", "pheatmap"
), update = FALSE, ask = FALSE)
