# Run the full analysis from the project root:  Rscript run_all.R
steps <- c("scripts/01_load_and_qc.R",
           "scripts/02_differential_expression.R",
           "scripts/03_visualization.R",
           "scripts/04_pathway_enrichment.R")
for (s in steps) {
  message("\n========== ", s, " ==========")
  source(s, local = new.env())
}
message("\nAll steps finished. See figures/ and results/.")
