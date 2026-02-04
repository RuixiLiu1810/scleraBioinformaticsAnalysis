# Minimal entrypoint example
source("R/config.R")
source("R/io_table.R")
source("R/preprocess.R")
source("R/design.R")
source("R/de_deseq2.R")
source("R/de_limma.R")
source("R/idmap.R")
source("R/enrich.R")
source("R/viz.R")
source("R/export.R")
source("R/pipeline.R")

# Edit these paths
meta_file <- "meta.csv"
expr_file <- "tpm_matrix.tsv"    # TPM/FPKM
counts_file <- "raw_counts.tsv"  # Raw counts

# Example: TPM/FPKM
cfg <- make_cfg(
  species="hs",
  input_type="expr_table",
  outdir="outputs_tpm",
  contrasts=list(c("group","HL","L"))
)
cfg$expr <- list(scale="tpm", log_transform=TRUE)  # log2(x+1)

run_pipeline(cfg, input=list(file=expr_file, meta=meta_file, gene_col="SYMBOL"), engine="limma")
