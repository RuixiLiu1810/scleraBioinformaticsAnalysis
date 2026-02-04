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
counts_file <- "GSE280609_sclerarna_count.csv"  # Raw counts

# Example: TPM/FPKM
cfg <- make_cfg(
  species="cpoc",
  input_type="counts_table",
  outdir="outputs_count",
  contrasts=list(c("group","myopia","control"))
)
cfg$expr <- list(scale="raw", log_transform=FALSE)  # log2(x+1)

run_pipeline(cfg, input=list(file=counts_file, meta=meta_file, gene_col="gene_id"), engine="deseq2")
