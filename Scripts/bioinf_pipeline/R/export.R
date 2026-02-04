ensure_outdir <- function(outdir) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  outdir
}

write_deg <- function(deg, outdir, prefix) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table")
  ensure_outdir(outdir)
  fn <- file.path(outdir, paste0(prefix, "_DEG.tsv"))
  data.table::fwrite(deg, fn, sep = "\t")
  fn
}

save_rds <- function(obj, outdir, name) {
  ensure_outdir(outdir)
  fn <- file.path(outdir, paste0(name, ".rds"))
  saveRDS(obj, fn)
  fn
}

save_plot <- function(p, outdir, name, w = 7, h = 5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
  ensure_outdir(outdir)
  fn <- file.path(outdir, paste0(name, ".pdf"))
  ggplot2::ggsave(fn, plot = p, width = w, height = h)
  fn
}
