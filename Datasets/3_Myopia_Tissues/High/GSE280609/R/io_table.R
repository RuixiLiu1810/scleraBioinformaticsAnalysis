read_meta <- function(meta_file) {
  meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot("sample_id" %in% colnames(meta))
  meta
}

align_by_meta <- function(mat, meta) {
  stopifnot(!is.null(colnames(mat)))
  keep <- intersect(meta$sample_id, colnames(mat))
  if (length(keep) == 0) stop("No overlap between meta$sample_id and matrix colnames.")
  meta2 <- meta[match(keep, meta$sample_id), , drop = FALSE]
  mat2  <- mat[, keep, drop = FALSE]
  list(mat = mat2, meta = meta2)
}

io_counts_table <- function(counts_file, meta_file, gene_col = "SYMBOL") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table")
  df <- data.table::fread(counts_file, data.table = FALSE)
  stopifnot(gene_col %in% colnames(df))
  rn <- df[[gene_col]]
  mat <- as.matrix(df[, setdiff(colnames(df), gene_col), drop = FALSE])
  rownames(mat) <- rn

  meta <- read_meta(meta_file)
  aligned <- align_by_meta(mat, meta)

  list(counts = aligned$mat, expr = NULL, meta = aligned$meta, feature = NULL)
}

io_expr_table <- function(expr_file, meta_file, gene_col = "SYMBOL") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table")
  df <- data.table::fread(expr_file, data.table = FALSE)
  stopifnot(gene_col %in% colnames(df))
  rn <- df[[gene_col]]
  mat <- as.matrix(df[, setdiff(colnames(df), gene_col), drop = FALSE])
  rownames(mat) <- rn

  meta <- read_meta(meta_file)
  aligned <- align_by_meta(mat, meta)

  list(counts = NULL, expr = aligned$mat, meta = aligned$meta, feature = NULL)
}
