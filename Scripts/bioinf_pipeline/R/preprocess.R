collapse_duplicated_genes <- function(mat, method = c("sum","mean")) {
  method <- match.arg(method)
  if (!any(duplicated(rownames(mat)))) return(mat)
  idx <- split(seq_len(nrow(mat)), rownames(mat))
  agg <- lapply(idx, function(ii) {
    if (length(ii) == 1) return(mat[ii, , drop=FALSE])
    if (method == "sum") colSums(mat[ii, , drop=FALSE]) else colMeans(mat[ii, , drop=FALSE])
  })
  out <- do.call(rbind, agg)
  out
}

prep_counts <- function(obj, min_count = 10, min_samples = 2) {
  stopifnot(!is.null(obj$counts))
  mat <- obj$counts
  mat <- collapse_duplicated_genes(mat, "sum")

  keep <- rowSums(mat >= min_count) >= min_samples
  mat <- mat[keep, , drop = FALSE]

  if (any(mat < 0, na.rm=TRUE)) stop("Counts contains negative values.")
  frac <- mean(abs(mat - round(mat)) > 1e-6)
  if (frac > 0.05) warning("Counts has many non-integers; ensure this is raw counts, not normalized.")
  obj$counts <- mat
  obj
}

needs_log2 <- function(expr) {
  q <- stats::quantile(expr, probs=c(0.5, 0.9), na.rm=TRUE)
  (q[[2]] > 50 && q[[1]] > 10)
}

prep_expr <- function(obj, force_log = NA) {
  stopifnot(!is.null(obj$expr) || !is.null(obj$counts))
  if (is.null(obj$expr) && !is.null(obj$counts)) return(obj)

  expr <- obj$expr
  expr <- collapse_duplicated_genes(expr, "mean")

  do_log <- if (is.na(force_log)) needs_log2(expr) else isTRUE(force_log)
  if (do_log) expr <- log2(expr + 1)

  obj$expr <- expr
  obj
}
