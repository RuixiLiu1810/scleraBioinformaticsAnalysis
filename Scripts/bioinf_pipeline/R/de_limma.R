de_limma <- function(obj, cfg, contrast) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Please install limma")

  expr <- obj$expr
  if (is.null(expr)) stop("obj$expr is required for limma.")

  coldata <- make_coldata(obj, cfg$design$group_col)
  validate_contrast(coldata, contrast, cfg$design$group_col)

  group <- coldata[[cfg$design$group_col]]
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  A <- contrast[[2]]; B <- contrast[[3]]
  cont_mat <- limma::makeContrasts(contrasts = sprintf("%s-%s", A, B), levels = design)

  fit <- limma::lmFit(expr, design)
  fit2 <- limma::contrasts.fit(fit, cont_mat)
  fit2 <- limma::eBayes(fit2, trend = TRUE)

  tab <- limma::topTable(fit2, number = Inf, sort.by = "P")
  tab$SYMBOL <- rownames(tab)

  out <- data.frame(
    SYMBOL = tab$SYMBOL,
    logFC = tab$logFC,
    P.Value = tab$P.Value,
    padj = tab$adj.P.Val,
    stringsAsFactors = FALSE
  )
  out$change <- "NOT"
  out$change[!is.na(out$padj) & out$padj < cfg$cutoff$padj & out$logFC >= cfg$cutoff$logFC] <- "UP"
  out$change[!is.na(out$padj) & out$padj < cfg$cutoff$padj & out$logFC <= -cfg$cutoff$logFC] <- "DOWN"

  list(deg = out, fit = fit2)
}
