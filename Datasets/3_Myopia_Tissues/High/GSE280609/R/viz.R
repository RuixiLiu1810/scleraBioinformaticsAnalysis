plot_volcano <- function(deg, cfg, top_n = 10) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Please install ggrepel")

  d <- deg
  d$neglog10 <- -log10(pmax(d$padj, 1e-300))
  d$label <- ""

  idx <- order(d$padj, -abs(d$logFC))
  idx <- idx[!is.na(idx)]
  idx <- head(idx, top_n)
  d$label[idx] <- d$SYMBOL[idx]

  ggplot2::ggplot(d, ggplot2::aes(x=logFC, y=neglog10)) +
    ggplot2::geom_point(ggplot2::aes(color=change), alpha=0.8) +
    ggrepel::geom_text_repel(ggplot2::aes(label=label), max.overlaps = Inf, size=3) +
    ggplot2::theme_bw() +
    ggplot2::labs(x="log2FC", y="-log10(padj)")
}

plot_heatmap_top <- function(expr, deg, top_n = 50) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Please install pheatmap")
  stopifnot(!is.null(expr))
  d <- deg[deg$change != "NOT" & !is.na(deg$padj), , drop=FALSE]
  d <- d[order(d$padj), , drop=FALSE]
  sel <- unique(d$SYMBOL)[1:min(top_n, nrow(d))]
  mat <- expr[intersect(sel, rownames(expr)), , drop=FALSE]
  if (nrow(mat) < 2) stop("Too few genes for heatmap.")
  pheatmap::pheatmap(mat, scale="row", show_rownames=TRUE, show_colnames=TRUE)
}

plot_enrich_dot <- function(enrich_res, show = 15, title = NULL) {
  if (!requireNamespace("enrichplot", quietly = TRUE)) stop("Please install enrichplot")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
  p <- enrichplot::dotplot(enrich_res, showCategory = show)
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}
