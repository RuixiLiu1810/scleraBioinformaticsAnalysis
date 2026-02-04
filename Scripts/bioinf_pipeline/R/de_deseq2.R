de_deseq2 <- function(obj, cfg, contrast) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Please install DESeq2")
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("Please install SummarizedExperiment")

  counts <- obj$counts
  if (is.null(counts)) stop("obj$counts is required for DESeq2.")

  coldata <- make_coldata(obj, cfg$design$group_col)
  validate_contrast(coldata, contrast, cfg$design$group_col)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = coldata,
    design = stats::as.formula(cfg$design$formula)
  )
  dds <- DESeq2::DESeq(dds)

  res <- DESeq2::results(dds, contrast = contrast)
  tab <- as.data.frame(res)
  tab$SYMBOL <- rownames(tab)

  out <- data.frame(
    SYMBOL = tab$SYMBOL,
    logFC = tab$log2FoldChange,
    P.Value = tab$pvalue,
    padj = tab$padj,
    stringsAsFactors = FALSE
  )
  out$change <- "NOT"
  out$change[!is.na(out$padj) & out$padj < cfg$cutoff$padj & out$logFC >= cfg$cutoff$logFC] <- "UP"
  out$change[!is.na(out$padj) & out$padj < cfg$cutoff$padj & out$logFC <= -cfg$cutoff$logFC] <- "DOWN"

  vst <- DESeq2::vst(dds, blind = TRUE)
  expr_vst <- SummarizedExperiment::assay(vst)

  list(deg = out, expr_vst = expr_vst, dds = dds)
}
