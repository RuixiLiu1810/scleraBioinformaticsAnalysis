`%||%` <- function(x, y) if (!is.null(x)) x else y

run_pipeline <- function(cfg, input, engine = NULL) {
  if (is.null(engine)) {
    engine <- if (cfg$input$type == "counts_table") "deseq2" else "limma"
  }
  outdir <- ensure_outdir(cfg$outdir)

  obj <- switch(cfg$input$type,
    counts_table = io_counts_table(input$file, input$meta, gene_col = input$gene_col %||% "SYMBOL"),
    expr_table   = io_expr_table(input$file, input$meta, gene_col = input$gene_col %||% "SYMBOL"),
    stop("Unsupported input type")
  )

  if (cfg$input$type == "counts_table") obj <- prep_counts(obj)
  if (cfg$input$type == "expr_table") {
    # TPM/FPKM default: log2(x+1)
    force_log <- isTRUE(cfg$expr$log_transform)
    obj <- prep_expr(obj, force_log = force_log)
  }

  results <- list()
  for (i in seq_along(cfg$design$contrasts)) {
    ct <- cfg$design$contrasts[[i]]
    tag <- paste0(ct[[2]], "_vs_", ct[[3]])

    if (engine == "deseq2") {
      res <- de_deseq2(obj, cfg, ct)
      deg <- res$deg
      save_rds(res$expr_vst, outdir, paste0(tag, "_expr_vst"))
    } else if (engine == "limma") {
      res <- de_limma(obj, cfg, ct)
      deg <- res$deg
    } else {
      stop("Unsupported engine: ", engine)
    }

    write_deg(deg, outdir, paste0(tag))

    pvol <- plot_volcano(deg, cfg, top_n = 10)
    save_plot(pvol, outdir, paste0(tag, "_volcano"))

    # Enrichment attempts (non-fatal)
    try({
      ego_up <- enrich_go(deg, cfg, "UP")
      save_rds(ego_up, outdir, paste0(tag, "_GO_UP"))
      save_plot(plot_enrich_dot(ego_up, 15, paste0(tag, " GO UP")), outdir, paste0(tag, "_GO_UP_dot"))
    }, silent = TRUE)

    try({
      ek_up <- enrich_kegg(deg, cfg, "UP")
      save_rds(ek_up, outdir, paste0(tag, "_KEGG_UP"))
      save_plot(plot_enrich_dot(ek_up, 15, paste0(tag, " KEGG UP")), outdir, paste0(tag, "_KEGG_UP_dot"))
    }, silent = TRUE)

    results[[tag]] <- list(deg = deg)
  }

  invisible(results)
}
