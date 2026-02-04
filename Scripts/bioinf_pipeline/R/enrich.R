deg_to_entrez <- function(deg, cfg, direction = c("UP","DOWN","ALL"), allow_alias = TRUE) {
  direction <- match.arg(direction)
  d <- deg
  if (direction != "ALL") d <- d[d$change == direction, , drop=FALSE]
  syms <- unique(d$SYMBOL)
  mp <- symbol_to_entrez(syms, cfg, allow_alias = allow_alias)
  unique(mp$ENTREZID)
}

enrich_go <- function(deg, cfg, direction = c("UP","DOWN"), ont = "ALL") {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("Please install clusterProfiler")
  direction <- match.arg(direction)

  if (!isTRUE(cfg$has_orgdb)) {
    stop("GO enrichment requires an OrgDb. For Cavia porcellus (cpoc), provide/build an OrgDb or skip GO.")
  }

  OrgDb <- ensure_orgdb(cfg)
  entrez <- deg_to_entrez(deg, cfg, direction = direction)

  if (length(entrez) < 5) stop("Too few genes for GO enrichment.")
  clusterProfiler::enrichGO(
    gene = entrez,
    OrgDb = OrgDb,
    keyType = "ENTREZID",
    ont = ont,
    pAdjustMethod = "BH",
    readable = TRUE
  )
}

enrich_kegg <- function(deg, cfg, direction = c("UP","DOWN")) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("Please install clusterProfiler")
  direction <- match.arg(direction)
  OrgDb <- ensure_orgdb(cfg)
  entrez <- deg_to_entrez(deg, cfg, direction = direction)

  if (length(entrez) < 5) stop("Too few genes for KEGG enrichment.")
  ek <- clusterProfiler::enrichKEGG(
    gene = entrez,
    organism = cfg$kegg_org,
    pAdjustMethod = "BH"
  )
  clusterProfiler::setReadable(ek, OrgDb = OrgDb, keyType = "ENTREZID")
}

make_rank_from_deg <- function(deg, score = c("logFC","signed_logp")) {
  score <- match.arg(score)
  d <- deg[!is.na(deg$padj) & !is.na(deg$logFC), , drop=FALSE]
  if (nrow(d) == 0) stop("No valid rows to build rank.")
  if (score == "logFC") {
    r <- d$logFC
  } else {
    r <- sign(d$logFC) * (-log10(pmax(d$padj, 1e-300)))
  }
  names(r) <- d$SYMBOL
  sort(r, decreasing = TRUE)
}
