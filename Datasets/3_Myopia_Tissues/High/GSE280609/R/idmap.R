ensure_orgdb <- function(cfg) {
  if (isTRUE(cfg$has_orgdb) && !is.null(cfg$orgdb_pkg)) {
    if (!requireNamespace(cfg$orgdb_pkg, quietly = TRUE)) {
      stop(sprintf("Missing OrgDb package: %s", cfg$orgdb_pkg))
    }
    return(get(cfg$orgdb_pkg))
  }
  return(NULL)
}

symbol_to_entrez <- function(symbols, cfg, allow_alias = TRUE) {

  # ---- Guinea pig branch ----
  if (cfg$species == "cpoc") {
    # biomaRt only (SYMBOL -> ENTREZ)
    return(symbol_to_entrez_biomart(symbols))
  }

  # ---- hs/mm original branch ----
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) stop("Please install AnnotationDbi")
  OrgDb <- ensure_orgdb(cfg)

  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]

  m1 <- AnnotationDbi::select(
    OrgDb,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )
  m1 <- m1[!is.na(m1$ENTREZID), c("SYMBOL","ENTREZID"), drop=FALSE]

  if (!allow_alias) return(unique(m1))

  mapped_sym <- unique(m1$SYMBOL)
  rest <- setdiff(symbols, mapped_sym)
  if (length(rest) == 0) return(unique(m1))

  m2 <- AnnotationDbi::select(
    OrgDb,
    keys = rest,
    keytype = "ALIAS",
    columns = c("ALIAS", "ENTREZID")
  )
  colnames(m2)[colnames(m2)=="ALIAS"] <- "SYMBOL"
  m2 <- m2[!is.na(m2$ENTREZID), c("SYMBOL","ENTREZID"), drop=FALSE]

  unique(rbind(m1, m2))
}

symbol_to_entrez_biomart <- function(symbols) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Please install biomaRt for Cavia porcellus mapping.")
  }

  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]

  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "cporcellus_gene_ensembl")
  # Ensembl supports Cavia porcellus and provides its annotation resources. :contentReference[oaicite:3]{index=3}

  # 常见属性：external_gene_name(=SYMBOL), entrezgene_id
  mp <- biomaRt::getBM(
    attributes = c("external_gene_name", "entrezgene_id"),
    filters    = "external_gene_name",
    values     = symbols,
    mart       = mart
  )

  colnames(mp) <- c("SYMBOL", "ENTREZID")
  mp <- mp[!is.na(mp$ENTREZID) & nzchar(mp$SYMBOL), , drop=FALSE]
  unique(mp)
}