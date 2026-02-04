get_species_map <- function(species = c("hs","mm","cpoc")) {
  species <- match.arg(species)

  if (species == "hs") {
    list(species="hs", orgdb_pkg="org.Hs.eg.db", kegg="hsa", has_orgdb=TRUE)
  } else if (species == "mm") {
    list(species="mm", orgdb_pkg="org.Mm.eg.db", kegg="mmu", has_orgdb=TRUE)
  } else if (species == "cpoc") {
    # Guinea pig: KEGG uses 'cpoc' in pathway IDs and organism pages. :contentReference[oaicite:1]{index=1}
    list(species="cpoc", orgdb_pkg=NULL, kegg="cpoc", has_orgdb=FALSE)
  }
}

make_cfg <- function(
  species = c("hs","mm","cpoc"),
  input_type = c("geo","counts_table","expr_table","seurat"),
  outdir = "outputs",
  group_col = "group",
  design = "~ group",
  contrasts = list(c("group","A","B")),
  cutoff = list(padj = 0.05, logFC = 1)
){
  sp <- get_species_map(match.arg(species))
  cfg <- list(
    species = sp$species,
    orgdb_pkg = sp$orgdb_pkg,
    kegg_org = sp$kegg,
    input = list(type = match.arg(input_type)),
    design = list(group_col = group_col, formula = design, contrasts = contrasts),
    cutoff = cutoff,
    outdir = outdir
  )
  cfg
}
