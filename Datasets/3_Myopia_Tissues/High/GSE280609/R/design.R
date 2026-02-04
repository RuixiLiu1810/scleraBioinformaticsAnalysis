make_coldata <- function(obj, group_col = "group") {
  meta <- obj$meta
  stopifnot(group_col %in% colnames(meta))
  cd <- meta
  cd[[group_col]] <- factor(cd[[group_col]])
  rownames(cd) <- cd$sample_id
  cd
}

validate_contrast <- function(coldata, contrast, group_col = "group") {
  stopifnot(length(contrast) == 3)
  stopifnot(contrast[[1]] == group_col)
  lv <- levels(coldata[[group_col]])
  if (!(contrast[[2]] %in% lv && contrast[[3]] %in% lv)) {
    stop(sprintf("Contrast levels not in group levels. group levels: %s", paste(lv, collapse=", ")))
  }
  TRUE
}
