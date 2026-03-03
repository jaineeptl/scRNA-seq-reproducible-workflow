load_rna_with_optional_soupx <- function(cfg) {
  cr_dir <- here::here(cfg$tenx$cellranger_dir)
  stop_if_missing(cr_dir, "CellRanger directory")

  if (isTRUE(cfg$soupx$enabled)) {
    soup_ch <- SoupX::load10X(cr_dir)
    soup_ch <- SoupX::autoEstCont(soup_ch)
    rna_counts <- SoupX::adjustCounts(soup_ch, roundToInt = FALSE)
  } else {
    # works for 10x multiome output folder too
    x <- Seurat::Read10X(data.dir = cr_dir)
    # multiome can return list; RNA matrix usually named "Gene Expression"
    if (is.list(x)) {
      rna_counts <- x[["Gene Expression"]] %||% x[[1]]
    } else {
      rna_counts <- x
    }
  }

  obj <- Seurat::CreateSeuratObject(counts = rna_counts, assay = "RNA")
  obj
}