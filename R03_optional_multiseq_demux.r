apply_barcode_map <- function(cfg, obj, meta_col = "MULTIseq_call",
                              map_path = "config/barcode_map.csv",
                              out_col = "sample") {
  p <- here::here(map_path)
  if (!file.exists(p)) return(obj)

  map <- read.csv(p, stringsAsFactors = FALSE)
  stopifnot(all(c("barcode_id", "sample_label") %in% colnames(map)))

  m <- setNames(map$sample_label, map$barcode_id)
  obj[[out_col]] <- unname(m[obj[[meta_col, drop = TRUE]]])

  obj
}