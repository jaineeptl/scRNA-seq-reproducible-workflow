run_atac_lsi <- function(cfg, obj, assay = "ATAC") {
  Seurat::DefaultAssay(obj) <- assay
  obj <- Signac::FindTopFeatures(obj, min.cutoff = 5)
  obj <- Signac::RunTFIDF(obj)
  obj <- Signac::RunSVD(obj)
  obj
}

run_wnn <- function(cfg, obj) {
  dims_rna <- cfg$analysis$dims_rna[[1]]:cfg$analysis$dims_rna[[2]]
  dims_lsi <- cfg$analysis$dims_lsi[[1]]:cfg$analysis$dims_lsi[[2]]

  obj <- Seurat::FindMultiModalNeighbors(
    object = obj,
    reduction.list = list("pca", "lsi"),
    dims.list = list(dims_rna, dims_lsi),
    modality.weight.name = "RNA.weight"
  )

  obj <- Seurat::RunUMAP(obj, nn.name = "weighted.nn", assay = "RNA",
                         umap.method = "uwot", reduction.name = "wnn.umap",
                         reduction.key = "wnnUMAP_")

  obj <- Seurat::FindClusters(obj, graph.name = "wsnn", algorithm = 3,
                              resolution = cfg$analysis$resolution_wnn, verbose = FALSE)
  obj
}