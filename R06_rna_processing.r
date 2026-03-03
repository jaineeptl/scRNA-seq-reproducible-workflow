run_rna_sct <- function(cfg, obj) {
  Seurat::DefaultAssay(obj) <- "RNA"
  obj <- Seurat::NormalizeData(obj)

  # Cell cycle scoring (species-aware optional: user can override later)
  s.genes   <- Seurat::cc.genes$s.genes
  g2m.genes <- Seurat::cc.genes$g2m.genes
  obj <- Seurat::CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

  vars <- c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","S.Score","G2M.Score")
  vars <- vars[vars %in% colnames(obj@meta.data)]

  obj <- Seurat::SCTransform(obj, vars.to.regress = vars, verbose = FALSE)
  obj
}

run_rna_umap_cluster <- function(cfg, obj) {
  obj <- Seurat::RunPCA(obj)
  dims <- cfg$analysis$dims_rna[[1]]:cfg$analysis$dims_rna[[2]]
  obj <- Seurat::FindNeighbors(obj, dims = dims)
  obj <- Seurat::FindClusters(obj, resolution = cfg$analysis$resolution_rna)
  obj <- Seurat::RunUMAP(obj, dims = dims, umap.method = "uwot",
                         reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
  obj
}