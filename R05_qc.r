calc_mad_bounds <- function(x, k = 3) {
  med <- stats::median(x, na.rm = TRUE)
  md  <- stats::mad(x, na.rm = TRUE)
  c(lower = med - k * md, upper = med + k * md)
}

rna_qc <- function(cfg, obj) {
  Seurat::DefaultAssay(obj) <- "RNA"

  obj[["percent.mt"]]   <- Seurat::PercentageFeatureSet(obj, pattern = cfg$species$mito_pattern)
  obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(obj, pattern = cfg$species$ribo_pattern)

  md <- cfg$qc$rna$mad_k
  mito_cap <- cfg$qc$rna$mito_cap %||% Inf

  df <- Seurat::FetchData(obj, vars = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"))
  b_nf <- calc_mad_bounds(df$nFeature_RNA, md)
  b_nc <- calc_mad_bounds(df$nCount_RNA,   md)
  b_mt <- calc_mad_bounds(df$percent.mt,   md)
  b_rb <- calc_mad_bounds(df$percent.ribo, md)
  b_mt["upper"] <- min(b_mt["upper"], mito_cap)

  p <- ggpubr::ggarrange(
    Seurat::VlnPlot(obj, "percent.ribo", pt.size = 0.1) + ggplot2::geom_hline(yintercept = b_rb, linetype="dashed"),
    Seurat::VlnPlot(obj, "percent.mt",   pt.size = 0.1) + ggplot2::geom_hline(yintercept = b_mt, linetype="dashed"),
    Seurat::VlnPlot(obj, "nFeature_RNA", pt.size = 0.1) + ggplot2::geom_hline(yintercept = b_nf, linetype="dashed"),
    Seurat::VlnPlot(obj, "nCount_RNA",   pt.size = 0.1) + ggplot2::geom_hline(yintercept = b_nc, linetype="dashed"),
    ncol = 2, nrow = 2
  )
  save_plot(cfg, "qc_rna.png", p, width = 10, height = 8, dpi = 300)

  obj <- subset(obj, subset =
    nFeature_RNA > b_nf["lower"] & nFeature_RNA < b_nf["upper"] &
    nCount_RNA   > b_nc["lower"] & nCount_RNA   < b_nc["upper"] &
    percent.mt   > b_mt["lower"] & percent.mt   < b_mt["upper"] &
    percent.ribo > b_rb["lower"] & percent.ribo < b_rb["upper"]
  )

  obj
}

atac_qc <- function(cfg, obj) {
  Signac::DefaultAssay(obj) <- "ATAC"
  obj <- Signac::NucleosomeSignal(obj)
  obj <- Signac::TSSEnrichment(obj, fast = TRUE)

  md <- cfg$qc$atac$mad_k
  df <- Seurat::FetchData(obj, vars = c("nCount_ATAC","TSS.enrichment","nucleosome_signal","peak_region_fragments"))

  b_nc  <- calc_mad_bounds(df$nCount_ATAC, md)
  b_tss <- calc_mad_bounds(df$TSS.enrichment, md)
  b_nuc <- calc_mad_bounds(df$nucleosome_signal, md)
  b_prf <- calc_mad_bounds(df$peak_region_fragments, md)

  p <- ggpubr::ggarrange(
    Seurat::VlnPlot(obj, "TSS.enrichment", pt.size = 0.1) + ggplot2::geom_hline(yintercept = b_tss, linetype="dashed"),
    Seurat::VlnPlot(obj, "nucleosome_signal", pt.size = 0.1) + ggplot2::geom_hline(yintercept = b_nuc, linetype="dashed"),
    Seurat::VlnPlot(obj, "peak_region_fragments", pt.size = 0.1) + ggplot2::geom_hline(yintercept = b_prf, linetype="dashed"),
    Seurat::VlnPlot(obj, "nCount_ATAC", pt.size = 0.1) + ggplot2::geom_hline(yintercept = b_nc, linetype="dashed"),
    ncol = 2, nrow = 2
  )
  save_plot(cfg, "qc_atac.png", p, width = 10, height = 8, dpi = 300)

  obj <- subset(obj, subset =
    nCount_ATAC         > b_nc["lower"]  & nCount_ATAC         < b_nc["upper"] &
    TSS.enrichment      > b_tss["lower"] & TSS.enrichment      < b_tss["upper"] &
    nucleosome_signal   > b_nuc["lower"] & nucleosome_signal   < b_nuc["upper"] &
    peak_region_fragments > b_prf["lower"] & peak_region_fragments < b_prf["upper"]
  )

  obj
}