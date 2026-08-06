#!/usr/bin/env Rscript
# =============================================================================
# OA_RA_Bcells.R
#
# OA vs RA synovium B-cell analysis:
#   1. Re-embed and cluster B cells
#   2. Marker visualization (including CXCL8)
#
# BCR processing is in OA_RA_BCR.R
#
# Custom packages developed for this study:
#   bioscrna  Single-cell RNA-seq analysis toolkit built on Bioconductor
#             (data loading, QC, dimension reduction, clustering, and plotting
#             for SingleCellExperiment objects).
#             https://github.com/zhangzmr/bioscrna
#             remotes::install_github("zhangzmr/bioscrna")
# =============================================================================

suppressPackageStartupMessages({
  library(bioscrna)
  library(BiocParallel)
  library(patchwork)
})

# Project root: set this to the repository root before running
# setwd("/path/to/aNav_code")

dir.create("object/OA_RA", showWarnings = FALSE, recursive = TRUE)
dir.create("plot/OA_RA/bcells", showWarnings = FALSE, recursive = TRUE)
dir.create("tables/OA_RA", showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1. Re-embed and cluster B cells
# -----------------------------------------------------------------------------
sce <- readRDS("object/OA_RA/bcells.rds")
sce <- scranNormalize(sce, graph.fun = "louvain")
sce <- dimensionReduction(sce, n = 2000, rank = 50)
sce <- merging(
  sce, batch = "batch", dim = 50,
  method = c("FastMNN", "Harmony"), cores = 16
)
saveRDS(sce, file = "object/OA_RA/bcells_reduced.rds")

k_range <- seq(10, 50, by = 2)
sce <- clustering(
  sce, k_range, method = "louvain",
  dimred = "batch_corrected", min_cl = 3, cores = 16, top = 20
)
saveRDS(sce, file = "object/OA_RA/bcells_clustered.rds")

sce <- plotCluster(
  sce, dimred = "UMAP_corrected",
  file = "plot/OA_RA/bcells/UMAP_clusters.pdf",
  width = 8, height = 8, label = TRUE
)

pdf("plot/OA_RA/bcells/UMAP_batch.pdf", width = 8, height = 8)
reducedPlot(sce, var = "batch", reduction.use = "UMAP_corrected",
            size = 0.5, do.raster = TRUE, do.label = FALSE, main = "Samples")
if ("disease" %in% colnames(colData(sce))) {
  reducedPlot(sce, var = "disease", reduction.use = "UMAP_corrected",
              size = 0.5, do.raster = TRUE, do.label = FALSE, main = "Disease")
}
dev.off()


# -----------------------------------------------------------------------------
# 2. Marker visualization
# -----------------------------------------------------------------------------
bcell_genes <- c(
  "IGHM", "CD27", "JCHAIN", "CD38", "MZB1", "XBP1", "PDIA6", "IGHG1", "IGHG3",
  "CD24", "ITGAM", "IGHD", "TCL1A", "CD1C", "AICDA", "BCL6", "CD69", "JUN",
  "ITGAX", "LAMP1", "TBX21", "CXCR5", "CR2", "FCRL5", "CD19", "CXCL8"
)
intergene <- intersect(bcell_genes, rownames(sce))

p <- plot_density(sce, intergene, size = 0.5, reduction = "UMAP_corrected")
pdf("plot/OA_RA/bcells/marker_density.pdf", width = 30, height = 20)
print(p + plot_layout(ncol = 6))
dev.off()

if ("CXCL8" %in% rownames(sce)) {
  pdf("plot/OA_RA/bcells/UMAP_CXCL8.pdf", width = 6, height = 6)
  print(plot_density(sce, "CXCL8", size = 0.5, reduction = "UMAP_corrected"))
  dev.off()
}

message("Done. B-cell object: object/OA_RA/bcells_clustered.rds")
message("Done. Plots: plot/OA_RA/bcells/")
