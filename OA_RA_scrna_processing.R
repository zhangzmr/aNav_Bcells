#!/usr/bin/env Rscript
# =============================================================================
# OA_RA_scrna_processing.R
#
# OA vs RA synovium scRNA-seq processing (all cells):
#   1. Load 10X data
#   2. Quality control
#   3. Normalization, dimension reduction, and batch integration
#   4. Clustering and disease labeling
#   5. Subset B cells
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
  library(bluster)
  library(patchwork)
  library(dittoSeq)
})

# Project root: set this to the repository root before running
# setwd("/path/to/aNav_code")

dir.create("object/OA_RA", showWarnings = FALSE, recursive = TRUE)
dir.create("plot/OA_RA", showWarnings = FALSE, recursive = TRUE)
dir.create("tables/OA_RA", showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# Expected layout: data/OA_RA/rna/<sample>/sample_filtered_feature_bc_matrix/
sample <- list.files("data/OA_RA/rna")
directory <- paste0("data/OA_RA/rna/", sample, "/sample_filtered_feature_bc_matrix")
sample_sheet <- data.frame(sample = sample, directory = directory)

sce <- loading_sce(directory = sample_sheet$directory, sample = sample_sheet$sample)
saveRDS(sce, file = "object/OA_RA/allcells_raw.rds")


# -----------------------------------------------------------------------------
# 2. Quality control
# -----------------------------------------------------------------------------
sce <- plotQC(sce, file = "plot/OA_RA/Quality_Control_before.pdf", width = 20, height = 5)
sce <- qualityfilter(
  sce,
  filter_outlier = TRUE,
  sum_max = 50000,
  sum_min = 800,
  gene_max = 8000,
  gene_min = 500,
  mito_percent = 20,
  filter_gene = TRUE,
  min_counts = 5,
  batch = sce$batch
)
sce <- plotQC(sce, file = "plot/OA_RA/Quality_Control_after.pdf", width = 17, height = 5)
saveRDS(sce, file = "object/OA_RA/allcells_filtered.rds")


# -----------------------------------------------------------------------------
# 3. Normalization, dimension reduction, and batch integration
# -----------------------------------------------------------------------------
sce <- scranNormalize(sce, graph.fun = "louvain")
sce <- dimensionReduction(sce, n = 2000, rank = 50)
sce <- merging(
  sce, batch = "batch", dim = 50,
  method = c("FastMNN", "Harmony"), cores = 16
)
saveRDS(sce, file = "object/OA_RA/allcells_reduced.rds")


# -----------------------------------------------------------------------------
# 4. Clustering and disease labels
# -----------------------------------------------------------------------------
sce$cluster_MNN <- clusterCells(
  sce,
  use.dimred = "batch_corrected",
  BLUSPARAM = NNGraphParam(
    k = 50, cluster.fun = "louvain", BPPARAM = MulticoreParam(16)
  )
)
sce$cluster_HARMONY <- clusterCells(
  sce,
  use.dimred = "HARMONY",
  BLUSPARAM = NNGraphParam(
    k = 10, cluster.fun = "louvain", BPPARAM = MulticoreParam(16)
  )
)

sce$disease[sce$batch %in% c(
  "OA103", "OA103B", "OA121", "OA128", "OA128B",
  "OA47", "OA47B", "OA55", "OA55T"
)] <- "OA"
sce$disease[sce$batch %in% c(
  "RA401B", "RA401T", "RA443B", "RA443T", "RA470",
  "RA510B", "RA510T", "RA563B", "RA563T"
)] <- "RA"
sce$disease <- factor(sce$disease, levels = c("OA", "RA"))

saveRDS(sce, file = "object/OA_RA/allcells_clustered.rds")

pdf("plot/OA_RA/allcells_UMAP_batch.pdf", width = 8, height = 8)
reducedPlot(sce, var = "batch", reduction.use = "UMAP", size = 0.3,
            do.raster = TRUE, do.label = FALSE, main = "UMAP")
reducedPlot(sce, var = "batch", reduction.use = "UMAP_HARMONY", size = 0.3,
            do.raster = TRUE, do.label = FALSE, main = "HARMONY")
reducedPlot(sce, var = "batch", reduction.use = "UMAP_corrected", size = 0.3,
            do.raster = TRUE, do.label = FALSE, main = "FastMNN")
dev.off()

pdf("plot/OA_RA/allcells_UMAP_clusters.pdf", width = 8, height = 8)
reducedPlot(sce, var = "cluster_MNN", reduction.use = "UMAP_corrected",
            size = 0.3, do.raster = TRUE, do.label = TRUE, main = "MNN Clusters")
reducedPlot(sce, var = "cluster_HARMONY", reduction.use = "UMAP_HARMONY",
            size = 0.3, do.raster = TRUE, do.label = TRUE, main = "HARMONY Clusters")
reducedPlot(sce, var = "disease", reduction.use = "UMAP_corrected",
            size = 0.3, do.raster = TRUE, do.label = FALSE, main = "Disease")
dev.off()

immune_genes <- c(
  "CD19", "IGHM", "IGHD", "CD27", "CD38", "CD24", "CR2", "CD3D", "CD4", "CD8A",
  "JCHAIN", "MS4A1", "CD3E", "CD14", "XBP1", "IGHG1", "IGHG3", "FCER1A",
  "FCGR3A", "LYZ", "GZMK", "GZMB", "CXCL8", "CXCR1", "CXCR2"
)
intergene <- intersect(immune_genes, rownames(sce))
p <- plot_density(sce, intergene, size = 0.5, reduction = "UMAP_corrected")
pdf("plot/OA_RA/allcells_marker_density.pdf", width = 30, height = 20)
print(p + plot_layout(ncol = 6))
dev.off()


# -----------------------------------------------------------------------------
# 5. Subset B cells
# -----------------------------------------------------------------------------
# Clusters identified from marker inspection (FastMNN / louvain)
bcells <- sce[, sce$cluster_MNN %in% c(7, 11, 17)]
saveRDS(bcells, file = "object/OA_RA/bcells.rds")

message("Done. All-cell object: object/OA_RA/allcells_clustered.rds")
message("Done. B-cell subset: object/OA_RA/bcells.rds")
