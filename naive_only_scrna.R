#!/usr/bin/env Rscript
# =============================================================================
# naive_only_scrna.R
#
# Naive-only cohort scRNA-seq processing:
#   1. Load, QC, reduce, and cluster all cells
#   2. Subset naive B cells and re-embed
#   3. Annotate cell types and visualize
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
  library(dittoSeq)
  library(patchwork)
})

# Project root: set this to the repository root before running
# setwd("/path/to/aNav_code")

dir.create("object/naive_only", showWarnings = FALSE, recursive = TRUE)
dir.create("plot/naive_only", showWarnings = FALSE, recursive = TRUE)
dir.create("tables/naive_only", showWarnings = FALSE, recursive = TRUE)

celltype_colors <- c(
  "rNaive" = "#FF7F0E",
  "aNaive1" = "#1F77B4",
  "aNaive2" = "#D62728",
  "aNaive3" = "#2CA02C"
)

exclude_ig <- function(sce) {
  ig_filter <- c("IGLV", "IGHV", "IGKV", "IGLC")
  ig_index <- unlist(lapply(ig_filter, function(x) grep(x, rownames(sce))))
  sce[!rownames(sce) %in% rownames(sce)[ig_index], ]
}


# -----------------------------------------------------------------------------
# 1. Load, QC, reduce, and cluster all cells
# -----------------------------------------------------------------------------
# Expected layout: data/naive_only/gex/<sample>/filtered_feature_bc_matrix/
sample <- list.files("data/naive_only/gex")
directory <- paste0("data/naive_only/gex/", sample, "/filtered_feature_bc_matrix")

sce <- loading_sce(directory = directory, sample = sample)
saveRDS(sce, file = "object/naive_only/allcells_raw.rds")

sce <- plotQC(sce, file = "plot/naive_only/Quality_Control_before.pdf", width = 20, height = 5)
sce <- qualityfilter(
  sce,
  filter_outlier = TRUE,
  sum_max = 20000,
  sum_min = 800,
  gene_max = 7000,
  gene_min = 500,
  mito_percent = 10,
  filter_gene = TRUE,
  min_counts = 5,
  batch = sce$batch
)
sce <- plotQC(sce, file = "plot/naive_only/Quality_Control_after.pdf", width = 20, height = 5)
saveRDS(sce, file = "object/naive_only/allcells_filtered.rds")

sce <- scranNormalize(sce, graph.fun = "louvain")
temp <- exclude_ig(sce)
temp <- dimensionReduction(temp, n = 2000, rank = 50)
temp <- merging(
  temp, batch = "batch", dim = 50,
  method = c("FastMNN", "Harmony"), cores = 8
)

for (rd in c("PCA", "UMAP", "batch_corrected", "UMAP_corrected", "HARMONY", "UMAP_HARMONY")) {
  if (rd %in% reducedDimNames(temp)) reducedDim(sce, rd) <- reducedDim(temp, rd)
}
saveRDS(sce, file = "object/naive_only/allcells_reduced.rds")

sce$cluster_MNN <- clusterCells(
  sce,
  use.dimred = "batch_corrected",
  BLUSPARAM = NNGraphParam(k = 18, cluster.fun = "louvain", BPPARAM = MulticoreParam(8))
)
sce$cluster_HARMONY <- clusterCells(
  sce,
  use.dimred = "HARMONY",
  BLUSPARAM = NNGraphParam(k = 18, cluster.fun = "louvain", BPPARAM = MulticoreParam(8))
)
saveRDS(sce, file = "object/naive_only/allcells_clustered.rds")

pdf("plot/naive_only/allcells_UMAP.pdf", width = 7, height = 6)
reducedPlot(sce, var = "batch", reduction.use = "UMAP_corrected", size = 0.3,
            do.raster = TRUE, do.label = FALSE, main = "Samples")
reducedPlot(sce, var = "cluster_MNN", reduction.use = "UMAP_corrected", size = 0.3,
            do.raster = TRUE, do.label = TRUE, main = "Clusters")
dev.off()

# Remove non-naive clusters (from marker inspection)
naive <- sce[, !sce$cluster_MNN %in% c(6, 8)]
saveRDS(naive, file = "object/naive_only/naive_raw.rds")


# -----------------------------------------------------------------------------
# 2. Re-embed and cluster naive B cells
# -----------------------------------------------------------------------------
sce <- readRDS("object/naive_only/naive_raw.rds")
sce <- scranNormalize(sce, graph.fun = "louvain")
temp <- exclude_ig(sce)
temp <- dimensionReduction(temp, n = 2000, rank = 50)
temp <- merging(
  temp, batch = "batch", dim = 50,
  method = c("FastMNN", "Harmony"), cores = 8
)
for (rd in c("PCA", "UMAP", "batch_corrected", "UMAP_corrected", "HARMONY", "UMAP_HARMONY")) {
  if (rd %in% reducedDimNames(temp)) reducedDim(sce, rd) <- reducedDim(temp, rd)
}
saveRDS(sce, file = "object/naive_only/naive_reduced.rds")

k_range <- seq(10, 50, by = 5)
sce <- clustering(
  sce, k_range, method = "louvain",
  dimred = "batch_corrected", min_cl = 3, cores = 8, top = length(k_range)
)
if (is.list(sce) && !is.null(sce$sce)) sce <- sce$sce
saveRDS(sce, file = "object/naive_only/naive_clustered.rds")

sce <- plotCluster(
  sce, dimred = "UMAP_corrected",
  file = "plot/naive_only/naive_UMAP_clusters.pdf",
  width = 7, height = 6, label = TRUE
)


# -----------------------------------------------------------------------------
# 3. Annotate cell types
# -----------------------------------------------------------------------------
sce$cellcluster <- sce$k.50_cluster.fun.louvain
sce$celltype <- NA_character_
sce$celltype[sce$cellcluster %in% c(1, 2, 4)] <- "aNaive1"
sce$celltype[sce$cellcluster %in% c(5)] <- "rNaive"
sce$celltype[sce$cellcluster %in% c(3)] <- "aNaive2"
sce$celltype[sce$cellcluster %in% c(6)] <- "aNaive3"
sce$celltype <- factor(
  sce$celltype, levels = c("rNaive", "aNaive1", "aNaive2", "aNaive3")
)

# Keep CD27-negative cells
if ("CD27" %in% rownames(sce)) {
  cd27 <- as.numeric(assay(sce, "logcounts")["CD27", ])
  sce <- sce[, cd27 <= 0]
}

sce$disease <- ifelse(grepl("RA", sce$batch), "RA", "HC")
sce$disease <- factor(sce$disease, levels = c("HC", "RA"))
saveRDS(sce, file = "object/naive_only/naive_annotated.rds")

pdf("plot/naive_only/naive_UMAP_celltype.pdf", width = 7, height = 6)
print(
  reducedPlot(
    sce, var = "celltype", reduction.use = "UMAP_corrected",
    size = 0.3, do.raster = TRUE, do.label = FALSE
  ) + scale_color_manual(values = celltype_colors)
)
print(
  reducedPlot(
    sce, var = "disease", reduction.use = "UMAP_corrected",
    size = 0.3, do.raster = TRUE, do.label = FALSE
  )
)
dev.off()

pdf("plot/naive_only/naive_barplot_celltype.pdf", width = 3, height = 6)
print(
  dittoBarPlot(sce, "celltype", group.by = "disease") +
    scale_fill_manual(values = celltype_colors)
)
dev.off()

marker_genes <- c(
  "TCL1A", "CD27", "IGHM", "IGHD", "CXCR5", "CD69", "JUN",
  "TBX21", "HLA-DRA", "CR2", "ITGAX"
)
intergene <- intersect(marker_genes, rownames(sce))
pdf("plot/naive_only/naive_dotplot_celltype.pdf", width = 6, height = 4)
print(dittoDotPlot(sce, vars = intergene, group.by = "celltype"))
dev.off()

message("Done. Annotated naive object: object/naive_only/naive_annotated.rds")
message("Done. Plots: plot/naive_only/")
