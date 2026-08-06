#!/usr/bin/env Rscript
# =============================================================================
# scrna_processing.R
#
# scRNA-seq processing pipeline for RA B-cell analysis:
#   1. Load 10X data
#   2. Quality control
#   3. Normalization, dimension reduction, and batch integration
#   4. UMAP
#   5. Clustering and cell-type annotation
#   6. Trajectory analysis (diffusion map + slingshot)
#
# Expected directory layout (relative to project root):
#   data/rna_seq/<sample>/sample_filtered_feature_bc_matrix/
#   object/   (intermediate Seurat/SCE objects)
#   plot/     (QC and UMAP figures)
#   tables/   (cluster markers)
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
  library(batchelor)
  library(patchwork)
  library(dittoSeq)
  library(destiny)
  library(slingshot)
  library(ggplot2)
  library(dplyr)
})

# Project root: set this to the repository root before running
# setwd("/path/to/aNav_code")

dir.create("object", showWarnings = FALSE, recursive = TRUE)
dir.create("plot", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
sample <- list.files("data/rna_seq")
directory <- paste0("data/rna_seq/", sample, "/sample_filtered_feature_bc_matrix")
sample_sheet <- data.frame(sample = sample, directory = directory)

sce <- loading_sce(directory = sample_sheet$directory, sample = sample_sheet$sample)
saveRDS(sce, file = "object/allcells_raw.rds")


# -----------------------------------------------------------------------------
# 2. Quality control
# -----------------------------------------------------------------------------
sce <- plotQC(sce, file = "plot/Quality_Control_before.pdf", width = 7, height = 5)

sce <- qualityfilter(
  sce,
  filter_outlier = TRUE,
  sum_max = 20000,
  sum_min = 1000,
  gene_max = 6000,
  gene_min = 500,
  mito_percent = 10,
  filter_gene = TRUE,
  min_counts = 5,
  batch = sce$batch
)

sce <- plotQC(sce, file = "plot/Quality_Control_after.pdf", width = 7, height = 5)
saveRDS(sce, file = "object/allcells_filtered.rds")


# -----------------------------------------------------------------------------
# 3. Normalization, dimension reduction, and batch integration
# -----------------------------------------------------------------------------
sce <- scranNormalize(sce, graph.fun = "louvain")

# Exclude Ig variable/constant genes from clustering features
ig_filter <- c("IGLV", "IGHV", "IGKV", "IGLC")
ig_index <- unlist(lapply(ig_filter, function(x) grep(x, rownames(sce))))
ig_genes <- rownames(sce)[ig_index]

temp <- sce[!rownames(sce) %in% ig_genes, ]
temp <- dimensionReduction(temp, n = 2000, rank = 50)

# fastMNN batch correction
temp_corrected <- fastMNN(temp, batch = temp$batch, d = 50, BPPARAM = MulticoreParam())
reducedDim(temp, "batch_corrected") <- reducedDim(temp_corrected, "corrected")
temp <- runUMAP(temp, dimred = "batch_corrected", name = "UMAP_corrected",
                BPPARAM = MulticoreParam())

# Harmony batch correction
temp <- RunHarmony(temp, "batch", ncores = 9, verbose = TRUE)
temp <- runUMAP(temp, dimred = "HARMONY", name = "UMAP_HARMONY",
                BPPARAM = MulticoreParam())

# Transfer embeddings back to full gene set
for (rd in c("PCA", "UMAP", "batch_corrected", "UMAP_corrected", "HARMONY", "UMAP_HARMONY")) {
  reducedDim(sce, rd) <- reducedDim(temp, rd)
}
saveRDS(sce, file = "object/allcells_reduced.rds")


# -----------------------------------------------------------------------------
# 4. UMAP (batch and marker visualization)
# -----------------------------------------------------------------------------
pdf("plot/allcells_UMAP_batch.pdf", width = 8, height = 8)
reducedPlot(sce, var = "batch", reduction.use = "UMAP", size = 0.5,
            do.raster = TRUE, main = "UMAP_original")
reducedPlot(sce, var = "batch", reduction.use = "UMAP_HARMONY", size = 0.5,
            do.raster = TRUE, main = "UMAP_HARMONY")
reducedPlot(sce, var = "batch", reduction.use = "UMAP_corrected", size = 0.5,
            do.raster = TRUE, main = "UMAP_FastMNN")
dev.off()

immune_genes <- c(
  "CD19", "IGHM", "IGHD", "CD27", "CD38", "CD24", "CR2", "CD3D", "CD4", "CD8A",
  "JCHAIN", "MS4A1", "CD3E", "CD14", "XBP1", "IGHG1", "IGHG3", "FCER1A",
  "FCGR3A", "LYZ"
)
intergene <- intersect(immune_genes, rownames(sce))

p <- plot_density(sce, intergene, size = 0.5, reduction = "UMAP_corrected")
pdf("plot/allcells_immuneGenes_density_FastMNN.pdf", width = 30, height = 16)
print(p + plot_layout(ncol = 6))
dev.off()


# -----------------------------------------------------------------------------
# 5. Clustering (all cells) and B-cell filtering
# -----------------------------------------------------------------------------
sce$cluster_HARMONY <- clusterCells(
  sce,
  use.dimred = "HARMONY",
  BLUSPARAM = SNNGraphParam(k = 16, cluster.fun = "louvain")
)
sce$cluster_batchCorrect <- clusterCells(
  sce,
  use.dimred = "batch_corrected",
  BLUSPARAM = SNNGraphParam(k = 30, cluster.fun = "louvain")
)
saveRDS(sce, file = "object/allcells_clustered.rds")

pdf("plot/allcells_UMAP_clusters.pdf", width = 8, height = 8)
reducedPlot(sce, var = "cluster_HARMONY", reduction.use = "UMAP_HARMONY",
            size = 0.5, do.raster = TRUE, do.label = TRUE, main = "Clusters (HARMONY)")
reducedPlot(sce, var = "cluster_batchCorrect", reduction.use = "UMAP_corrected",
            size = 0.5, do.raster = TRUE, do.label = TRUE, main = "Clusters (FastMNN)")
dev.off()

# Remove non-B clusters identified from marker inspection (clusters 4 and 9)
sce <- sce[, !sce$cluster_batchCorrect %in% c(4, 9)]
saveRDS(sce, file = "object/bcells_cleaned.rds")


# -----------------------------------------------------------------------------
# 6. Re-process B cells: normalize, reduce, integrate, UMAP
# -----------------------------------------------------------------------------
sce <- scranNormalize(sce, graph.fun = "louvain")

ig_index <- unlist(lapply(ig_filter, function(x) grep(x, rownames(sce))))
ig_genes <- rownames(sce)[ig_index]
temp <- sce[!rownames(sce) %in% ig_genes, ]

temp <- dimensionReduction(temp, n = 2000, rank = 50)
temp_corrected <- fastMNN(temp, batch = temp$batch, d = 50, BPPARAM = MulticoreParam())
reducedDim(temp, "batch_corrected") <- reducedDim(temp_corrected, "corrected")
temp <- runUMAP(temp, dimred = "batch_corrected", name = "UMAP_corrected",
                BPPARAM = MulticoreParam())

for (rd in c("PCA", "UMAP", "batch_corrected", "UMAP_corrected")) {
  reducedDim(sce, rd) <- reducedDim(temp, rd)
}
saveRDS(sce, file = "object/bcells_reduced.rds")

pdf("plot/bcells_UMAP_batch.pdf", width = 8, height = 8)
reducedPlot(sce, var = "batch", reduction.use = "UMAP", size = 0.5,
            do.raster = TRUE, main = "UMAP_original")
reducedPlot(sce, var = "batch", reduction.use = "UMAP_corrected", size = 0.5,
            do.raster = TRUE, main = "UMAP_FastMNN")
dev.off()

bcell_genes <- c(
  "CD38", "MZB1", "XBP1", "IGLL5", "PDIA6", "IGHG3", "LILRA4",
  "FCRL3", "FGR", "RHOB", "IFITM1", "IFIT2", "IFIT3", "CRIP2",
  "CRIP1", "S100A10", "ITGB1", "HOPX", "JCHAIN", "CCDC50", "COTL1",
  "IGHG1", "IGHA1", "CD40", "CD83", "CD27", "IGHM", "IGHD", "TCL1A",
  "MEF2C", "PLPP5", "FCER2", "CXCR4", "IL4R", "CXCR5", "CD24", "ITGAM",
  "CD1C", "AICDA", "BCL6", "ITGAX", "LAMP1", "TBX21", "CR2", "FCRL5", "CD19"
)
intergene <- intersect(bcell_genes, rownames(sce))
p <- plot_density(sce, intergene, size = 0.5, reduction = "UMAP_corrected")
pdf("plot/bcells_markerGenes_density_FastMNN.pdf", width = 30, height = 36)
print(p + plot_layout(ncol = 6))
dev.off()


# -----------------------------------------------------------------------------
# 7. B-cell clustering and marker detection
# -----------------------------------------------------------------------------
k_range <- seq(10, 50, by = 2)
sce <- clustering(
  sce,
  k_range,
  method = "louvain",
  dimred = "batch_corrected",
  min_cl = 5,
  cores = 8,
  top = 10
)
saveRDS(sce, file = "object/bcells_clustered.rds")

sce <- plotCluster(
  sce,
  dimred = "UMAP_corrected",
  file = "plot/bcells_UMAP_clusters.pdf",
  width = 8,
  height = 8,
  label = TRUE
)

# Cluster resolution used for annotation (k = 36, louvain)
cluster_col <- "k.36_cluster.fun.louvain"
markers <- clusterMarker(sce, top = 10, sce[[cluster_col]])
write.csv(markers, file = "tables/cluster_markers.csv", row.names = FALSE)

selected_markers <- c(
  "ARHGAP15", "AFF3", "FCHSD2", "EBF1", "ANKRD44",
  "ARHGAP24", "ARID1B", "BACH2", "CAMK1D", "PRKCB"
)
intergene <- intersect(c(bcell_genes, selected_markers), rownames(sce))

pdf("plot/bcells_genes_heatmap.pdf", width = 5, height = 8)
averageHeatmap(
  sce, intergene, sce[[cluster_col]],
  annot.by = c("ids"), show_colnames = TRUE,
  cluster_cols = TRUE, scale = "row", main = "FastMNN Clusters"
)
print(dittoDotPlot(sce, vars = intergene, group.by = cluster_col) +
        coord_flip() + ggtitle("FastMNN Clusters"))
dev.off()


# -----------------------------------------------------------------------------
# 8. Cell-type annotation
# -----------------------------------------------------------------------------
sce$cellcluster <- sce[[cluster_col]]

# Broad (celltype1) and fine (celltype2) labels; cluster order = IDs 1..n
celltype1 <- c(
  "naiveB", "naiveB", "switchedMEM_CXCR5+", "unswitchedMEM",
  "unswitchedMEM_CXCR5+", "Exhausted_ABC", "switchedMEM", "naiveB"
)
celltype2 <- c(
  "aNaive_CXCR5+", "aNaive_CXCR5+_CD38hi", "switchedMEM_CXCR5+", "unswitchedMEM",
  "unswitchedMEM_CXCR5+", "Exhausted_ABC", "switchedMEM", "rNaive"
)

annotation <- data.frame(
  cluster = seq_len(length(unique(sce$cellcluster))),
  celltype1 = celltype1,
  celltype2 = celltype2
)

sce <- addCellType(sce, annotation = annotation)
sce$cluster_type <- paste0(sce$celltype1, "_", sce$cellcluster)

# Final manuscript labels (used in figures)
sce$celltype <- as.character(sce$celltype2)
sce$celltype[sce$celltype2 %in% c("unswitchedMEM_CXCR5+")] <- "TransCXCR5+"
sce$celltype[sce$celltype2 %in% c("unswitchedMEM")] <- "MZ_like"
sce$celltype[sce$celltype2 %in% c("switchedMEM_CXCR5+")] <- "unswitchedMEM_CXCR5+"
sce$celltype <- factor(
  sce$celltype,
  levels = c(
    "rNaive", "aNaive_CXCR5+", "aNaive_CXCR5+_CD38hi",
    "TransCXCR5+", "MZ_like", "unswitchedMEM_CXCR5+", "switchedMEM",
    "Exhausted_ABC"
  )
)

# Sample metadata: condition and timepoint from batch IDs
sce$condition[sce$batch %in% c(
  "PC1A", "PC1B", "PC2A", "PC2B", "PC3A", "PC3B", "PC4A", "PC4B", "PC5A", "PC5B"
)] <- "Convertor"
sce$condition[sce$batch %in% c(
  "PN1A", "PN1B", "PN2A", "PN2B", "PN3A", "PN3B", "PN4A", "PN4B", "PN5A", "PN5B"
)] <- "Non_convertor"

sce$timepoint[sce$batch %in% c(
  "PC1A", "PC2A", "PC3A", "PC4A", "PC5A", "PN1A", "PN2A", "PN3A", "PN4A", "PN5A"
)] <- "Before"
sce$timepoint[sce$batch %in% c(
  "PC1B", "PC2B", "PC3B", "PC4B", "PC5B", "PN1B", "PN2B", "PN3B", "PN4B", "PN5B"
)] <- "After"

sce$condition_time <- paste0(sce$condition, "_", sce$timepoint)
sce$condition <- factor(sce$condition, levels = c("Non_convertor", "Convertor"))
sce$timepoint <- factor(sce$timepoint, levels = c("Before", "After"))
sce$condition_time <- factor(
  sce$condition_time,
  levels = c(
    "Non_convertor_Before", "Non_convertor_After",
    "Convertor_Before", "Convertor_After"
  )
)

saveRDS(sce, file = "object/bcells_annotated.rds")

pdf("plot/bcells_UMAP_celltype.pdf", width = 8, height = 8)
reducedPlot(sce, var = "celltype1", reduction.use = "UMAP_corrected",
            size = 0.5, do.raster = TRUE, main = "UMAP Celltype (broad)")
reducedPlot(sce, var = "celltype2", reduction.use = "UMAP_corrected",
            size = 0.5, do.raster = TRUE, main = "UMAP Celltype (fine)")
reducedPlot(sce, var = "celltype", reduction.use = "UMAP_corrected",
            size = 0.5, do.raster = TRUE, main = "UMAP Celltype (final)")
dev.off()


# -----------------------------------------------------------------------------
# 9. Trajectory analysis (diffusion map + slingshot)
# -----------------------------------------------------------------------------
# Diffusion map on batch-corrected embedding
mat <- as.matrix(reducedDim(sce, "batch_corrected"))
mat_diff <- DiffusionMap(mat, verbose = TRUE, n_pcs = NULL, distance = "euclidean")
reducedDim(sce, "DM") <- mat_diff@eigenvectors

pdf("plot/bcells_DM_celltype.pdf", width = 8, height = 8)
reducedPlot(sce, var = "celltype", reduction.use = "DM", size = 0.5,
            do.raster = TRUE, do.label = TRUE, main = "Diffusion map")
dev.off()

# Slingshot lineages
sce <- slingshot(
  sce,
  clusterLabels = "celltype",
  reducedDim = "DM",
  approx_points = 150
)
slsCrv <- slingCurves(sce, as.df = TRUE)

# Average pseudotime across lineages
pt_cols <- grep("^slingPseudotime_", colnames(colData(sce)), value = TRUE)
pt_mat <- as.matrix(colData(sce)[, pt_cols, drop = FALSE])
sce$dpt_average <- rowMeans(pt_mat, na.rm = TRUE)

saveRDS(sce, file = "object/bcells_trajectory.rds")

# Trajectory curves on diffusion map
pdf("plot/bcells_DM_trajectory.pdf", width = 8, height = 8)
gra <- reducedPlot(sce, var = "celltype", reduction.use = "DM", size = 0.5,
                   do.raster = TRUE, do.label = TRUE, main = "Slingshot lineages")
print(
  gra + theme_classic() +
    geom_path(
      data = slsCrv %>% arrange(Order),
      aes(DC1, DC2, group = Lineage),
      linewidth = 1.2
    )
)

gra <- reducedPlot(sce, var = "dpt_average", reduction.use = "DM", size = 0.5,
                   do.raster = TRUE, do.label = FALSE, main = "Pseudotime")
print(
  gra + scale_color_continuous(type = "viridis") + theme_classic() +
    geom_path(
      data = slsCrv %>% arrange(Order),
      aes(DC1, DC2, group = Lineage),
      linewidth = 1.2
    )
)
dev.off()

# Pseudotime on UMAP
pdf("plot/bcells_UMAP_pseudotime.pdf", width = 8, height = 8)
print(
  reducedPlot(sce, var = "dpt_average", reduction.use = "UMAP_corrected",
              size = 0.5, do.raster = TRUE, do.label = FALSE) +
    scale_color_continuous(type = "viridis") + ggtitle("Pseudotime")
)
dev.off()

message("Done. Annotated object: object/bcells_annotated.rds")
message("Done. Trajectory object: object/bcells_trajectory.rds")
