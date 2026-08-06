#!/usr/bin/env Rscript
# =============================================================================
# proportion_PCA.R
#
# Cell-type proportion and pseudobulk PCA analysis for RA B cells:
#   1. Cell-type composition by condition
#   2. Per-sample proportion boxplots (Before / After)
#   3. Pseudobulk PCA by cell type (RUV-corrected)
#
# Custom packages developed for this study:
#   bioscrna  Single-cell RNA-seq analysis toolkit built on Bioconductor
#             (data loading, QC, dimension reduction, clustering, and plotting
#             for SingleCellExperiment objects).
#             https://github.com/zhangzmr/bioscrna
#             remotes::install_github("zhangzmr/bioscrna")
#
#   scgtest   Differential expression toolkit for single-cell and pseudobulk
#             RNA-seq (edgeR, RUV correction, pathway enrichment, and volcano
#             plots).
#             https://github.com/zhangzmr/scgtest
#             remotes::install_github("zhangzmr/scgtest")
# =============================================================================

suppressPackageStartupMessages({
  library(bioscrna)
  library(scgtest)
  library(dittoSeq)
  library(ggplot2)
  library(ggforce)
  library(stringr)
  library(edgeR)
  library(limma)
  library(MatrixGenerics)
})

# Project root: set this to the repository root before running
# setwd("/path/to/aNav_code")

dir.create("tables/proportion", showWarnings = FALSE, recursive = TRUE)
dir.create("plot/proportion", showWarnings = FALSE, recursive = TRUE)
dir.create("plot/PCA", showWarnings = FALSE, recursive = TRUE)
dir.create("object", showWarnings = FALSE, recursive = TRUE)

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

celltype_colors <- c(
  "#A469BD", "#D5695D", "#F5B041", "#FFBCA7",
  "#5DADE2", "#484F98", "#8A7067", "#52BE80"
)
pca_colors <- c(
  "dodgerblue", "blue", "green", "darkgreen", "deeppink", "darkorchid"
)


# -----------------------------------------------------------------------------
# 1. Load annotated B cells and set metadata
# -----------------------------------------------------------------------------
sce <- readRDS("object/bcells_annotated.rds")

batch_id <- gsub("^B_", "", as.character(sce$batch))

sce$condition[batch_id %in% c(
  "PC1A", "PC1B", "PC2A", "PC2B", "PC3A", "PC3B", "PC4A", "PC4B", "PC5A", "PC5B"
)] <- "Convertor"
sce$condition[batch_id %in% c(
  "PN1A", "PN1B", "PN2A", "PN2B", "PN3A", "PN3B", "PN4A", "PN4B", "PN5A", "PN5B"
)] <- "Non_convertor"
sce$condition[batch_id %in% c(
  "NN1A", "NN1B", "NN2A", "NN2B", "NN3A", "NN3B"
)] <- "Negative_control"

sce$timepoint[batch_id %in% c(
  "NN1A", "NN2A", "NN3A", "PC1A", "PC2A", "PC3A", "PC4A", "PC5A",
  "PN1A", "PN2A", "PN3A", "PN4A", "PN5A"
)] <- "Before"
sce$timepoint[batch_id %in% c(
  "NN1B", "NN2B", "NN3B", "PC1B", "PC2B", "PC3B", "PC4B", "PC5B",
  "PN1B", "PN2B", "PN3B", "PN4B", "PN5B"
)] <- "After"

sce$condition_time <- paste0(sce$condition, "_", sce$timepoint)
sce$condition <- factor(
  sce$condition,
  levels = c("Negative_control", "Non_convertor", "Convertor")
)
sce$timepoint <- factor(sce$timepoint, levels = c("Before", "After"))
sce$condition_time <- factor(
  sce$condition_time,
  levels = c(
    "Negative_control_Before", "Negative_control_After",
    "Non_convertor_Before", "Non_convertor_After",
    "Convertor_Before", "Convertor_After"
  )
)

sce$celltype <- factor(
  sce$celltype,
  levels = c(
    "rNaive", "aNaive_CXCR5+", "aNaive_CXCR5+_CD38hi",
    "TransCXCR5+", "MZ_like", "unswitchedMEM_CXCR5+", "switchedMEM",
    "Exhausted_ABC"
  )
)


# -----------------------------------------------------------------------------
# 2. Cell-type composition by condition
# -----------------------------------------------------------------------------
pdf("plot/proportion/barplot_by_condition.pdf", width = 6, height = 6)
print(
  dittoBarPlot(
    sce, "celltype", group.by = "condition",
    retain.factor.levels = TRUE, color.panel = celltype_colors
  ) +
    theme(
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12),
      title = element_blank()
    )
)
dev.off()

pdf("plot/proportion/barplot_by_condition_time.pdf", width = 8, height = 6)
print(
  dittoBarPlot(
    sce, "celltype", group.by = "condition_time",
    retain.factor.levels = TRUE, color.panel = celltype_colors
  ) +
    theme(
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12),
      title = element_blank()
    )
)
dev.off()


# -----------------------------------------------------------------------------
# 3. Per-sample proportions and boxplots
# -----------------------------------------------------------------------------
sce$group <- paste0(sce$condition_time, "--", sce$batch)
plot_out <- dittoBarPlot(
  sce, "celltype", group.by = "group",
  retain.factor.levels = TRUE, data.out = TRUE
)
plot_data <- plot_out$data
split_meta <- data.frame(str_split_fixed(plot_data$grouping, pattern = "--", n = 2))
colnames(split_meta) <- c("condition_time", "sample")
plot_data <- cbind(plot_data, split_meta)
plot_data$condition_time <- factor(
  plot_data$condition_time,
  levels = levels(sce$condition_time)
)

write.csv(
  plot_data,
  file = "tables/proportion/Bcell_proportion_by_sample.csv",
  row.names = FALSE
)

pdf("plot/proportion/boxplot_by_condition_time.pdf", width = 17, height = 4)
print(
  ggplot(
    plot_data,
    aes(x = label, y = percent, fill = condition_time)
  ) +
    theme_classic() +
    geom_boxplot() +
    geom_jitter(
      shape = 21, position = position_jitterdodge(),
      size = 1.5, stroke = 1
    ) +
    theme(
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    ) +
    ylab("Percentage")
)
dev.off()


# -----------------------------------------------------------------------------
# 4. Pseudobulk PCA by cell type (RUV-corrected)
# -----------------------------------------------------------------------------
add_bulk_metadata <- function(bulk) {
  metadata <- bulk$metadata
  colnames(metadata) <- c("sample", "nCells")
  sid <- gsub("^B_", "", as.character(metadata$sample))

  metadata$condition[sid %in% c(
    "PC1A", "PC1B", "PC2A", "PC2B", "PC3A", "PC3B", "PC4A", "PC4B", "PC5A", "PC5B"
  )] <- "Convertor"
  metadata$condition[sid %in% c(
    "PN1A", "PN1B", "PN2A", "PN2B", "PN3A", "PN3B", "PN4A", "PN4B", "PN5A", "PN5B"
  )] <- "Non_convertor"
  metadata$condition[sid %in% c(
    "NN1A", "NN1B", "NN2A", "NN2B", "NN3A", "NN3B"
  )] <- "Negative_control"

  metadata$timepoint[sid %in% c(
    "NN1A", "NN2A", "NN3A", "PC1A", "PC2A", "PC3A", "PC4A", "PC5A",
    "PN1A", "PN2A", "PN3A", "PN4A", "PN5A"
  )] <- "Before"
  metadata$timepoint[sid %in% c(
    "NN1B", "NN2B", "NN3B", "PC1B", "PC2B", "PC3B", "PC4B", "PC5B",
    "PN1B", "PN2B", "PN3B", "PN4B", "PN5B"
  )] <- "After"

  metadata$condition_time <- paste0(metadata$condition, "_", metadata$timepoint)
  metadata$condition <- factor(
    metadata$condition,
    levels = c("Negative_control", "Non_convertor", "Convertor")
  )
  metadata$timepoint <- factor(metadata$timepoint, levels = c("Before", "After"))
  metadata$condition_time <- factor(
    metadata$condition_time,
    levels = c(
      "Negative_control_Before", "Negative_control_After",
      "Non_convertor_Before", "Non_convertor_After",
      "Convertor_Before", "Convertor_After"
    )
  )
  bulk$metadata <- metadata
  bulk
}

design_levels <- c(
  "Negative_control_Before", "Negative_control_After",
  "Non_convertor_Before", "Non_convertor_After",
  "Convertor_Before", "Convertor_After"
)

celltypes <- levels(sce$celltype)
pca_list <- list()

for (type in celltypes) {
  scesub <- sce[, sce$celltype %in% type]
  if (ncol(scesub) < 10) next

  bulk <- convertBulkFilter(scesub, group = scesub$batch, filter = FALSE)
  bulk <- add_bulk_metadata(bulk)

  counts <- bulk$counts
  metadata <- bulk$metadata
  rowcounts <- rowSums2(counts)
  keep_genes <- names(rowcounts)[rowcounts >= 5]
  counts <- counts[rownames(counts) %in% keep_genes, , drop = FALSE]

  # Need enough samples for PCA
  if (ncol(counts) < 3) next

  design <- model.matrix(~ 0 + condition_time, data = metadata)
  colnames(design) <- design_levels

  set2 <- RUVcorrection(
    counts = counts, metadata = metadata, design = design,
    group = "condition_time",
    contrasts = "Convertor_Before-Non_convertor_Before",
    k = 1, nGenes = 5000
  )

  dge <- DGEList(counts = counts, samples = metadata, group = metadata$condition_time)
  dge <- calcNormFactors(dge, method = "TMM")
  cpm <- edgeR::cpm(dge, log = FALSE)
  cpm <- log2(cpm + 1)

  removed <- removeBatchEffect(cpm, covariates = set2$W, design = design)
  pca_res <- prcomp(t(removed), scale. = TRUE)
  pca <- pca_res$x

  plot_df <- cbind(as.data.frame(pca[, 1:2]), metadata)
  gra <- ggplot(plot_df, aes(x = PC1, y = PC2, color = condition_time)) +
    geom_point(size = 3) +
    geom_mark_ellipse(aes(color = condition_time), linetype = 2, size = 0.7) +
    scale_color_manual(values = pca_colors) +
    theme_classic() +
    ggtitle(type) +
    theme(
      legend.position = "right",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )

  out <- paste0("plot/PCA/", safe_name(type), "_PCA.pdf")
  pdf(out, height = 5, width = 6)
  print(gra)
  dev.off()

  pca_list[[type]] <- list(
    removed = removed,
    metadata = metadata,
    pca = pca_res
  )
  message("PCA: ", type)
}

saveRDS(pca_list, file = "object/bcells_PCA_list.rds")

message("Done. Proportion plots: plot/proportion/")
message("Done. PCA plots: plot/PCA/")
message("Done. Proportion table: tables/proportion/Bcell_proportion_by_sample.csv")
