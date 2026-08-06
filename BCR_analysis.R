#!/usr/bin/env Rscript
# =============================================================================
# BCR_analysis.R
#
# BCR repertoire analysis for RA B cells:
#   1. Load and clean BCR annotations
#   2. Heavy / light chain cloning, germline, and SHM
#   3. Merge BCR with annotated scRNA-seq
#   4. Clonality classification and SHM visualization
#   5. aNaive-related clonal fate analysis
#   6. Phylogenetic trees (screceptor)
#
# Custom packages developed for this study:
#   bioscrna   Single-cell RNA-seq analysis toolkit built on Bioconductor
#              (data loading, QC, dimension reduction, clustering, and plotting
#              for SingleCellExperiment objects).
#              https://github.com/zhangzmr/bioscrna
#              remotes::install_github("zhangzmr/bioscrna")
#
#   screceptor Process single-cell BCR/TCR repertoire data (loading, cleaning,
#              clonal clustering, germline reconstruction, SHM, merging with
#              gene expression, and phylogenetic tree helpers).
#              https://github.com/zhangzmr/screceptor
#              remotes::install_github("zhangzmr/screceptor")
# =============================================================================

suppressPackageStartupMessages({
  library(bioscrna)
  library(screceptor)
  library(dplyr)
  library(ggplot2)
  library(dittoSeq)
  library(viridis)
  library(ggsci)
})

# Project root: set this to the repository root before running
# setwd("/path/to/aNav_code")

# Path to IMGT / Dowser germline reference (required for germlineBCR)
# references <- readRDS("data/reference/human_vdj/vdj_references.rds")
references <- NULL

dir.create("data/bcr", showWarnings = FALSE, recursive = TRUE)
dir.create("bcr_result", showWarnings = FALSE, recursive = TRUE)
dir.create("plot/bcr", showWarnings = FALSE, recursive = TRUE)
dir.create("tables/bcr", showWarnings = FALSE, recursive = TRUE)
dir.create("tree/data", showWarnings = FALSE, recursive = TRUE)
dir.create("object", showWarnings = FALSE, recursive = TRUE)

add_sample_metadata <- function(df, sample_col = "sample_id") {
  sid <- gsub("^B_", "", as.character(df[[sample_col]]))
  df$condition[sid %in% c(
    "PC1A", "PC1B", "PC2A", "PC2B", "PC3A", "PC3B", "PC4A", "PC4B", "PC5A", "PC5B"
  )] <- "Convertor"
  df$condition[sid %in% c(
    "PN1A", "PN1B", "PN2A", "PN2B", "PN3A", "PN3B", "PN4A", "PN4B", "PN5A", "PN5B"
  )] <- "Non_convertor"
  df$condition[sid %in% c(
    "NN1A", "NN1B", "NN2A", "NN2B", "NN3A", "NN3B"
  )] <- "Negative_control"

  df$timepoint[sid %in% c(
    "NN1A", "NN2A", "NN3A", "PC1A", "PC2A", "PC3A", "PC4A", "PC5A",
    "PN1A", "PN2A", "PN3A", "PN4A", "PN5A"
  )] <- "Before"
  df$timepoint[sid %in% c(
    "NN1B", "NN2B", "NN3B", "PC1B", "PC2B", "PC3B", "PC4B", "PC5B",
    "PN1B", "PN2B", "PN3B", "PN4B", "PN5B"
  )] <- "After"

  df$condition_time <- paste0(df$condition, "_", df$timepoint)
  df$condition <- factor(
    df$condition,
    levels = c("Negative_control", "Non_convertor", "Convertor")
  )
  df$timepoint <- factor(df$timepoint, levels = c("Before", "After"))
  df$condition_time <- factor(
    df$condition_time,
    levels = c(
      "Negative_control_Before", "Negative_control_After",
      "Non_convertor_Before", "Non_convertor_After",
      "Convertor_Before", "Convertor_After"
    )
  )
  df
}

assign_clonality <- function(x, size_col) {
  out <- rep(NA_character_, length(x[[size_col]]))
  sz <- x[[size_col]]
  out[sz == 1] <- "Singleton"
  out[sz >= 2 & sz < 5] <- "Small"
  out[sz >= 5 & sz < 10] <- "Medium"
  out[sz >= 10] <- "Large"
  factor(out, levels = c("Large", "Medium", "Small", "Singleton"))
}


# -----------------------------------------------------------------------------
# 1. Load BCR annotations
# -----------------------------------------------------------------------------
# Expected layout:
#   data/bcr/<sample>/BCR_data_sequences_igblast_db-pass.tsv
sample <- list.files("data/bcr")
directory <- file.path(
  "data/bcr", sample, "BCR_data_sequences_igblast_db-pass.tsv"
)
sample_sheet <- data.frame(sample = sample, directory = directory)
write.csv(sample_sheet, file = "bcr_result/bcr_sheet.csv", row.names = FALSE)

bcrlist <- loadBCR(directory = sample_sheet$directory, sample = sample_sheet$sample)
bcr <- bind_rows(bcrlist)
bcr <- cleanBCR(bcr)
saveRDS(bcr, file = "bcr_result/bcr_cleaned.rds")


# -----------------------------------------------------------------------------
# 2. Clone heavy / light chains, resolve pairs, germline, SHM
# -----------------------------------------------------------------------------
heavy_list <- list()
light_list <- list()

for (s in unique(bcr$sample_id)) {
  tmp <- bcr[bcr$sample_id %in% s, ]

  heavy <- subsetHC(tmp)
  heavy <- cloneHC(heavy, similarity = 0.75)
  heavy_list[[s]] <- heavy

  light <- subsetLC(tmp)
  light <- cloneLC(light, bcrlist_heavy = heavy_list)
  light_list[[s]] <- light

  message("Cloned sample: ", s)
}

combined <- combineHL(heavy_list, light_list, output = "combine")
combined <- subCloneLC(combined, nproc = 6)

if (!is.null(references)) {
  combined <- germlineBCR(combined, reference = references, nproc = 6)
  combined <- hypermuBCR(combined, nproc = 6)
} else {
  warning(
    "Germline reference not set (`references` is NULL). ",
    "Skipping germlineBCR / hypermuBCR. ",
    "Set `references` to the IMGT/Dowser germline object to enable SHM."
  )
}

combined <- add_sample_metadata(combined)
saveRDS(combined, file = "bcr_result/BCR_combined.rds")
write.csv(combined, file = "bcr_result/BCR_combined.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 3. Merge BCR with annotated scRNA-seq
# -----------------------------------------------------------------------------
sce <- readRDS("object/bcells_annotated.rds")

# Ensure batch_barcode exists for merging
if (!"batch_barcode" %in% colnames(colData(sce))) {
  sce$batch_barcode <- paste0(sce$batch, "_", colnames(sce))
}

sce_bcr <- mergeBCRtoGEX(sce, bcr = combined, intersect = TRUE)

gex_vars <- intersect(
  c("celltype", "celltype1", "celltype2", "condition", "timepoint",
    "condition_time", "EBV_score_KNN"),
  colnames(colData(sce))
)
BCRgex <- mergeGEXtoBCR(
  sce, bcr = combined, intersect = FALSE, var = gex_vars
)
BCRgex <- add_sample_metadata(BCRgex)


# -----------------------------------------------------------------------------
# 4. Clonality classification
# -----------------------------------------------------------------------------
if ("light_cloneSizeHL" %in% colnames(colData(sce_bcr))) {
  sce_bcr$clonality <- assign_clonality(colData(sce_bcr), "light_cloneSizeHL")
} else if ("cloneSizeHL" %in% colnames(colData(sce_bcr))) {
  sce_bcr$clonality <- assign_clonality(colData(sce_bcr), "cloneSizeHL")
}

if ("cloneSizeHL" %in% colnames(BCRgex)) {
  BCRgex$clonality <- assign_clonality(BCRgex, "cloneSizeHL")
}

saveRDS(sce_bcr, file = "object/bcells_annotated_bcr.rds")
saveRDS(BCRgex, file = "bcr_result/bcr_combined_annotated.rds")
write.csv(
  BCRgex,
  file = "bcr_result/bcr_combined_annotated.csv",
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 5. SHM and clonality plots
# -----------------------------------------------------------------------------
pdf("plot/bcr/UMAP_SHM.pdf", width = 8, height = 8)
print(
  reducedPlot(
    sce_bcr, var = "celltype", reduction.use = "UMAP_corrected",
    size = 0.5, do.raster = TRUE, do.label = TRUE
  )
)
if ("heavy_SHM" %in% colnames(colData(sce_bcr))) {
  print(
    reducedPlot(
      sce_bcr, var = "heavy_SHM", reduction.use = "UMAP_corrected",
      size = 0.5, do.raster = TRUE
    ) + scale_color_viridis(option = "magma")
  )
}
if ("light_SHM" %in% colnames(colData(sce_bcr))) {
  print(
    reducedPlot(
      sce_bcr, var = "light_SHM", reduction.use = "UMAP_corrected",
      size = 0.5, do.raster = TRUE
    ) + scale_color_viridis(option = "magma")
  )
}
dev.off()

clonality_colors <- c("#982c2c", "#f8984f", "#f3d78a", "#f2f1e6")
if ("clonality" %in% colnames(colData(sce_bcr))) {
  pdf("plot/bcr/barplot_clonality.pdf", width = 6, height = 4)
  print(
    dittoBarPlot(
      sce_bcr, var = "clonality", group.by = "condition_time",
      retain.factor.levels = TRUE, color.panel = clonality_colors
    )
  )
  print(
    dittoBarPlot(
      sce_bcr, var = "clonality", group.by = "celltype",
      retain.factor.levels = TRUE, color.panel = clonality_colors
    )
  )
  dev.off()
}

# SHM by cell type (heavy chain table)
heavy <- BCRgex[BCRgex$locus %in% "IGH", ]
shm_col <- intersect(c("SHM_freq", "SHM", "mu_freq"), colnames(heavy))
if (length(shm_col) > 0 && "celltype" %in% colnames(heavy)) {
  shm_col <- shm_col[1]
  pdf("plot/bcr/boxplot_SHM_by_celltype.pdf", width = 8, height = 5)
  print(
    ggplot(heavy, aes(x = celltype, y = .data[[shm_col]], fill = celltype)) +
      geom_boxplot(outlier.size = 0.3) +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      ) +
      ylab("Heavy-chain SHM") + xlab(NULL)
  )
  dev.off()
}


# -----------------------------------------------------------------------------
# 6. aNaive-related clonal fate
# -----------------------------------------------------------------------------
# Clones containing IGHM+ aNaive_CXCR5+ cells; track related non-naive fates
bcr_table <- BCRgex
bcr_table <- bcr_table[!is.na(bcr_table$celltype), ]

if (all(c("cloneSizeHC", "c_call", "clone_id", "celltype") %in% colnames(bcr_table))) {
  # IGHM+ related clones
  igm <- bcr_table[bcr_table$c_call %in% "IGHM" & bcr_table$cloneSizeHC > 1, ]
  aNav <- igm[igm$celltype %in% "aNaive_CXCR5+", ]
  clone_id <- unique(aNav$clone_id)

  related <- bcr_table[bcr_table$clone_id %in% clone_id, ]
  related$temp <- paste0(related$condition_time, "--", related$celltype)
  related_nn <- related[
    !related$celltype %in% c("aNaive_CXCR5+", "aNaive_CXCR5+_CD38hi", "rNaive"),
  ]
  write.csv(
    related_nn,
    file = "tables/bcr/aNaive_related_IGHM.csv",
    row.names = FALSE
  )

  if (nrow(related_nn) > 0) {
    freq <- as.data.frame(table(related_nn$temp))
    freq$condition <- gsub("--.*", "", freq$Var1)
    freq$celltype <- gsub(".*--", "", freq$Var1)
    freq$Freq <- freq$Freq / ave(freq$Freq, freq$condition, FUN = sum)

    pdf("plot/bcr/aNaive_related_IGHM_fate.pdf", width = 11, height = 6)
    print(
      ggplot(freq, aes(x = condition, y = Freq, fill = celltype)) +
        geom_bar(stat = "identity") +
        scale_fill_bmj() +
        theme_classic() +
        ylab("Proportion") + xlab(NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    )
    dev.off()
  }

  # Class-switched related clones
  switched <- bcr_table[
    bcr_table$cloneSizeHC > 1 &
      bcr_table$locus %in% "IGH" &
      bcr_table$clone_id %in% clone_id &
      !bcr_table$c_call %in% "IGHM" &
      !bcr_table$celltype %in% c("aNaive_CXCR5+", "aNaive_CXCR5+_CD38hi", "rNaive"),
  ]
  write.csv(
    switched,
    file = "tables/bcr/aNaive_related_class_switched.csv",
    row.names = FALSE
  )

  if (nrow(switched) > 0) {
    switched$temp <- paste0(switched$condition_time, "--", switched$celltype)
    freq <- as.data.frame(table(switched$temp))
    freq$condition <- gsub("--.*", "", freq$Var1)
    freq$celltype <- gsub(".*--", "", freq$Var1)
    freq$Freq <- freq$Freq / ave(freq$Freq, freq$condition, FUN = sum)

    pdf("plot/bcr/aNaive_related_class_switched_fate.pdf", width = 11, height = 6)
    print(
      ggplot(freq, aes(x = condition, y = Freq, fill = celltype)) +
        geom_bar(stat = "identity") +
        scale_fill_bmj() +
        theme_classic() +
        ylab("Proportion") + xlab(NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    )
    dev.off()
  }
}

# -----------------------------------------------------------------------------
# 7. Phylogenetic trees (screceptor)
# -----------------------------------------------------------------------------
# Write tree inputs; after external RAxML / IgPhyML, assemble and plot.
bcr_table <- readRDS("bcr_result/bcr_combined_annotated.rds")
if (!"SHM" %in% colnames(bcr_table)) bcr_table$SHM <- bcr_table$SHM_freq

for (s in unique(bcr_table$sample_id)) {
  data <- bcr_table[bcr_table$sample_id %in% s, ]
  cutted <- cutoffBCR(data, clone_size = 2, clonesize_col = "cloneSizeHL")
  out_dir <- file.path("tree/data", s)
  germlineTree(cutted, sample = s, out_file_dir = out_dir)
  smallTreeTable(cutted, sample = s, clone_col = "clone_subgroup_id",
                 out_file_dir = out_dir)
}

# Example: plot one sample after RAxML / IgPhyML outputs are available
# dir <- "tree/tree_files/B_PC1A"
# tip_bcr <- addClonality(readRDS(file.path(dir, "bcr_table.rds")))
# tree <- assembleBigTree(file.path(dir, "RAxML_bestTree.germline_tree"),
#                         dir = dir, tip_ids = tip_bcr$sequence_id)
# layers <- list(
#   list(column = "celltype", type = "tile", scale = "manual",
#        values = screceptor_celltype_colors, name = "Celltype"),
#   list(column = "SHM", type = "tile", scale = "viridis",
#        viridis_option = "C", limits = c(0, 0.25), name = "SHM"),
#   list(column = "c_call", type = "star", scale = "manual",
#        values = screceptor_isotype_colors,
#        shapes = screceptor_isotype_shapes, name = "Isotypes")
# )
# plotBigTree(tree, tip_bcr, layers = layers)
# plotSmallTrees(file.path(dir, "small_tree"), tip_bcr, layers = layers)

message("Done. BCR combined: bcr_result/BCR_combined.rds")
message("Done. Annotated BCR table: bcr_result/bcr_combined_annotated.rds")
message("Done. GEX+BCR object: object/bcells_annotated_bcr.rds")
message("Done. Tree inputs: tree/data/")
message("Done. Plots: plot/bcr/")
