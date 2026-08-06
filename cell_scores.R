#!/usr/bin/env Rscript
# =============================================================================
# cell_scores.R
#
# Pathway / signature scoring for RA B cells:
#   1. Build KEGG pathway gene signatures
#   2. UCell scoring and kNN smoothing
#   3. Score visualization by cell type and condition
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
  library(UCell)
  library(BiocParallel)
  library(limma)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(SingleCellExperiment)
  library(dittoSeq)
  library(ggplot2)
})

# Project root: set this to the repository root before running
# setwd("/path/to/aNav_code")

dir.create("plot/scores", showWarnings = FALSE, recursive = TRUE)
dir.create("tables/scores", showWarnings = FALSE, recursive = TRUE)
dir.create("object", showWarnings = FALSE, recursive = TRUE)

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)


# -----------------------------------------------------------------------------
# 1. Load annotated B cells
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


# -----------------------------------------------------------------------------
# 2. Build KEGG pathway signatures
# -----------------------------------------------------------------------------
pwy_names <- c(
  "Rheumatoid_arthritis", "MAPK", "NFkappaB", "PDL1_PD1",
  "Antigen_processing", "Protein_endoplasmic_reticulum",
  "EBV_infection", "Tcell_leukemia_virus1"
)
pwy_ids <- c(
  "hsa05323", "hsa04010", "hsa04064", "hsa05235",
  "hsa04612", "hsa04141", "hsa05169", "hsa05166"
)

kegg <- getGeneKEGGLinks(species.KEGG = "hsa")
signatures <- list()

for (i in seq_along(pwy_ids)) {
  pwy_ref <- kegg[kegg$PathwayID %in% pwy_ids[i], ]
  pwy_ref$genes <- mapIds(
    org.Hs.eg.db,
    keys = pwy_ref$GeneID,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  intergene <- intersect(pwy_ref$genes, rownames(sce))
  signatures[[pwy_names[i]]] <- intergene
  message(pwy_names[i], " has ", length(intergene), " genes")
}

# Autoimmune signature: genes shared across autoimmune-related KEGG pathways
auto_ids <- c("hsa05323", "hsa05322", "hsa04940", "hsa05321", "hsa05320")
pwy_ref <- kegg[kegg$PathwayID %in% auto_ids[1], ]
pwy_ref$genes <- mapIds(
  org.Hs.eg.db,
  keys = pwy_ref$GeneID,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
auto_genes <- pwy_ref$genes

for (i in 2:length(auto_ids)) {
  pwy_ref <- kegg[kegg$PathwayID %in% auto_ids[i], ]
  pwy_ref$genes <- mapIds(
    org.Hs.eg.db,
    keys = pwy_ref$GeneID,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  auto_genes <- intersect(auto_genes, pwy_ref$genes)
}
signatures[["Autoimmune"]] <- intersect(auto_genes, rownames(sce))
message("Autoimmune has ", length(signatures[["Autoimmune"]]), " genes")

# Save signature gene lists
sig_df <- do.call(rbind, lapply(names(signatures), function(n) {
  data.frame(signature = n, gene = signatures[[n]], stringsAsFactors = FALSE)
}))
write.csv(sig_df, file = "tables/scores/signature_genes.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 3. UCell scoring and kNN smoothing
# -----------------------------------------------------------------------------
sce <- ScoreSignatures_UCell(
  sce,
  features = signatures,
  assay = "counts",
  name = NULL,
  BPPARAM = MulticoreParam()
)
sce <- SmoothKNN(
  sce,
  signature.names = names(signatures),
  reduction = "batch_corrected",
  BPPARAM = MulticoreParam()
)

# Add smoothed scores to colData for easy plotting
ucell_knn <- altExp(sce, "UCell_kNN")
ucell_knn <- data.frame(t(assay(ucell_knn, "UCell_kNN")))
# Align column names (UCell may append suffixes)
colnames(ucell_knn) <- gsub("_UCell$", "", colnames(ucell_knn))
colData(sce) <- cbind(colData(sce), ucell_knn)

# Convenience alias for EBV score used in figures
if ("EBV_infection_kNN" %in% colnames(colData(sce))) {
  sce$EBV_score_KNN <- sce$EBV_infection_kNN
} else if ("EBV_infection" %in% colnames(colData(sce))) {
  sce$EBV_score_KNN <- sce$EBV_infection
}

saveRDS(sce, file = "object/bcells_scored.rds")


# -----------------------------------------------------------------------------
# 4. Score plots by cell type
# -----------------------------------------------------------------------------
score_cols <- intersect(
  c(paste0(names(signatures), "_kNN"), names(signatures), "EBV_score_KNN"),
  colnames(colData(sce))
)
# Prefer smoothed (_kNN) scores when available
score_cols <- unique(score_cols)
score_plot <- grep("_kNN$|EBV_score_KNN", score_cols, value = TRUE)
if (length(score_plot) == 0) score_plot <- score_cols

pdf("plot/scores/scores_by_celltype.pdf", height = 6, width = 10)
for (p in score_plot) {
  print(
    dittoPlot(
      sce, p, group.by = "celltype",
      plots = c("jitter", "boxplot"), do.raster = TRUE
    ) + ggtitle(p)
  )
}
dev.off()


# -----------------------------------------------------------------------------
# 5. Score plots by condition (Before samples)
# -----------------------------------------------------------------------------
pdf("plot/scores/scores_by_condition_Before.pdf", height = 6, width = 10)
temp <- sce[, sce$timepoint %in% "Before"]
for (p in score_plot) {
  print(
    dittoPlot(
      temp, p, group.by = "condition",
      plots = c("jitter", "boxplot"), do.raster = TRUE, legend.show = FALSE
    ) + ggtitle(paste0(p, " (Before)"))
  )
}
dev.off()

pdf("plot/scores/scores_by_condition_time.pdf", height = 6, width = 12)
for (p in score_plot) {
  print(
    dittoPlot(
      sce, p, group.by = "condition_time",
      plots = c("jitter", "boxplot"), do.raster = TRUE, legend.show = FALSE
    ) + ggtitle(p)
  )
}
dev.off()


# -----------------------------------------------------------------------------
# 6. EBV score by cell type within each condition_time
# -----------------------------------------------------------------------------
if ("EBV_score_KNN" %in% colnames(colData(sce)) ||
    "EBV_infection_kNN" %in% colnames(colData(sce))) {
  ebv_col <- if ("EBV_score_KNN" %in% colnames(colData(sce))) {
    "EBV_score_KNN"
  } else {
    "EBV_infection_kNN"
  }

  pdf("plot/scores/EBV_score_by_celltype.pdf", height = 6, width = 8)
  print(
    dittoPlot(
      sce, ebv_col, group.by = "celltype",
      plots = c("jitter", "boxplot"), do.raster = TRUE
    )
  )
  dev.off()

  celltypes <- unique(as.character(sce$celltype))
  for (type in celltypes) {
    temp <- sce[, sce$celltype %in% type]
    out <- paste0("plot/scores/EBV_score_", safe_name(type), ".pdf")
    pdf(out, height = 6, width = 8)
    print(
      dittoPlot(
        temp, ebv_col, group.by = "condition_time",
        plots = c("jitter", "boxplot"), do.raster = TRUE
      ) + ggtitle(type)
    )
    dev.off()
  }
}

message("Done. Scored object: object/bcells_scored.rds")
message("Done. Score plots: plot/scores/")
message("Done. Signature genes: tables/scores/signature_genes.csv")
