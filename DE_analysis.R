#!/usr/bin/env Rscript
# =============================================================================
# DE_analysis.R
#
# Differential expression analysis for RA B cells:
#   1. Pseudobulk aggregation by cell type and sample
#   2. RUV-corrected edgeR DE (condition x time contrasts)
#   3. Pathway enrichment (GSEA)
#   4. Volcano plots
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
  library(org.Hs.eg.db)
  library(clusterProfiler)
})

# Project root: set this to the repository root before running
# setwd("/path/to/aNav_code")

dir.create("tables/DE", showWarnings = FALSE, recursive = TRUE)
dir.create("plot/DE/volcano", showWarnings = FALSE, recursive = TRUE)
dir.create("object", showWarnings = FALSE, recursive = TRUE)

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)


# -----------------------------------------------------------------------------
# 1. Load annotated B cells
# -----------------------------------------------------------------------------
sce <- readRDS("object/bcells_annotated.rds")

# Ensure condition / timepoint metadata from batch IDs
# Supports batch names with or without a "B_" prefix
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

celltypes <- levels(sce$celltype)
if (is.null(celltypes)) celltypes <- unique(as.character(sce$celltype))


# -----------------------------------------------------------------------------
# 2. Pseudobulk by cell type
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

cluster_list <- list()
for (type in celltypes) {
  scesub <- sce[, sce$celltype %in% type]
  bulk <- convertBulkFilter(scesub, group = scesub$batch, filter = FALSE)
  cluster_list[[type]] <- add_bulk_metadata(bulk)
  message("Pseudobulk: ", type)
}
saveRDS(cluster_list, file = "object/bcell_bulk_list.rds")


# -----------------------------------------------------------------------------
# 3. RUV-corrected pseudobulk DE (edgeR LRT)
# -----------------------------------------------------------------------------
contrasts_pairs <- c(
  "Convertor_After-Convertor_Before",
  "Non_convertor_After-Non_convertor_Before",
  "Negative_control_After-Negative_control_Before",
  "Non_convertor_Before-Negative_control_Before",
  "Convertor_Before-Negative_control_Before",
  "Convertor_Before-Non_convertor_Before",
  "Non_convertor_After-Negative_control_After",
  "Convertor_After-Negative_control_After",
  "Convertor_After-Non_convertor_After"
)

design_levels <- c(
  "Negative_control_Before", "Negative_control_After",
  "Non_convertor_Before", "Non_convertor_After",
  "Convertor_Before", "Convertor_After"
)

marker_all <- list()
marker_top <- list()

for (type in celltypes) {
  bulk <- cluster_list[[type]]
  counts <- bulk$counts
  metadata <- bulk$metadata

  design <- model.matrix(~ 0 + condition_time, data = metadata)
  colnames(design) <- design_levels

  type_all <- list()
  type_top <- list()

  for (pair in contrasts_pairs) {
    topgenes <- fullDEtest(
      counts, metadata, design,
      group = "condition_time",
      contrasts = pair,
      method = "lrt",
      filter = TRUE,
      normalize = "TMM",
      RUV = TRUE,
      k = 1,
      nGenes = 5000,
      model_out = FALSE
    )
    type_all[[pair]] <- topgenes$top
    type_top[[pair]] <- topgenes$topSign

    out <- paste0(
      "tables/DE/RUV_", safe_name(type), "_", safe_name(pair), ".csv"
    )
    write.csv(topgenes$top, file = out, row.names = FALSE)
  }

  marker_all[[type]] <- type_all
  marker_top[[type]] <- type_top
  message("RUV DE: ", type)
}


# -----------------------------------------------------------------------------
# 4. Pathway enrichment (GSEA) on RUV DE results
# -----------------------------------------------------------------------------
pathway_pairs <- c(
  "Convertor_Before-Non_convertor_Before",
  "Convertor_After-Non_convertor_After"
)

pathway_all <- list()
pathway_top <- list()

for (type in celltypes) {
  gene_all <- marker_all[[type]]
  temp_all <- list()
  temp_top <- list()

  for (pair in pathway_pairs) {
    if (!pair %in% names(gene_all)) next
    top <- gene_all[[pair]]
    results <- profilerTest(
      top, method = "gse", lfc = 0.5, minGSSize = 5, cutoff = 1
    )
    sign <- results[results$p.adjust < 0.05, ]

    out <- paste0(
      "tables/DE/pathway_", safe_name(type), "_", safe_name(pair), ".csv"
    )
    write.csv(results, file = out, row.names = TRUE)

    temp_all[[pair]] <- results
    temp_top[[pair]] <- sign
  }

  pathway_all[[type]] <- temp_all
  pathway_top[[type]] <- temp_top
  message("Pathway: ", type)
}


# -----------------------------------------------------------------------------
# 5. Volcano plots (RUV DE results)
# -----------------------------------------------------------------------------
files <- list.files(
  path = "tables/DE",
  pattern = "^RUV_.*\\.csv$",
  full.names = TRUE
)

for (f in files) {
  top <- read.csv(f)
  if (!"genes" %in% colnames(top) && !"logFC" %in% colnames(top)) next
  if (!"genes" %in% colnames(top) && !is.null(rownames(top))) {
    top$genes <- rownames(top)
  }
  # Drop uninformative contig-like genes if present
  if ("genes" %in% colnames(top)) {
    drop <- grep("AC00", top$genes)
    if (length(drop) > 0) top <- top[-drop, ]
  }

  tt <- gsub("^RUV_|\\.csv$", "", basename(f))
  gra <- plot_volcano(top, title = tt)

  out <- file.path("plot/DE/volcano", paste0("volcano_", basename(f)))
  out <- sub("\\.csv$", ".pdf", out)
  pdf(out, width = 10, height = 7)
  print(gra)
  dev.off()
}


# -----------------------------------------------------------------------------
# 6. Save combined results
# -----------------------------------------------------------------------------
save(
  marker_all, marker_top, pathway_all, pathway_top, cluster_list,
  file = "object/bcells_DE_results.RData"
)

message("Done. DE tables: tables/DE/")
message("Done. Volcano plots: plot/DE/volcano/")
message("Done. Combined results: object/bcells_DE_results.RData")
