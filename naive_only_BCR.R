#!/usr/bin/env Rscript
# =============================================================================
# naive_only_BCR.R
#
# Naive-only cohort BCR analysis:
#   1. Load, clean, clone, germline, and SHM
#   2. Merge with annotated naive scRNA-seq
#   3. Clonality barplots by cell type
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
  library(dowser)
  library(dplyr)
  library(ggplot2)
})

# Project root: set this to the repository root before running
# setwd("/path/to/aNav_code")

similarity <- 0.75
ncores <- 6
# Path to IMGT VDJ reference directory (required for germlineBCR)
# vdj_reference <- "data/reference/human_vdj/vdj/"
vdj_reference <- NULL

dir.create("bcr_result/naive_only", showWarnings = FALSE, recursive = TRUE)
dir.create("plot/naive_only", showWarnings = FALSE, recursive = TRUE)
dir.create("object/naive_only", showWarnings = FALSE, recursive = TRUE)
out <- "bcr_result/naive_only"


# -----------------------------------------------------------------------------
# 1. Load and process BCR
# -----------------------------------------------------------------------------
# Expected layout:
#   data/naive_only/bcr/vdj/<sample>/BCR_data_sequences_filtered_igblast_db-pass.tsv
sample <- list.files("data/naive_only/bcr/vdj")
directory <- file.path(
  "data/naive_only/bcr/vdj", sample,
  "BCR_data_sequences_filtered_igblast_db-pass.tsv"
)

bcrlist <- loadBCR(directory = directory, sample = sample)
saveRDS(bcrlist, file = file.path(out, "BCR_raw.rds"))

bcrlist <- lapply(bcrlist, cleanBCR)
bcrlist_heavy <- lapply(lapply(bcrlist, subsetHC), cloneHC, similarity = similarity)
bcrlist_light <- lapply(
  lapply(bcrlist, subsetLC), cloneLC,
  bcrlist_heavy = bcrlist_heavy, col = "sample_id"
)

bcrTogether <- combineHL(
  bcrlist_heavy, bcrlist_light,
  col = "sample_id", clean_seq = FALSE, output = "list"
)
bcrTogether <- lapply(bcrTogether, subCloneLC, nproc = ncores, minseq = 1)

if (is.null(vdj_reference)) {
  stop("Set `vdj_reference` to the IMGT VDJ reference directory.")
}
reference <- readIMGT(dir = vdj_reference)
bcrTogether <- lapply(
  bcrTogether, germlineBCR,
  reference = reference, clone = "clone_subgroup_id",
  nproc = ncores, trim_lengths = TRUE
)

bcrlist_heavy <- lapply(lapply(bcrTogether, subsetHC), hypermuBCR, nproc = ncores)
bcrlist_light <- lapply(lapply(bcrTogether, subsetLC), hypermuBCR, nproc = ncores)

bcrCombine <- combineHL(
  bcrlist_heavy, bcrlist_light,
  clean_seq = FALSE, col = "sample_id", output = "combine"
)
saveRDS(bcrCombine, file = file.path(out, "BCR_combined.rds"))


# -----------------------------------------------------------------------------
# 2. Merge with annotated naive GEX
# -----------------------------------------------------------------------------
sce <- readRDS("object/naive_only/naive_annotated.rds")
sce_bcr <- mergeBCRtoGEX(sce, bcr = bcrCombine, intersect = TRUE)
saveRDS(sce_bcr, file = "object/naive_only/naive_bcr_merged.rds")

BCRgex <- mergeGEXtoBCR(
  sce, bcr = bcrCombine, intersect = FALSE,
  var = c("cellcluster", "celltype")
)
BCRgex <- BCRgex[BCRgex$sample_id %in% unique(sce$batch), ]
BCRgex$disease <- ifelse(grepl("RA", BCRgex$sample_id), "RA", "HC")
BCRgex$disease <- factor(BCRgex$disease, levels = c("HC", "RA"))
BCRgex <- addClonality(BCRgex, size_col = "cloneSizeHL", breaks = c(1, 5, 20))
saveRDS(BCRgex, file = file.path(out, "bcr_combined_annotated.rds"))
write.csv(BCRgex, file = file.path(out, "combined_bcr_annotated.csv"), row.names = FALSE)


# -----------------------------------------------------------------------------
# 3. Clonality barplots by cell type
# -----------------------------------------------------------------------------
heavy <- BCRgex[BCRgex$locus %in% "IGH" & !is.na(BCRgex$celltype), ]
col <- c(
  Large = "#982c2c", Medium = "#f8984f",
  Small = "#f3d78a", Singleton = "#f2f1e6"
)

pdf("plot/naive_only/BCR_clonal_by_celltype.pdf", width = 4, height = 5)
for (ct in levels(droplevels(factor(heavy$celltype)))) {
  heavy_sub <- heavy[heavy$celltype %in% ct, ]
  if (nrow(heavy_sub) == 0) next

  plot_df <- heavy_sub %>%
    count(sample_id, clonality) %>%
    group_by(sample_id) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  print(
    ggplot(plot_df, aes(x = sample_id, y = pct, fill = clonality)) +
      geom_col(position = "stack", width = 0.85) +
      scale_fill_manual(values = col) +
      theme_classic() +
      labs(x = NULL, y = "Percentage (%)", title = ct)
  )
}
dev.off()

message("Done. BCR combined: ", file.path(out, "BCR_combined.rds"))
message("Done. Annotated BCR: ", file.path(out, "bcr_combined_annotated.rds"))
message("Done. Plots: plot/naive_only/")
