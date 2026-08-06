#!/usr/bin/env Rscript
# =============================================================================
# OA_RA_BCR.R
#
# OA vs RA synovium BCR repertoire processing:
#   1. Load BCR annotations
#   2. Clean sequences
#   3. Clone heavy and light chains
#   4. Germline reconstruction
#   5. Somatic hypermutation (SHM)
#   6. Export combined BCR tables
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
})

# Project root: set this to the repository root before running
# setwd("/path/to/aNav_code")

# Parameters
similarity <- 0.75
ncores <- 6
# Path to IMGT VDJ reference directory (required for germlineBCR)
# vdj_reference <- "data/reference/human_vdj/vdj/"
vdj_reference <- NULL

dir.create("bcr_result/OA_RA", showWarnings = FALSE, recursive = TRUE)
out <- "bcr_result/OA_RA"


# -----------------------------------------------------------------------------
# 1. Load BCR annotations
# -----------------------------------------------------------------------------
# Expected layout:
#   data/OA_RA/bcr/<sample>/BCR_data_sequences_igblast_db-pass.tsv
sample <- list.files("data/OA_RA/bcr")
directory <- file.path(
  "data/OA_RA/bcr", sample, "BCR_data_sequences_igblast_db-pass.tsv"
)
sample_sheet <- data.frame(sample = sample, directory = directory)
write.csv(sample_sheet, file = file.path(out, "bcr_sheet.csv"), row.names = FALSE)

bcrlist <- loadBCR(directory = sample_sheet$directory, sample = sample_sheet$sample)
saveRDS(bcrlist, file = file.path(out, "BCR_raw.rds"))


# -----------------------------------------------------------------------------
# 2. Clean BCR
# -----------------------------------------------------------------------------
bcrlist <- lapply(bcrlist, cleanBCR)
saveRDS(bcrlist, file = file.path(out, "BCR_cleaned.rds"))


# -----------------------------------------------------------------------------
# 3. Clone heavy chains; transfer clone IDs to light chains
# -----------------------------------------------------------------------------
bcrlist_heavy <- lapply(bcrlist, subsetHC)
bcrlist_heavy <- lapply(bcrlist_heavy, cloneHC, similarity = similarity)

bcrlist_light <- lapply(bcrlist, subsetLC)
bcrlist_light <- lapply(
  bcrlist_light, cloneLC,
  bcrlist_heavy = bcrlist_heavy, col = "sample_id"
)

bcrTogether <- combineHL(
  bcrlist_heavy = bcrlist_heavy,
  bcrlist_light = bcrlist_light,
  col = "sample_id",
  clean_seq = FALSE,
  output = "list"
)
saveRDS(bcrTogether, file = file.path(out, "BCR_heavy_cloned.rds"))


# -----------------------------------------------------------------------------
# 4. Resolve light-chain subclones
# -----------------------------------------------------------------------------
bcrTogether <- lapply(bcrTogether, subCloneLC, nproc = ncores, minseq = 1)
saveRDS(bcrTogether, file = file.path(out, "BCR_heavy_light_cloned.rds"))


# -----------------------------------------------------------------------------
# 5. Germline reconstruction
# -----------------------------------------------------------------------------
if (is.null(vdj_reference)) {
  stop(
    "Set `vdj_reference` to the IMGT VDJ reference directory ",
    "before running germlineBCR / hypermuBCR."
  )
}
reference <- readIMGT(dir = vdj_reference)
bcrTogether <- lapply(
  bcrTogether, germlineBCR,
  reference = reference,
  clone = "clone_subgroup_id",
  nproc = ncores,
  trim_lengths = TRUE
)
saveRDS(bcrTogether, file = file.path(out, "BCR_heavy_light_cloned_germline.rds"))


# -----------------------------------------------------------------------------
# 6. SHM on heavy and light chains
# -----------------------------------------------------------------------------
bcrlist_heavy <- lapply(bcrTogether, subsetHC)
bcrlist_light <- lapply(bcrTogether, subsetLC)

bcrlist_heavy <- lapply(
  bcrlist_heavy, hypermuBCR,
  nproc = ncores, cloneColumn = "clone_subgroup_id"
)
bcrlist_light <- lapply(
  bcrlist_light, hypermuBCR,
  nproc = ncores, cloneColumn = "clone_subgroup_id"
)

bcrTogether <- combineHL(
  bcrlist_heavy = bcrlist_heavy,
  bcrlist_light = bcrlist_light,
  clean_seq = FALSE,
  col = "sample_id",
  output = "list"
)
bcrCombine <- combineHL(
  bcrlist_heavy = bcrlist_heavy,
  bcrlist_light = bcrlist_light,
  clean_seq = FALSE,
  col = "sample_id",
  output = "combine"
)

saveRDS(bcrTogether, file = file.path(out, "BCR_heavy_light_cloned_germline_SHM.rds"))
saveRDS(bcrTogether, file = file.path(out, "BCR_list.rds"))
saveRDS(bcrCombine, file = file.path(out, "BCR_combined.rds"))


# -----------------------------------------------------------------------------
# 7. Export per-sample and combined tables
# -----------------------------------------------------------------------------
for (n in names(bcrTogether)) {
  write.csv(
    bcrTogether[[n]],
    file = file.path(out, paste0(n, "_bcr_table.csv")),
    row.names = FALSE
  )
}
write.csv(
  bcrCombine,
  file = file.path(out, "combined_bcr_table.csv"),
  row.names = FALSE
)

message("Done. BCR list: ", file.path(out, "BCR_list.rds"))
message("Done. BCR combined: ", file.path(out, "BCR_combined.rds"))
message("Done. Tables: ", out, "/")
