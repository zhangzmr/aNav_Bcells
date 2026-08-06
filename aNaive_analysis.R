#!/usr/bin/env Rscript
# =============================================================================
# aNaive_analysis.R
#
# Activated naive (aNaive) B-cell focused analysis:
#   1. Subset naive populations and re-embed
#   2. Composition and marker visualization
#   3. Pseudobulk PCA of aNaive_CXCR5+
#   4. Differential expression (aNaive vs rNaive; conditions)
#   5. BCR clonal fate and alluvial plots
#
# Custom packages developed for this study:
#   bioscrna   Single-cell RNA-seq analysis toolkit built on Bioconductor
#              (data loading, QC, dimension reduction, clustering, and plotting
#              for SingleCellExperiment objects).
#              https://github.com/zhangzmr/bioscrna
#              remotes::install_github("zhangzmr/bioscrna")
#
#   scgtest    Differential expression toolkit for single-cell and pseudobulk
#              RNA-seq (edgeR, RUV correction, pathway enrichment, and volcano
#              plots).
#              https://github.com/zhangzmr/scgtest
#              remotes::install_github("zhangzmr/scgtest")
#
#   screceptor Process single-cell BCR/TCR repertoire data (loading, cleaning,
#              clonal clustering, germline reconstruction, SHM, merging with
#              gene expression, and phylogenetic tree helpers).
#              https://github.com/zhangzmr/screceptor
#              remotes::install_github("zhangzmr/screceptor")
# =============================================================================

suppressPackageStartupMessages({
  library(bioscrna)
  library(scgtest)
  library(BiocParallel)
  library(batchelor)
  library(dittoSeq)
  library(ggplot2)
  library(ggforce)
  library(ggalluvial)
  library(dplyr)
  library(tidyr)
  library(edgeR)
  library(limma)
  library(MatrixGenerics)
  library(viridis)
  library(ggsci)
})

# Project root: set this to the repository root before running
# setwd("/path/to/aNav_code")

dir.create("object/aNaive", showWarnings = FALSE, recursive = TRUE)
dir.create("plot/aNaive", showWarnings = FALSE, recursive = TRUE)
dir.create("tables/aNaive", showWarnings = FALSE, recursive = TRUE)

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

naive_levels <- c("rNaive", "aNaive_CXCR5+", "aNaive_CXCR5+_CD38hi")
naive_colors <- c("#A469BD", "#D5695D", "#F5B041")
pca_colors <- c(
  "dodgerblue", "blue", "green", "darkgreen", "deeppink", "darkorchid"
)

add_condition_metadata <- function(x, sample_col = "batch") {
  sid <- gsub("^B_", "", as.character(x[[sample_col]]))
  x$condition[sid %in% c(
    "PC1A", "PC1B", "PC2A", "PC2B", "PC3A", "PC3B", "PC4A", "PC4B", "PC5A", "PC5B"
  )] <- "Convertor"
  x$condition[sid %in% c(
    "PN1A", "PN1B", "PN2A", "PN2B", "PN3A", "PN3B", "PN4A", "PN4B", "PN5A", "PN5B"
  )] <- "Non_convertor"
  x$condition[sid %in% c(
    "NN1A", "NN1B", "NN2A", "NN2B", "NN3A", "NN3B"
  )] <- "Negative_control"

  x$timepoint[sid %in% c(
    "NN1A", "NN2A", "NN3A", "PC1A", "PC2A", "PC3A", "PC4A", "PC5A",
    "PN1A", "PN2A", "PN3A", "PN4A", "PN5A"
  )] <- "Before"
  x$timepoint[sid %in% c(
    "NN1B", "NN2B", "NN3B", "PC1B", "PC2B", "PC3B", "PC4B", "PC5B",
    "PN1B", "PN2B", "PN3B", "PN4B", "PN5B"
  )] <- "After"

  x$condition_time <- paste0(x$condition, "_", x$timepoint)
  x$condition <- factor(
    x$condition,
    levels = c("Negative_control", "Non_convertor", "Convertor")
  )
  x$timepoint <- factor(x$timepoint, levels = c("Before", "After"))
  x$condition_time <- factor(
    x$condition_time,
    levels = c(
      "Negative_control_Before", "Negative_control_After",
      "Non_convertor_Before", "Non_convertor_After",
      "Convertor_Before", "Convertor_After"
    )
  )
  x
}


# -----------------------------------------------------------------------------
# 1. Subset naive populations and re-embed
# -----------------------------------------------------------------------------
sce <- readRDS("object/bcells_annotated.rds")
sce <- add_condition_metadata(sce, sample_col = "batch")
sce <- sce[, sce$celltype %in% naive_levels]
sce$celltype <- factor(sce$celltype, levels = naive_levels)

# Exclude Ig variable genes from dimension reduction features
ig_filter <- c("IGLV", "IGHV", "IGKV", "IGLC")
ig_index <- unlist(lapply(ig_filter, function(x) grep(x, rownames(sce))))
ig_genes <- rownames(sce)[ig_index]
temp <- sce[!rownames(sce) %in% ig_genes, ]

temp <- dimensionReduction(temp, n = 2000, rank = 50)
temp_corrected <- fastMNN(temp, batch = temp$batch, d = 50, BPPARAM = MulticoreParam())
reducedDim(temp, "batch_corrected") <- reducedDim(temp_corrected, "corrected")
temp <- runUMAP(
  temp, dimred = "batch_corrected", name = "UMAP_corrected",
  BPPARAM = MulticoreParam()
)

for (rd in c("PCA", "UMAP", "batch_corrected", "UMAP_corrected")) {
  if (rd %in% reducedDimNames(temp)) {
    reducedDim(sce, rd) <- reducedDim(temp, rd)
  }
}
saveRDS(sce, file = "object/aNaive/aNaive_updated.rds")


# -----------------------------------------------------------------------------
# 2. UMAP and composition
# -----------------------------------------------------------------------------
pdf("plot/aNaive/UMAP_celltype.pdf", height = 6, width = 7)
print(
  reducedPlot(
    sce, var = "celltype", reduction.use = "UMAP_corrected",
    size = 0.5, do.raster = TRUE, do.label = FALSE
  ) +
    scale_color_manual(values = naive_colors) +
    theme_classic()
)
dev.off()

pdf("plot/aNaive/UMAP_by_condition_time.pdf", height = 5, width = 7)
for (ct in levels(sce$condition_time)) {
  temp <- sce[, sce$condition_time %in% ct]
  if (ncol(temp) == 0) next
  print(
    reducedPlot(
      temp, var = "celltype", reduction.use = "UMAP_corrected",
      size = 1, do.raster = TRUE, do.label = FALSE
    ) +
      scale_color_manual(values = naive_colors) +
      theme_classic() +
      ggtitle(ct)
  )
}
dev.off()

pdf("plot/aNaive/barplot_composition.pdf", height = 5, width = 6)
print(
  dittoBarPlot(
    sce, "celltype", group.by = "condition_time",
    color.panel = naive_colors, retain.factor.levels = TRUE
  )
)
print(
  dittoBarPlot(
    sce, "celltype", group.by = "condition",
    color.panel = naive_colors, retain.factor.levels = TRUE
  )
)
dev.off()

prop <- dittoBarPlot(
  sce, "celltype", group.by = "condition_time",
  retain.factor.levels = TRUE, data.out = TRUE
)$data
write.csv(prop, file = "tables/aNaive/celltype_proportion.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 3. Marker visualization
# -----------------------------------------------------------------------------
# Before-treatment rNaive vs aNaive grouping
before <- sce[, sce$timepoint %in% "Before"]
before$group <- ifelse(
  before$celltype %in% "rNaive", "rNaive", "aNaive"
)
before$group <- paste0(before$group, "_", before$condition)
before$group <- factor(
  before$group,
  levels = c(
    "rNaive_Negative_control", "rNaive_Non_convertor", "rNaive_Convertor",
    "aNaive_Negative_control", "aNaive_Non_convertor", "aNaive_Convertor"
  )
)

act_genes <- c(
  "CD38", "CD69", "JUN", "JUND", "CD83", "FOS", "FOSB", "JUNB",
  "NFKB1", "NFKB2", "LYN", "SYK", "CXCR5", "IGHM", "IGHD", "TCL1A"
)
intergene <- intersect(act_genes, rownames(before))

pdf("plot/aNaive/heatmap_activation_markers.pdf", height = 5, width = 4)
averageHeatmap(
  before, genes = intergene, groups = before$celltype,
  show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE
)
dev.off()

pdf("plot/aNaive/dotplot_activation_markers.pdf", height = 6, width = 8)
print(
  dittoDotPlot(before, vars = intergene, group.by = "group") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()


# -----------------------------------------------------------------------------
# 4. Pseudobulk PCA of aNaive_CXCR5+
# -----------------------------------------------------------------------------
scesub <- sce[, sce$celltype %in% "aNaive_CXCR5+"]
bulk <- convertBulkFilter(scesub, group = scesub$batch, filter = FALSE)

metadata <- bulk$metadata
colnames(metadata) <- c("sample", "nCells")
metadata <- add_condition_metadata(metadata, sample_col = "sample")
bulk$metadata <- metadata

counts <- bulk$counts
rowcounts <- rowSums2(counts)
keep_genes <- names(rowcounts)[rowcounts >= 5]
counts <- counts[rownames(counts) %in% keep_genes, , drop = FALSE]
logcounts <- log2(counts + 1)

design <- model.matrix(~ 0 + condition_time, data = metadata)
colnames(design) <- levels(metadata$condition_time)

removed <- removeBatchEffect(logcounts, design = design)
pca_res <- prcomp(t(removed), scale. = TRUE)
plot_df <- cbind(as.data.frame(pca_res$x[, 1:2]), metadata)

pdf("plot/aNaive/aNaive_CXCR5_PCA.pdf", height = 5, width = 6)
print(
  ggplot(plot_df, aes(x = PC1, y = PC2, color = condition_time)) +
    geom_point(size = 3) +
    geom_mark_ellipse(aes(color = condition_time), linetype = 2, size = 0.7) +
    scale_color_manual(values = pca_colors) +
    theme_classic() +
    ggtitle("aNaive_CXCR5+ pseudobulk PCA")
)
dev.off()
saveRDS(list(pca = pca_res, metadata = metadata, removed = removed),
        file = "object/aNaive/aNaive_CXCR5_PCA.rds")


# -----------------------------------------------------------------------------
# 5. Differential expression
# -----------------------------------------------------------------------------
# aNaive vs rNaive within each condition_time
ct_pairs <- c("aNaive_CXCR5+--rNaive", "aNaive_CXCR5+_CD38hi--rNaive")
for (ct in levels(sce$condition_time)) {
  temp <- sce[, sce$condition_time %in% ct]
  if (length(unique(temp$celltype)) < 2) next
  markers <- wilcoxTest(
    sce = temp, group_vector = temp$celltype, pairs = ct_pairs
  )
  for (p in names(markers)) {
    out <- paste0(
      "tables/aNaive/DE_", safe_name(ct), "_", safe_name(p), ".csv"
    )
    write.csv(markers[[p]], file = out, row.names = FALSE)
  }
  message("DE celltype: ", ct)
}

# Condition contrasts within aNaive_CXCR5+
scesub <- sce[, sce$celltype %in% "aNaive_CXCR5+"]
cond_pairs <- c(
  "Convertor_Before--Non_convertor_Before",
  "Convertor_After--Non_convertor_After",
  "Convertor_Before--Negative_control_Before",
  "Convertor_After--Negative_control_After",
  "Non_convertor_Before--Negative_control_Before",
  "Non_convertor_After--Negative_control_After"
)
markers <- wilcoxTest(
  sce = scesub, group_vector = scesub$condition_time, pairs = cond_pairs
)
for (p in names(markers)) {
  out <- paste0("tables/aNaive/DE_aNaive_CXCR5_", safe_name(p), ".csv")
  write.csv(markers[[p]], file = out, row.names = FALSE)
}


# -----------------------------------------------------------------------------
# 6. BCR clonal fate and alluvial plots
# -----------------------------------------------------------------------------
# Prefer annotated BCR table from BCR_analysis.R
bcr_path <- "bcr_result/bcr_combined_annotated.rds"
if (!file.exists(bcr_path)) {
  warning("BCR table not found at ", bcr_path, "; skipping clonal fate / alluvial.")
} else {
  bcr_table <- readRDS(bcr_path)
  bcr_table <- bcr_table[!is.na(bcr_table$celltype), ]
  if (!"batch_clone" %in% colnames(bcr_table) &&
      all(c("sample_id", "clone_id") %in% colnames(bcr_table))) {
    bcr_table$batch_clone <- paste0(bcr_table$sample_id, "_", bcr_table$clone_id)
  }

  # IGHM+ aNaive_CXCR5+ clones at Before
  expanded <- bcr_table
  if ("cloneSizeHC" %in% colnames(bcr_table)) {
    expanded <- bcr_table[bcr_table$cloneSizeHC > 1, ]
  }

  aNav <- expanded[
    expanded$celltype %in% "aNaive_CXCR5+" &
      expanded$c_call %in% "IGHM" &
      expanded$timepoint %in% "Before",
  ]
  clone_ids <- unique(aNav$clone_id)
  cell_id_before <- unique(aNav$batch_barcode)

  before <- expanded[expanded$batch_barcode %in% cell_id_before, ]
  after <- expanded[
    expanded$clone_id %in% clone_ids & expanded$timepoint %in% "After",
  ]
  mucl <- rbind(before, after)
  write.csv(mucl, file = "tables/aNaive/aNaive_related_clones.csv", row.names = FALSE)

  # Related clone percentage per sample
  if ("batch_clone" %in% colnames(expanded)) {
    related_clones <- unique(aNav$batch_clone)
    expanded$aNav_related <- ifelse(
      expanded$batch_clone %in% related_clones, "aNav_related", "Others"
    )
    expanded$aNav_related[expanded$celltype %in% "aNaive_CXCR5+"] <- "aNaive"

    percent <- lapply(unique(expanded$sample_id), function(i) {
      temp <- expanded[expanded$sample_id %in% i, ]
      total <- length(unique(temp$batch_clone))
      related <- length(unique(
        temp$batch_clone[temp$aNav_related %in% c("aNav_related", "aNaive")]
      ))
      data.frame(sample_id = i, percent = related / total)
    })
    percent <- bind_rows(percent)
    write.csv(
      percent,
      file = "tables/aNaive/aNav_related_clone_percent.csv",
      row.names = FALSE
    )
  }

  # Alluvial: Before vs After cell-type composition of related clones
  type <- c("Non_convertor", "Convertor", "Negative_control")
  pdf("plot/aNaive/aNaive_alluvial_plot.pdf", width = 8, height = 6)
  for (t in type) {
    dat <- mucl[mucl$condition %in% t, ]
    dat <- dat[dat$locus %in% "IGH", ]
    if (nrow(dat) == 0) next

    pair <- c("Before", "After")
    inter <- unique(as.character(dat$celltype))
    freq_time <- data.frame(type = inter, Before = 0, After = 0)
    for (c in inter) {
      temp <- dat[dat$celltype %in% c, ]
      freq_time[freq_time$type %in% c, "Before"] <-
        sum(temp$timepoint %in% "Before")
      freq_time[freq_time$type %in% c, "After"] <-
        sum(temp$timepoint %in% "After")
    }

    temp1 <- data.frame(
      type = freq_time$type, size = freq_time$Before, Sample = "Before"
    )
    temp2 <- data.frame(
      type = freq_time$type, size = freq_time$After, Sample = "After"
    )
    Con.df <- rbind(temp1, temp2)
    Con.df$Sample <- factor(Con.df$Sample, levels = pair)
    Con.df <- Con.df[Con.df$size > 0, ]
    if (nrow(Con.df) == 0) next

    print(
      ggplot(
        Con.df,
        aes(
          x = Sample, fill = type, stratum = type, alluvium = type,
          y = size, label = type
        )
      ) +
        geom_stratum() +
        geom_flow(stat = "alluvium") +
        theme_classic() +
        ggtitle(t) +
        theme(axis.title.x = element_blank())
    )
  }
  dev.off()

  # Fate barplots (non-naive members of related clones)
  related_nn <- mucl[
    !mucl$celltype %in% c("aNaive_CXCR5+", "aNaive_CXCR5+_CD38hi", "rNaive"),
  ]
  if (nrow(related_nn) > 0) {
    related_nn$temp <- paste0(related_nn$condition_time, "--", related_nn$celltype)
    freq <- as.data.frame(table(related_nn$temp))
    freq$condition <- gsub("--.*", "", freq$Var1)
    freq$celltype <- gsub(".*--", "", freq$Var1)
    freq$Freq <- freq$Freq / ave(freq$Freq, freq$condition, FUN = sum)

    pdf("plot/aNaive/aNaive_related_fate_barplot.pdf", width = 11, height = 6)
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

message("Done. aNaive object: object/aNaive/aNaive_updated.rds")
message("Done. Plots: plot/aNaive/")
message("Done. Tables: tables/aNaive/")
