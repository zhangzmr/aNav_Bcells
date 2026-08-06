# aNav B-cell analysis

Analysis scripts for single-cell RNA-seq and BCR repertoire profiling of activated naive B cells in rheumatoid arthritis.

## Manuscript

**Senescent activated naive B cells promote anti-citrullinated antigen T cell responses and the transition to clinical rheumatoid arthritis**

Xiaohao Wu, Mengrui Zhang, Jae-Seung Moon, Eun Kyung Song, Joshua C Abrams, Laura S van Dam, Orr Sharpe, Marie Feser, Laurie Moss, Peggy P Ho, Melanie H Smith, Laura T Donlin, Jane H Buckner, Eddie A James, Gary S Firestein, Yuko Okamoto, Tobiaz V Lanz, Eric Meffre, V Michael Holers, Kevin D Deane, William H Robinson

## Custom packages

These scripts rely on packages developed for this study:

| Package | Description | Install |
|---------|-------------|---------|
| [bioscrna](https://github.com/zhangzmr/bioscrna) | Single-cell RNA-seq toolkit (loading, QC, reduction, clustering, plotting) | `remotes::install_github("zhangzmr/bioscrna")` |
| [scgtest](https://github.com/zhangzmr/scgtest) | Pseudobulk DE, RUV correction, pathway enrichment, volcano plots | `remotes::install_github("zhangzmr/scgtest")` |
| [screceptor](https://github.com/zhangzmr/screceptor) | BCR/TCR processing (cloning, germline, SHM, GEX merge, trees) | `remotes::install_github("zhangzmr/screceptor")` |

## Scripts

### Main blood cohort (RA B cells)

| Script | Description |
|--------|-------------|
| `scrna_processing.R` | Load 10X data, QC, normalize, integrate, cluster, annotate B-cell types, trajectory |
| `DE_analysis.R` | Pseudobulk DE (RUV + edgeR), pathway enrichment, volcano plots |
| `proportion_PCA.R` | Cell-type proportions and pseudobulk PCA |
| `cell_scores.R` | Pathway / UCell scores and visualization |
| `BCR_analysis.R` | BCR cloning, germline, SHM, merge with GEX, clonality, aNaive fate, trees |
| `aNaive_analysis.R` | Focused aNaive subset analysis (composition, markers, PCA, DE, BCR fate) |

### OA vs RA synovium

| Script | Description |
|--------|-------------|
| `OA_RA_scrna_processing.R` | All-cell processing and B-cell subset |
| `OA_RA_Bcells.R` | B-cell re-embedding, clustering, marker visualization |
| `OA_RA_BCR.R` | Synovial BCR load → clean → clone → germline → SHM |

### Naive-only cohort

| Script | Description |
|--------|-------------|
| `naive_only_scrna.R` | All-cell + naive subset processing and cell-type annotation |
| `naive_only_BCR.R` | Naive-cohort BCR processing, GEX merge, clonality plots |

## Usage

Set the working directory to the project root before running (paths in scripts are relative):

```r
setwd("/path/to/project")
source("scrna_processing.R")
```

Or from the command line:

```bash
Rscript scrna_processing.R
```

Suggested run order for the main blood cohort:

1. `scrna_processing.R`
2. `DE_analysis.R` / `proportion_PCA.R` / `cell_scores.R` (independent after processing)
3. `BCR_analysis.R`
4. `aNaive_analysis.R`

For BCR scripts, set the IMGT VDJ reference path (`vdj_reference` or `references`) before running.

## Expected layout

```
project/
├── data/                 # input 10X / BCR files (see script headers)
├── object/               # intermediate SingleCellExperiment objects
├── plot/                 # figures
├── tables/               # markers, DE tables
├── bcr_result/           # BCR intermediates and exports
└── *.R                   # analysis scripts
```

Exact input paths are documented at the top of each script.
