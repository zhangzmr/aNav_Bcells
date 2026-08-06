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

| Script | Description |
|--------|-------------|
| `scrna_processing.R` | scRNA-seq processing, clustering, and annotation |
| `DE_analysis.R` | Pseudobulk differential expression and pathway enrichment |
| `proportion_PCA.R` | Cell-type proportions and pseudobulk PCA |
| `cell_scores.R` | Pathway / UCell scores |
| `BCR_analysis.R` | BCR repertoire analysis and merge with gene expression |
| `aNaive_analysis.R` | Activated naive B-cell focused analysis |
| `OA_RA_scrna_processing.R` | OA/RA synovium scRNA-seq processing |
| `OA_RA_Bcells.R` | OA/RA synovium B-cell analysis |
| `OA_RA_BCR.R` | OA/RA synovium BCR analysis |
| `naive_only_scrna.R` | Naive-only cohort scRNA-seq processing |
| `naive_only_BCR.R` | Naive-only cohort BCR analysis |

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

For BCR scripts, set the IMGT VDJ reference path (`vdj_reference` or `references`) before running. Input paths and outputs are documented in each script header.
