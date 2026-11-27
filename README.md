# Aging remodels the neuronal activity-dependent transcriptome

Analysis code for Um et al. (2026): single-cell RNA-seq of the mouse hippocampus across aging, with
an enhanced novel environment (eNE) and chemogenetic activation (DREADD).

Preprint: https://www.biorxiv.org/content/10.64898/2026.11.17.767442v1

## Contents

| Folder | Contents | Figures |
|---|---|---|
| `00_preprocessing` | Per-sample QC, ambient RNA removal (SoupX) and doublet removal (scDblFinder) | Methods |
| `01_DESeq2` | eNE pseudobulk DESeq2 (all comparisons), QC, cell-type annotation, ARG programs, IN and astrocyte subtypes | Fig. 2; Extended Data Figs. 2-6 |
| `02_Dynamic_Range_Compression` | Dynamic range compression, failed induction and basal alteration | Fig. 3; Extended Data Fig. 7 |
| `03_DREADD` | DREADD pseudobulk DESeq2 (covariate and downsampling approaches) and DREADD figures | Fig. 4; Extended Data Fig. 8 |
| `04_Buffering` | Buffering analysis (held-out cosine, buffering index, gene tiers), age-emergent ARGs, cell-cell communication | Fig. 5; Extended Data Figs. 9, 10 |
| `05_TF` | Transcription factor activity (decoupleR, CollecTRI) and eNE vs DREADD comparison | Fig. 6; Extended Data Fig. 11 |
| `demo` | Small simulated dataset that runs the core statistical steps | - |

Most script names start with the figure panels they produce (e.g. `Fig4e_ED8g_Compression.R`).

## 1. System requirements

**Operating system.** The full analysis was run on Linux (Ubuntu 22.04.5 LTS, Washington University RIS compute cluster). The demo was also tested on Linux with R 4.3.3 and DESeq2 1.42.0; any OS that runs R (Linux, macOS, Windows) should work.

**Software.**

| Package | Version tested |
|---|---|
| R | 4.4.0 |
| Seurat | 5.6.0 |
| SeuratObject | 5.4.0 |
| BPCells | 0.3.1 |
| SoupX | 1.6.2 |
| scDblFinder | 1.20.2 |
| DESeq2 | 1.46.0 |
| decoupleR | 2.12.0 |
| OmnipathR | 3.14.0 |
| clusterProfiler | 4.14.6 |
| org.Mm.eg.db | 3.20.0 |
| CellChat | 1.6.1 |
| NeuronChat | 1.0.0 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| ggplot2 | 4.0.3 |
| ggrepel | 0.9.8 |
| patchwork | 1.3.2 |

**Hardware.** The demo runs on any standard desktop or laptop (< 1 GB RAM). The full analysis loads Seurat
objects of 818,098 (eNE) and 322,130 (DREADD) cells with counts stored on disk by BPCells; we ran it on a
64-core cluster node and recommend at least 64 GB of RAM. No GPU or other non-standard hardware is required.

## 2. Installation

Install R (>= 4.3) from https://cran.r-project.org, then in R:

```r
install.packages(c("Seurat", "dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork", "remotes", "BiocManager"))
BiocManager::install(c("DESeq2", "scDblFinder", "decoupleR", "OmnipathR", "clusterProfiler", "org.Mm.eg.db"))
install.packages("SoupX")
remotes::install_github("bnprks/BPCells/r")          # needs a C++ compiler and HDF5
remotes::install_github("sqjin/CellChat")              # CellChat v1
remotes::install_github("Wei-BioMath/NeuronChat")
```

Then clone this repository:

```bash
git clone https://github.com/KoreanUm97/Um2026-experience-buffers-aging-transcriptome.git
```

**Typical install time** on a standard desktop: about 10 minutes for the demo (DESeq2 only) and 30-60 minutes
for all packages, most of it compiling BPCells and Bioconductor packages.

## 3. Demo

The demo simulates pseudobulk counts for 2 cell types x 3 ages x 2 conditions x 3 mice, with activity-regulated
genes (ARGs) whose induction contracts with age and aging genes that stimulation partially reverses. It then runs
the same core steps as the full pipelines: pseudobulk DESeq2 (design `~ grp`, genes with >= 10 counts in >= 2
samples), ARG calls (padj < 0.05, |log2FC| >= 0.25), the dynamic range of the young ARG program at each age, and
buffering (cosine between the aging and stimulation fold-change vectors, and the buffering index
BI = sum(OEEYHC^2) / sum(OHCYHC^2)).

Only DESeq2 is needed. From the repository root:

```bash
Rscript demo/run_demo.R
```

**Expected output** (`demo/demo_output/`): one DESeq2 table per cell type and comparison, `demo_summary.csv`
and `demo_plots.pdf`. The summary printed to the console should read:

| celltype | ARGs young | ARGs middle | ARGs old | median log2FC of young ARGs (young / middle / old) | aging genes | cosine aging vs EE | buffering index |
|---|---|---|---|---|---|---|---|
| CA1_ExNs | 150 | 97 | 2 | 1.44 / 0.92 / 0.42 | 124 | -0.87 | 0.39 |
| Astrocytes | 151 | 102 | 2 | 1.43 / 0.96 / 0.45 | 127 | -0.89 | 0.36 |

The simulated ARG induction contracts with age (1.5 -> 1.0 -> 0.5 log2FC), so the number and fold change of
induced genes fall with age; stimulation reverses 40% of the aging change, so the cosine is strongly negative
and the buffering index is close to the simulated (1 - 0.4)^2 = 0.36. Exact values can differ slightly with
the DESeq2 version.

**Expected run time:** under 1 minute on a standard desktop (about 15 seconds in our tests).

## 4. Running the analysis on the data

**Data.** Raw and processed data are in GEO under accession **GSEXXXXXX**. The analysis starts from the
processed Seurat objects. Rebuild each from the GEO files (shown for eNE; the same for DREADD):

```r
library(Seurat); library(BPCells)
counts <- open_matrix_hdf5("eNE_counts_bpcells.h5", group = "counts")
x      <- readRDS("eNE_metadata_reductions.rds")
obj    <- CreateSeuratObject(counts, meta.data = x$meta.data)
for (r in names(x$reductions)) obj[[r]] <- x$reductions[[r]]
obj    <- NormalizeData(obj)
saveRDS(obj, "eNE_object.rds")
```

`*_cell_metadata.csv.gz` holds the same per-cell annotation (sample, sex, age, condition, cell type, UMAP) as
plain text, and `GEO_sample_key.csv` lists each mouse with its probe barcode.

**Running the scripts.** Set the paths at the top of each script (`PATH_TO_...` placeholders), then run in
this order:

1. `00_preprocessing` (only needed to redo QC from the Cell Ranger h5 files; the processed objects on GEO already
   contain the quality-controlled cells and counts used in the paper)
2. `01_DESeq2/Fig2_ED4_Pipeline.R` (writes the DESeq2 tables used downstream, then Fig. 2), then `ED2_QC_Panel.R`,
   `ED3_ClusterAnnotation.R`, `ED5_IN_Subclusters.R`, `ED6_Astro_Subclusters.R`; `03_DREADD/DREADD_DESeq2.R`
3. `02_Dynamic_Range_Compression`, rest of `03_DREADD`
4. `04_Buffering`: `Buffering_Pseudobulk.R`, then `Buffering_Pipeline.R`, `DREADD_Buffering_HeldOut.R`,
   `Buffering_DREADD_vs_NE.R`, then the remaining scripts
5. `05_TF`: `TF_Buffering_Significance.R`, then `Fig6_ED11_TF_Figures.R`

Males and females are analysed separately throughout.

**Reproducing the figures.** Each folder's table row above lists the figures it produces, and each figure
script is named after its panels. Running the steps in order reproduces the quantitative results in the
manuscript; random steps (null distributions, downsampling) use the seeds set in each script.

## License

MIT (see `LICENSE`).

## Citation

Um J, Lei HC, Nakamura K, Tycksen E, Adelman JW, Pollina EA. Aging remodels the neuronal activity-dependent
transcriptome. (2026).
