# Preprocess_Samples.R
# Per-sample preprocessing of 10x Flex data, from the Cell Ranger filtered matrices
# to quality-filtered, ambient-RNA-corrected, doublet-free Seurat objects.
#
#   1. Load Cell Ranger filtered (cells) and raw (droplets) count matrices
#   2. QC filter: genes per cell between the lower limit (25th percentile, applied only
#      when the sample's minimum is below 200 genes) and Q3 + 1.5 x IQR;
#      >= 100 UMIs; <= 10% mitochondrial reads
#   3. Ambient RNA removal with SoupX (clusters: LogNormalize, 30 PCs, Louvain res 0.3)
#   4. Doublet removal with scDblFinder (singlets kept)
#   5. Cell IDs prefixed with the sample name; log normalization
#
# Input : sample sheet (CSV) with columns
#           sample, condition, sex, age, filtered_h5, raw_h5
# Output: <out_dir>/<sample>/<sample>.rds, <sample>_metadata.csv
#         <out_dir>/Preprocessing_Summary.csv (cells kept at each step)
# Existing per-sample outputs are skipped.

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(scDblFinder)
})

# ---- settings ----------------------------------------------------------------
if (!exists("sample_sheet")) sample_sheet <- "PATH_TO_SAMPLE_SHEET.csv"   # one row per sample (see header)
if (!exists("out_dir"))      out_dir      <- "PATH_TO_PREPROCESSING_OUTPUT"

RESOLUTION   <- 0.3
N_DIMS       <- 30
MITO_MAX     <- 10     # percent
UMI_MIN      <- 100
GENE_FLOOR   <- 200    # lower limit applied only if a sample's minimum genes/cell is below this
MIN_CELLS    <- 100    # samples with fewer cells are dropped
if (!exists("SEED")) SEED <- 42

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
options(future.globals.maxSize = 70e9)

# ---- helpers -----------------------------------------------------------------
read_counts <- function(h5) {
  m <- Read10X_h5(h5, use.names = TRUE, unique.features = TRUE)
  if (is.list(m)) m <- m[[1]]          # Gene Expression is the first modality
  m
}

cluster_cells <- function(obj, resolution = RESOLUTION, dims = N_DIMS) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
  obj
}

# ---- run ---------------------------------------------------------------------
sheet <- read.csv(sample_sheet, stringsAsFactors = FALSE)
stopifnot(all(c("sample", "condition", "sex", "age", "filtered_h5", "raw_h5") %in% names(sheet)))

summary_rows <- list()

for (i in seq_len(nrow(sheet))) {
  s        <- sheet$sample[i]
  s_dir    <- file.path(out_dir, s)
  rds_file <- file.path(s_dir, paste0(s, ".rds"))
  sum_file <- file.path(s_dir, paste0(s, "_step_counts.csv"))
  dir.create(s_dir, showWarnings = FALSE)

  if (file.exists(rds_file) && file.exists(sum_file)) {
    message(s, ": exists, skipping")
    summary_rows[[s]] <- read.csv(sum_file)
    next
  }
  message("\n", s, ": loading")

  counts <- read_counts(sheet$filtered_h5[i])
  raw    <- read_counts(sheet$raw_h5[i])

  obj <- CreateSeuratObject(counts = counts, project = s)
  obj$orig.ident <- s
  obj$condition  <- sheet$condition[i]
  obj$sex        <- sheet$sex[i]
  obj$age        <- sheet$age[i]
  n_loaded <- ncol(obj)

  if (n_loaded < MIN_CELLS) { message(s, ": < ", MIN_CELLS, " cells, dropped"); next }

  # 2. QC filter
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^mt-")
  nf <- obj$nFeature_RNA
  gene_lower <- if (min(nf) >= GENE_FLOOR) -Inf else base::unname(quantile(nf, 0.25))
  gene_upper <- base::unname(quantile(nf, 0.75)) + 1.5 * IQR(nf)
  obj <- subset(obj, subset = nFeature_RNA >= gene_lower & nFeature_RNA <= gene_upper &
                              percent.mt <= MITO_MAX & nCount_RNA >= UMI_MIN)
  n_qc <- ncol(obj)
  if (n_qc < MIN_CELLS) { message(s, ": < ", MIN_CELLS, " cells after QC, dropped"); next }

  # 3. Ambient RNA (SoupX)
  clus <- cluster_cells(obj)$seurat_clusters
  raw  <- raw[rownames(raw) %in% rownames(obj), , drop = FALSE]
  rho  <- NA
  sc   <- SoupChannel(tod = raw, toc = LayerData(obj, layer = "counts"), calcSoupProfile = TRUE)
  sc   <- setClusters(sc, setNames(as.character(clus), colnames(obj)))
  res  <- tryCatch({
    sc  <- autoEstCont(sc, doPlot = FALSE)
    rho <- sc$metaData$rho[1]
    adjustCounts(sc, roundToInt = TRUE)
  }, error = function(e) { message(s, ": SoupX failed (", conditionMessage(e), "), counts left uncorrected"); NULL })
  if (!is.null(res)) {
    obj[["original.counts"]] <- CreateAssayObject(counts = LayerData(obj, layer = "counts"))
    obj <- SetAssayData(obj, assay = "RNA", layer = "counts", new.data = res[, colnames(obj)])
  }
  rm(raw, sc, res); gc()

  # 4. Doublets (scDblFinder, no cluster information)
  set.seed(SEED)
  dbl <- scDblFinder(LayerData(obj, assay = "RNA", layer = "counts"))
  obj$scDblFinder.class <- dbl$scDblFinder.class
  obj$scDblFinder.score <- dbl$scDblFinder.score
  n_doublets <- sum(obj$scDblFinder.class == "doublet")
  obj <- subset(obj, subset = scDblFinder.class == "singlet")

  # 5. Unique cell IDs, normalization, per-sample clusters
  obj <- RenameCells(obj, add.cell.id = s)
  obj <- cluster_cells(obj)
  obj <- RunUMAP(obj, dims = 1:N_DIMS, verbose = FALSE)

  obj@misc$qc <- list(gene_lower = gene_lower, gene_upper = gene_upper,
                      mito_max = MITO_MAX, umi_min = UMI_MIN, soupx_rho = rho)

  step <- data.frame(sample = s, n_loaded = n_loaded, gene_lower = gene_lower, gene_upper = gene_upper,
                     n_after_qc = n_qc, soupx_rho = rho, n_doublets = n_doublets, n_final = ncol(obj))
  write.csv(obj[[]], file.path(s_dir, paste0(s, "_metadata.csv")))
  saveRDS(obj, rds_file)
  write.csv(step, sum_file, row.names = FALSE)
  summary_rows[[s]] <- step
  message(s, ": ", n_loaded, " loaded -> ", n_qc, " after QC -> ", ncol(obj), " final")
  rm(obj); gc()
}

summary_tab <- do.call(rbind, summary_rows)
write.csv(summary_tab, file.path(out_dir, "Preprocessing_Summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary_tab)
