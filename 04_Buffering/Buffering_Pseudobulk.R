################################################################################
### Buffering_Pseudobulk.R
###
### Builds the NE pseudobulk cache used by all buffering scripts:
###   <base_path>/Buffering_LOO_Split/pseudobulk_counts_and_meta.RDS
###   = list(pb_all = counts summed per cell type x mouse, pb_meta = sample metadata)
### Same aggregation as NE_DESeq2_Pipeline.R. Run once, before Buffering_Pipeline.R.
################################################################################

suppressPackageStartupMessages({ library(BPCells); library(Seurat) })

base_path <- "PATH_TO_NE_OUTPUT/"
NE_OBJECT <- "PATH_TO_NE_OBJECT.rds"        # integrated NE Seurat object (celltype, orig.ident, experimental.condition, age, sex)
out_dir   <- paste0(base_path, "Buffering_LOO_Split/")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
pb_cache  <- file.path(out_dir, "pseudobulk_counts_and_meta.RDS")

if (file.exists(pb_cache)) {
  cat("Pseudobulk cache already exists:", pb_cache, "\n")
} else {
  if (!exists("df_harmony")) df_harmony <- readRDS(NE_OBJECT)
  cat("=== Aggregating pseudobulk counts ===\n")
  pb_all <- suppressWarnings(
    AggregateExpression(
      object               = df_harmony,
      assays               = "RNA",
      return.seurat        = FALSE,
      group.by             = c("celltype","orig.ident",
                                "experimental.condition","age","sex"),
      normalization.method = NULL,
      slot                 = "counts",
      verbose              = TRUE
    )$RNA
  )
  split_names <- strsplit(colnames(pb_all),"_")
  parse_right <- function(parts,from_end) parts[length(parts)-from_end+1]
  pb_meta <- data.frame(
    sample     = colnames(pb_all),
    sex        = sapply(split_names,parse_right,from_end=1),
    age        = sapply(split_names,parse_right,from_end=2),
    condition  = sapply(split_names,parse_right,from_end=3),
    orig.ident = sapply(split_names,parse_right,from_end=4),
    celltype   = sapply(split_names,function(x)
                   paste(x[seq_len(length(x)-4)],collapse=" ")),
    stringsAsFactors=FALSE
  )
  pb_meta$age_norm <- tolower(pb_meta$age)
  pb_meta$sex_norm <- tolower(pb_meta$sex)
  saveRDS(list(pb_all=pb_all,pb_meta=pb_meta),pb_cache)
  cat("Saved:", pb_cache, "\n")
  print(table(pb_meta$age_norm, pb_meta$condition, pb_meta$sex_norm))
}
