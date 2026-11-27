################################################################################
### Fig2_ED4_Pipeline.R
###
### Figure 2 and Extended Data Fig. 4: from the integrated eNE Seurat object to every panel.
### Run top to bottom:   source("Fig2_ED4_Pipeline.R")
### Needs the integrated eNE object (df_harmony; deposited in GEO) for steps load, deseq2, 2b and 2c: either already in
### the session as df_harmony or read from OBJ_RDS (settings below).
###
###   Step     Figure / purpose                 What it does                                         Reads
###   -------  -------------------------------  ---------------------------------------------------  ---------------------------
###   load     Object loading                   reads df_harmony from OBJ_RDS if it is not already    OBJ_RDS (optional)
###                                             in the session
###   deseq2   Differential expression          pseudobulk DESeq2 per sex x cell type for all 11      df_harmony
###                                             comparisons -> DESeq2_Results/ (+ volcanoes).
###                                             SKIPPED AUTOMATICALLY when all result tables exist.
###   2a       Fig 2a  UMAP                     cell type + age / stimulation / sex (full.umap.rpca.sketch) df_harmony
###   2b       Fig 2b  IEG dotplot              canonical IEGs x cell types, Young / Middle / Old NE  DESeq2 (per cell type), df_harmony
###   2e       Fig 2b bars + Fig 2e             young NE program at middle / old age; lost genes      DESeq2 (ALL tables)
###                                             split into failed induction vs basal shift
###   2c       Fig 2c  NEUROeSTIMator           neuronal activity score: violins, trajectories,      saved scores (or computes them)
###                                             delta heatmaps, all on one shared scale              + df_harmony
###   2d       Fig 2d  Venn                     young / middle / old NE gene sets                     DESeq2 (ALL tables)
###   2f       Fig 2f  GO lost / maintained     GO BP + KEGG of lost and maintained young-NE genes    DESeq2 (ALL tables)
###   sv_all   Strip volcanoes, every comp.     YEEYHC, MEEMHC, OEEOHC, OEEYEE, MEEYEE, OEEMEE,      DESeq2 (ALL tables)
###                                             OHCYHC, MHCYHC; no log2FC ceiling
###   sv_mee   Strip volcano MEEYEE             middle NE vs young NE, basal-aging genes (MHCYHC) grey DESeq2 (ALL tables)
###   sv_oee   Strip volcano OEEYEE             old NE vs young NE, basal-aging genes (OHCYHC) grey   DESeq2 (ALL tables)
###   ed4cd    ED Fig 4c,d  YEEYHC pathways     GO BP + KEGG of young NE genes                        DESeq2 (ALL tables)
###
### Comparisons (DESeq2, each sex separately, ~ group, genes with >= 10 counts in >= 2 samples, BH padj < 0.05):
###   NE vs HC within age  YEEYHC, MEEMHC, OEEOHC      HC across age  OHCYHC, MHCYHC, OHCMHC
###   NE across age        OEEYEE, MEEYEE, OEEMEE      cross          OEEYHC, MEEYHC
###
### Every step runs in its own environment and writes to its own subfolder of <PATH_TO_OUTPUT>/Figures/Figure2/
### (2a_UMAP, 2b_IEG_Dotplot, 2c_NEUROeSTIMator, 2d_Venn, 2e_StackedBar_FailedVsBasal, 2f_Pathway_LostMaintained,
### StripVolcano_AllComparisons, StripVolcano_MEEYEE_ExclMHCYHC, StripVolcano_OEEYEE_ExclOHCYHC, ED4cd_Pathway_YEEYHC, Logs).
### DESeq2 tables go to <PATH_TO_OUTPUT>/DESeq2_Results/. If a figure step fails, the others still run; a summary at the
### end lists what finished and what failed. Steps deseq2, 2e, 2f and ed4cd skip files that already exist; the other
### figure steps redraw their (fast) plots. To run only some steps, set RUN below. Progress is printed as it goes and
### saved to Figures/Figure2/Logs/Figure2_<time>.log.
################################################################################

################################################################################
### SETTINGS
################################################################################

RUN <- c("load", "deseq2", "2a", "2b", "2e", "2c", "2d", "2f", "sv_all", "sv_mee", "sv_oee", "ed4cd")
# Strip volcanoes (sv_all / sv_mee / sv_oee): how each plot's log2FC axis is sized when Y_FIXED is NULL
#   "sig" = fits that plot's own most extreme SIGNIFICANT gene (padj < 0.05); non-significant genes beyond it are drawn
#           at the axis edge, so a few extreme low-count genes no longer flatten everything else
#   "all" = fits the most extreme gene of any kind (previous behaviour)
SV_Y_FIT <- "all"
SV_Y_QUANTILE <- 1      # 1 = axis reaches the single most extreme gene (no gene is moved); below 1, outliers sit at the edge
SV_Y_TRANS <- "break"        # "pseudo_log" = every gene at its true log2FC on a stretched axis: linear near 0, compressed
                             #   for large |log2FC|, so the bulk spreads out while outliers stay visible at their real values
                             # "linear"     = plain linear axis (previous behaviour)
SV_PSEUDO_SIGMA <- 0.5       # pseudo_log: axis is ~linear within +/- this log2FC; smaller = more stretch near 0
                             # "break"      = broken axis: linear up to a break point, then a short compressed segment
                             #   that still reaches the most extreme gene; the break is marked on the axis with //
SV_BREAK_Q      <- 0.99      # break: break point = this quantile of |log2FC| among significant genes (x 1.1)
SV_BREAK_SHRINK <- 0.15      # break: beyond the break, 1 log2FC takes this fraction of the space it would take below it
sv_ymax <- function(lfc, sig) {                                  # symmetric half-range for one strip volcano
  use <- if (identical(SV_Y_FIT, "sig") && any(sig, na.rm=TRUE)) lfc[sig %in% TRUE] else lfc
  y <- stats::quantile(abs(use), SV_Y_QUANTILE, na.rm=TRUE, names=FALSE) * 1.05
  n_edge <- sum(abs(lfc) > y, na.rm=TRUE)
  if (n_edge) cat(sprintf("  axis +/-%.2f; %d gene(s) beyond it drawn at the edge\n", y, n_edge))
  y
}
sv_break_at <- function(lfc, sig, ymax) {                        # break point for SV_Y_TRANS = "break" (NA = no break needed)
  if (!identical(SV_Y_TRANS, "break")) return(NA_real_)
  use <- if (any(sig, na.rm=TRUE)) lfc[sig %in% TRUE] else lfc
  B <- stats::quantile(abs(use), SV_BREAK_Q, na.rm=TRUE, names=FALSE) * 1.1
  if (!is.finite(B) || B >= ymax * 0.8) return(NA_real_)        # nothing far beyond the bulk: plain linear axis
  cat(sprintf("  broken axis at +/-%.2f (most extreme gene %.2f)\n", B, ymax / 1.05))
  B
}
sv_trans <- function(B=NA_real_) {
  if (identical(SV_Y_TRANS, "pseudo_log")) return(scales::pseudo_log_trans(sigma=SV_PSEUDO_SIGMA, base=2))
  if (identical(SV_Y_TRANS, "break") && is.finite(B)) {
    k <- SV_BREAK_SHRINK
    return(scales::trans_new("axis_break",
      transform=function(x) ifelse(abs(x) <= B, x, sign(x) * (B + (abs(x) - B) * k)),
      inverse  =function(x) ifelse(abs(x) <= B, x, sign(x) * (B + (abs(x) - B) / k))))
  }
  "identity"
}
sv_breaks <- function(ymax, B=NA_real_) {                        # axis ticks at real log2FC values
  if (identical(SV_Y_TRANS, "pseudo_log")) { b <- c(0.5, 1, 2, 4, 8, 16, 32); b <- b[b <= ymax]; return(c(-rev(b), 0, b)) }
  if (identical(SV_Y_TRANS, "break") && is.finite(B)) {
    inner <- scales::pretty_breaks(n=5)(c(-B, B)); inner <- inner[abs(inner) < B * 0.95]
    outer <- floor(ymax / 1.05)                                   # one tick near the most extreme gene
    return(sort(unique(c(inner, if (outer > B) c(-outer, outer)))))
  }
  scales::pretty_breaks(n=6)(c(-ymax, ymax))
}
sv_break_marks <- function(B=NA_real_) {                         # // marks on the y axis + faint guide lines at the break
  if (!identical(SV_Y_TRANS, "break") || !is.finite(B)) return(NULL)
  list(ggplot2::geom_hline(yintercept=c(-B, B), linetype="dotted", colour="grey55", linewidth=0.3),
       ggplot2::annotate("text", x=-Inf, y=c(-B, B), label="//", size=3.2, colour="grey20", hjust=0.5, vjust=0.5))
}
# -- Paths ----------------------------------------------------------------------------------------------------------
BASE    <- "PATH_TO_OUTPUT"            # DESeq2 tables and all figures are written here
OBJ_RDS    <- "PATH_TO_eNE_OBJECT.rds"     # integrated eNE Seurat object rebuilt from the GEO files (see README); ignored if df_harmony is in the session
COUNTS_H5  <- "PATH_TO_eNE_counts_bpcells.h5"  # counts file from GEO (eNE_counts_bpcells.h5);
                                               #   if left as is, eNE_counts_bpcells.h5 next to OBJ_RDS is used
if (grepl("^PATH_TO_", BASE)) stop("Set BASE at the top of this script")
DESEQ_DIR <- file.path(BASE, "DESeq2_Results")
FIG_ROOT  <- file.path(BASE, "Figures", "Figure2")   # every Figure 2 panel goes here, one subfolder per panel
dir.create(FIG_ROOT, recursive=TRUE, showWarnings=FALSE)

ALL_COMPS <- c("YEEYHC","MEEMHC","OEEOHC","OHCYHC","MHCYHC","OHCMHC","OEEYEE","MEEYEE","OEEMEE","OEEYHC","MEEYHC")
FIG_COMPS <- c("YEEYHC","MEEMHC","OEEOHC","OHCYHC","MHCYHC","OEEYEE","MEEYEE","OEEMEE")   # the tables Figure 2 reads

SECTION_INFO <- c(
  "load"   = "Load df_harmony",
  "deseq2" = "Step 1  pseudobulk DESeq2, all comparisons (skipped if results exist)",
  "2a"     = "Fig 2a  UMAP (cell type; age, stimulation state, sex)",
  "2b"     = "Fig 2b  IEG dotplot (Young / Middle / Old NE)",
  "2e"     = "Fig 2b bars + Fig 2e  young NE program at middle / old age (failed induction vs basal shift)",
  "2c"     = "Fig 2c  NEUROeSTIMator activity score (shared scale)",
  "2d"     = "Fig 2d  Venn of young / middle / old NE genes",
  "2f"     = "Fig 2f  GO / KEGG of lost and maintained genes",
  "sv_all" = "Strip volcanoes, every comparison (no ceiling)",
  "sv_mee" = "Strip volcano MEEYEE, basal-aging genes grey (no ceiling)",
  "sv_oee" = "Strip volcano OEEYEE, basal-aging genes grey (no ceiling)",
  "ed4cd"  = "ED Fig 4c,d  YEEYHC pathways (Sankey + dot)")
NEEDS_OBJECT <- c("load", "deseq2", "2a", "2b", "2c")    # steps that use df_harmony
CRITICAL     <- c("load", "deseq2")                # if these fail, stop (everything downstream depends on them)

LOG_DIR <- file.path(FIG_ROOT, "Logs"); dir.create(LOG_DIR, recursive=TRUE, showWarnings=FALSE)
LOG <- file.path(LOG_DIR, paste0("Figure2_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
say <- function(...) { m <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", ...); cat(m, "\n"); cat(m, "\n", file=LOG, append=TRUE) }
banner <- function(id) {
  line <- strrep("=", 90)
  cat("\n", line, "\n", sep=""); say(sprintf("STEP %s  |  %s", id, SECTION_INFO[[id]])); cat(line, "\n\n", sep="")
}
tables_ok <- function(comps) all(file.exists(file.path(DESEQ_DIR, paste0(rep(c("male","female"), each=length(comps)), "_ALL_", comps, ".csv"))))
deseq_complete <- function() tables_ok(ALL_COMPS)     # every comparison: step 1 is skipped
fig_ready      <- function() tables_ok(FIG_COMPS)     # what the figure steps need

################################################################################
### PRE-FLIGHT
################################################################################

say("Figure 2 pipeline started. Steps: ", paste(intersect(names(SECTION_INFO), RUN), collapse=", "), "   Log: ", LOG)
uses_obj <- intersect(RUN, NEEDS_OBJECT)
if ("deseq2" %in% RUN && deseq_complete()) uses_obj <- setdiff(uses_obj, "deseq2")
if (!length(setdiff(uses_obj, "load"))) uses_obj <- character(0)   # nothing downstream needs the object
NEED_OBJECT <- length(uses_obj) > 0
if (NEED_OBJECT && !exists("df_harmony") && is.null(OBJ_RDS))
  stop("Steps ", paste(uses_obj, collapse=", "), " need df_harmony: load it first or set OBJ_RDS")
say(if (NEED_OBJECT) paste0("df_harmony ", if (exists("df_harmony")) "in session" else paste("will be read from", OBJ_RDS),
                            "; used by: ", paste(uses_obj, collapse=", ")) else "df_harmony not needed for the selected steps")
say("DESeq2 results in ", DESEQ_DIR, ": ", if (deseq_complete()) "all comparisons present (step 1 will be skipped)" else
    paste0("missing ", paste(setdiff(paste0(rep(c("male","female"), each=length(ALL_COMPS)), "_ALL_", ALL_COMPS, ".csv"),
                                     list.files(DESEQ_DIR)), collapse=", "),
           if ("deseq2" %in% RUN) " (step 1 will fit only these)" else " (step 1 not in RUN)"))
if (!fig_ready() && !"deseq2" %in% RUN) stop("Figure 2 tables missing in ", DESEQ_DIR, ": add \"deseq2\" to RUN")

SECTIONS <- list()

################################################################################
### LOAD: df_harmony
################################################################################

SECTIONS[["load"]] <- function() {
  if (!NEED_OBJECT) { cat("df_harmony not needed for the selected steps: skipped\n"); return(invisible()) }
  if (!exists("df_harmony", envir=.GlobalEnv)) {
    if (grepl("^PATH_TO_", OBJ_RDS)) stop("Set OBJ_RDS to the integrated eNE object, or load it as df_harmony first")
    cat("Reading", OBJ_RDS, "...\n"); assign("df_harmony", readRDS(OBJ_RDS), envir=.GlobalEnv)
  } else cat("Using df_harmony already in the session\n")
  obj <- get("df_harmony", envir=.GlobalEnv)
  # point the counts at the downloaded GEO counts file (the object stores the file location it was saved with)
  ch5 <- if (!grepl("^PATH_TO_", COUNTS_H5)) COUNTS_H5 else file.path(dirname(OBJ_RDS), "eNE_counts_bpcells.h5")
  if (file.exists(ch5)) {
    cat("Counts from", ch5, "\n"); obj[["RNA"]]$counts <- BPCells::open_matrix_hdf5(ch5, group = "counts")
  }
  ok <- tryCatch({ invisible(as.matrix(LayerData(obj, assay="RNA", layer="counts")[1:3, 1:3])); TRUE }, error=function(e) FALSE)
  if (!ok) stop("Cannot read the counts: set COUNTS_H5 to eNE_counts_bpcells.h5 from GEO")
  if (!"data" %in% SeuratObject::Layers(obj[["RNA"]])) {     # the deposited object has counts only
    cat("Log-normalizing (NormalizeData) ...\n"); obj <- NormalizeData(obj, verbose=FALSE)
  }
  assign("df_harmony", obj, envir=.GlobalEnv)
}

################################################################################
### STEP deseq2: pseudobulk DESeq2. Skipped when every result table already exists in DESeq2_Results/.
################################################################################

SECTIONS[["deseq2"]] <- function() {
if (deseq_complete()) { cat("All", 2 * length(ALL_COMPS), "result tables exist in", DESEQ_DIR, "-> skipped\n"); return(invisible()) }
################################################################################
### DESeq2 PIPELINE - 36-SAMPLE LOGNORM INTEGRATED DATASET
### Merged IN cluster (PV IN + SST IN + VIP IN + Lamp5 IN = "IN")
### Includes MEEYHC (Middle EE vs Young HC), used in the buffering analysis
### Skip-if-exists throughout for safe reruns
### Writes DESeq2_Results/ and Volcano_Plots/ in BASE
################################################################################

library(conflicted)
library(BPCells)
library(Seurat)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyr)
library(tibble)

conflicts_prefer(base::intersect,     .quiet=TRUE)
conflicts_prefer(base::setdiff,       .quiet=TRUE)
conflicts_prefer(base::as.data.frame, .quiet=TRUE)
conflicts_prefer(dplyr::select,       .quiet=TRUE)
conflicts_prefer(dplyr::filter,       .quiet=TRUE)
conflicts_prefer(dplyr::arrange,      .quiet=TRUE)
conflicts_prefer(dplyr::mutate,       .quiet=TRUE)
conflicts_prefer(dplyr::pull,         .quiet=TRUE)
conflicts_prefer(DESeq2::vst,         .quiet=TRUE)
conflicts_prefer(dplyr::desc,         .quiet=TRUE)
conflicts_prefer(dplyr::rename,       .quiet=TRUE)

################################################################################
### SETUP
################################################################################

cat("Cell types:"); print(table(df_harmony$celltype))
cat("Sex:");        print(table(df_harmony$sex))
cat("Age:");        print(table(df_harmony$age))
cat("Condition:");  print(table(df_harmony$experimental.condition))

if (any(c("PV IN","SST IN","VIP IN","Lamp5 IN") %in% df_harmony$celltype))
  warning("IN subtype labels still present -- merge PV / SST / VIP / Lamp5 IN into \"IN\" first")

deseq_path   <- paste0(DESEQ_DIR, "/")
volcano_path <- paste0(file.path(BASE, "Volcano_Plots"), "/")
dir.create(deseq_path,   recursive=TRUE, showWarnings=FALSE)
dir.create(volcano_path, recursive=TRUE, showWarnings=FALSE)

PADJ_THRESHOLD <- 0.05
YOUNG_LABEL    <- "young"
OLD_LABEL      <- "old"

age_vals     <- unique(tolower(df_harmony@meta.data$age))
cat("Age labels in object:", paste(sort(age_vals),collapse=", "), "\n")
MIDDLE_LABEL <- age_vals[!age_vals %in% c("young","old")][1]
cat("Using middle age label:", MIDDLE_LABEL, "\n")
if (is.na(MIDDLE_LABEL)) stop("Could not detect middle-age label.")

################################################################################
### PSEUDOBULK AGGREGATION
################################################################################

cat("\n=== Aggregating pseudobulk counts ===\n")

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

cat("Pseudobulk samples:",nrow(pb_meta),"\n")
cat("Age x Condition x Sex:\n")
print(table(pb_meta$age_norm,pb_meta$condition,pb_meta$sex_norm))
celltypes_all <- sort(unique(pb_meta$celltype))
cat("Cell types:",paste(celltypes_all,collapse=", "),"\n")

################################################################################
### DESeq2 HELPER (with skip-if-exists)
################################################################################

run_deseq2 <- function(counts_mat, meta_df,
                        group_col, numerator, denominator,
                        celltype, comp_name,
                        out_dir=deseq_path, sex_label=current_sex) {
  # Skip if already done
  out_f <- file.path(out_dir,paste0(sex_label,"_",celltype,"_",comp_name,".csv"))
  if (file.exists(out_f)) {
    cat("  [skip existing]\n")
    return(read.csv(out_f,stringsAsFactors=FALSE))
  }

  ct_keep    <- meta_df$celltype==celltype
  sub_meta   <- meta_df[ct_keep,]
  sub_counts <- counts_mat[,ct_keep,drop=FALSE]

  grp_keep   <- sub_meta[[group_col]] %in% c(numerator,denominator)
  sub_meta   <- sub_meta[grp_keep,]
  sub_counts <- sub_counts[,grp_keep,drop=FALSE]
  if (ncol(sub_counts)==0) return(NULL)

  sub_meta$grp <- factor(sub_meta[[group_col]],
                          levels=c(denominator,numerator))
  grp_tab <- table(sub_meta$grp)
  cat("  Groups:",paste(names(grp_tab),grp_tab,sep="=",collapse=", "),"\n")

  if (length(grp_tab)<2 || any(grp_tab<2)) {
    cat("  Skipping -- need >= 2 samples per group\n"); return(NULL)
  }

  keep_genes <- rowSums(sub_counts>=10)>=2
  sub_counts <- sub_counts[keep_genes,,drop=FALSE]
  if (nrow(sub_counts)<10) {
    cat("  Skipping -- too few genes\n"); return(NULL)
  }

  tryCatch({
    dds    <- DESeqDataSetFromMatrix(sub_counts,sub_meta,design=~grp)
    dds    <- DESeq(dds,quiet=TRUE)
    res    <- results(dds,alpha=PADJ_THRESHOLD)
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::mutate(celltype=celltype,comparison=comp_name)
    n_sig <- sum(!is.na(res_df$padj)&res_df$padj<PADJ_THRESHOLD)
    cat("  DEGs:",n_sig,"\n")
    write.csv(res_df,out_f,row.names=FALSE)
    res_df
  }, error=function(e) { cat("  Error:",e$message,"\n"); NULL })
}

################################################################################
### CROSS-AGE HELPER (for OEEYHC, MEEYHC -- crosses condition and age)
################################################################################

run_cross_age <- function(pb_sex, pb_meta_sex, mask,
                           grp_col, ref_label, num_label,
                           comp_name, sex_label, celltypes) {
  all_f <- file.path(deseq_path,paste0(sex_label,"_ALL_",comp_name,".csv"))
  if (file.exists(all_f)) {
    cat("  ALL file exists -- loading\n")
    return(split(read.csv(all_f,stringsAsFactors=FALSE),
                 read.csv(all_f,stringsAsFactors=FALSE)$celltype))
  }

  sub_pb   <- pb_sex[,mask,drop=FALSE]
  sub_meta <- pb_meta_sex[mask,]
  sub_meta[[grp_col]] <- factor(sub_meta[[grp_col]],
                                  levels=c(ref_label,num_label))

  comp_results <- list()
  for (ct in celltypes) {
    ct_f <- file.path(deseq_path,paste0(sex_label,"_",ct,"_",comp_name,".csv"))
    if (file.exists(ct_f)) {
      cat(sprintf("  %-22s [skip]\n",ct))
      comp_results[[ct]] <- read.csv(ct_f,stringsAsFactors=FALSE)
      next
    }
    cat(sprintf("  %-22s: ",ct))
    ct_k  <- sub_meta$celltype==ct
    sm    <- sub_meta[ct_k,]
    sc    <- sub_pb[,ct_k,drop=FALSE]
    if (ncol(sc)==0) { cat("no samples\n"); next }
    gt    <- table(sm[[grp_col]])
    cat(paste(names(gt),gt,sep="=",collapse=" | ")," ")
    if (length(gt)<2||any(gt<2)) { cat("insufficient\n"); next }
    kg <- rowSums(sc>=10)>=2; sc <- sc[kg,,drop=FALSE]
    if (nrow(sc)<10) { cat("few genes\n"); next }
    res <- tryCatch({
      dds <- DESeqDataSetFromMatrix(sc,sm,design=as.formula(paste0("~",grp_col)))
      dds <- DESeq(dds,quiet=TRUE)
      r   <- results(dds,alpha=PADJ_THRESHOLD)
      as.data.frame(r) %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::mutate(celltype=ct,comparison=comp_name)
    },error=function(e){cat("Error:",e$message,"\n");NULL})
    if (!is.null(res)) {
      n_sig <- sum(!is.na(res$padj)&res$padj<PADJ_THRESHOLD)
      cat(sprintf("%d genes, %d sig\n",nrow(res),n_sig))
      write.csv(res,ct_f,row.names=FALSE)
      comp_results[[ct]] <- res
    } else cat("\n")
  }
  if (length(comp_results)>0) {
    all_df <- dplyr::bind_rows(comp_results)
    write.csv(all_df,all_f,row.names=FALSE)
    cat(sprintf("  Saved ALL_%s: %d rows\n",comp_name,nrow(all_df)))
  }
  comp_results
}

################################################################################
### VOLCANO PLOT FUNCTIONS
################################################################################

create_volcano <- function(deg_df, comparison_name, celltype,
                            out_dir, sex_label) {
  out_f <- file.path(out_dir,paste0(sex_label,"_",celltype,"_",
                                     comparison_name,"_Volcano.pdf"))
  if (file.exists(out_f)) return(invisible(NULL))
  df <- deg_df %>%
    dplyr::mutate(
      negLog10padj = -log10(padj),
      significant  = !is.na(padj)&padj<PADJ_THRESHOLD&abs(log2FoldChange)>0.5,
      color = dplyr::case_when(
        significant&log2FoldChange>0.5  ~ "Upregulated",
        significant&log2FoldChange< -0.5 ~ "Downregulated",
        TRUE ~ "Not Significant"))
  top_genes <- df %>% dplyr::filter(significant) %>%
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
    head(50) %>% dplyr::pull(gene)
  label_df <- df[df$gene %in% top_genes,]
  p <- ggplot(df,aes(x=log2FoldChange,y=negLog10padj)) +
    geom_point(aes(color=color),alpha=0.6,size=1.5) +
    scale_color_manual(values=c(Upregulated="red",Downregulated="blue",
                                 "Not Significant"="grey70"),name=NULL) +
    geom_vline(xintercept=c(-0.5,0.5),linetype="dashed",color="grey40") +
    geom_hline(yintercept=-log10(PADJ_THRESHOLD),linetype="dashed",color="grey40") +
    geom_text_repel(data=label_df,aes(label=gene),
                     size=2.5,box.padding=0.5,max.overlaps=50,seed=42) +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(face="bold",size=13)) +
    labs(title=paste0(celltype," (",sex_label,"): ",comparison_name),
         subtitle=paste0("Top 50 labeled | n=",nrow(df)),
         x="log2 Fold Change",y="-log10(adjusted p-value)")
  ggsave(out_f,p,width=10,height=8,device="pdf")
}

create_combined_volcano <- function(deg_list,comparison_name,out_dir,sex_label) {
  out_f <- file.path(out_dir,paste0("CombinedVolcano_",sex_label,"_",
                                     comparison_name,".pdf"))
  if (file.exists(out_f)) return(invisible(NULL))
  if (length(deg_list)==0) return(NULL)
  all_data <- dplyr::bind_rows(lapply(names(deg_list),function(ct)
    deg_list[[ct]] %>% dplyr::mutate(
      celltype=ct,
      sig=dplyr::case_when(
        !is.na(padj)&padj<PADJ_THRESHOLD&log2FoldChange>0 ~ "Up",
        !is.na(padj)&padj<PADJ_THRESHOLD&log2FoldChange<0 ~ "Down",
        TRUE ~ "NS"))))
  counts_df <- all_data %>% dplyr::filter(sig!="NS") %>%
    dplyr::group_by(celltype,sig) %>%
    dplyr::summarise(n=n(),.groups="drop") %>%
    tidyr::pivot_wider(names_from=sig,values_from=n,values_fill=0)
  if (!"Up"   %in% names(counts_df)) counts_df$Up   <- 0L
  if (!"Down" %in% names(counts_df)) counts_df$Down <- 0L
  ymax <- max(abs(all_data$log2FoldChange),na.rm=TRUE)
  y_pad <- ymax*0.20
  p <- ggplot() +
    geom_jitter(data=dplyr::filter(all_data,sig=="NS"),
                aes(x=celltype,y=log2FoldChange),
                width=0.25,alpha=0.4,size=1.2,color="gray70") +
    geom_jitter(data=dplyr::filter(all_data,sig=="Up"),
                aes(x=celltype,y=log2FoldChange),
                width=0.25,alpha=0.85,size=1.8,color="#D6604D") +
    geom_jitter(data=dplyr::filter(all_data,sig=="Down"),
                aes(x=celltype,y=log2FoldChange),
                width=0.25,alpha=0.85,size=1.8,color="#4393C3") +
    geom_hline(yintercept=0,linetype="dashed",color="black",linewidth=0.8) +
    geom_text(data=counts_df,
              aes(x=celltype,y=ymax+y_pad,label=paste0("",Up)),
              color="#D6604D",size=3.5,fontface="bold") +
    geom_text(data=counts_df,
              aes(x=celltype,y=-ymax-y_pad,label=paste0("",Down)),
              color="#4393C3",size=3.5,fontface="bold") +
    coord_cartesian(ylim=c(-ymax-y_pad*2.5,ymax+y_pad*2.5),clip="off") +
    theme_classic(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=9),
          plot.title=element_text(face="bold",size=13),
          plot.subtitle=element_text(size=10,color="grey40"),
          panel.grid.major.y=element_line(color="gray90",linewidth=0.3)) +
    labs(title=paste0(comparison_name," - ",tools::toTitleCase(sex_label)),
         subtitle="Red=Up (padj<0.05) | Blue=Down | Gray=NS",
         x="Cell Type",y="log2 Fold Change")
  ggsave(out_f,p,
         width=max(8,length(unique(all_data$celltype))*0.5+4),
         height=7,device="pdf")
  cat("  Saved combined volcano:",comparison_name,"\n")
  p
}

################################################################################
### COMPARISONS
################################################################################

EEvHC_comparisons <- list(
  list(name="YEEYHC",  age=YOUNG_LABEL,  num="EE", den="HC"),
  list(name="OEEOHC",  age=OLD_LABEL,    num="EE", den="HC"),
  list(name="MEEMHC",  age=MIDDLE_LABEL, num="EE", den="HC")
)
HC_age_comparisons <- list(
  list(name="OHCYHC", num=OLD_LABEL,    den=YOUNG_LABEL),
  list(name="OHCMHC", num=OLD_LABEL,    den=MIDDLE_LABEL),
  list(name="MHCYHC", num=MIDDLE_LABEL, den=YOUNG_LABEL)
)
EE_age_comparisons <- list(
  list(name="OEEYEE", num=OLD_LABEL,    den=YOUNG_LABEL),
  list(name="OEEMEE", num=OLD_LABEL,    den=MIDDLE_LABEL),
  list(name="MEEYEE", num=MIDDLE_LABEL, den=YOUNG_LABEL)
)

################################################################################
### MAIN LOOP
################################################################################

for (current_sex in c("male","female")) {

  cat(sprintf("\n%s\n### SEX: %s\n%s\n",
              strrep("#",80),toupper(current_sex),strrep("#",80)))

  sex_keep    <- pb_meta$sex_norm==current_sex
  pb_sex      <- pb_all[,sex_keep,drop=FALSE]
  pb_meta_sex <- pb_meta[sex_keep,]
  cat("Samples:",nrow(pb_meta_sex),"\n")
  celltypes <- sort(unique(pb_meta_sex$celltype))

  all_results <- list()

  # -- EE vs HC (within age) -------------------------------------------------
  for (comp in EEvHC_comparisons) {
    comp_name <- comp$name
    all_f     <- file.path(deseq_path,paste0(current_sex,"_ALL_",comp_name,".csv"))
    cat("\n===",comp_name,"===\n")
    if (file.exists(all_f)) {
      cat("  ALL file exists -- skip\n")
      all_results[[comp_name]] <- split(
        read.csv(all_f,stringsAsFactors=FALSE),
        read.csv(all_f,stringsAsFactors=FALSE)$celltype)
      next
    }
    age_keep    <- pb_meta_sex$age_norm==tolower(comp$age)
    pb_age      <- pb_sex[,age_keep,drop=FALSE]
    pb_meta_age <- pb_meta_sex[age_keep,]
    comp_results <- list()
    for (ct in celltypes) {
      cat("  Celltype:",ct,"\n")
      res <- run_deseq2(pb_age,pb_meta_age,"condition",
                         comp$num,comp$den,ct,comp_name)
      if (!is.null(res)) comp_results[[ct]] <- res
    }
    if (length(comp_results)>0) {
      all_results[[comp_name]] <- comp_results
      write.csv(dplyr::bind_rows(comp_results),all_f,row.names=FALSE)
    }
  }

  # -- HC age comparisons ----------------------------------------------------
  hc_keep    <- pb_meta_sex$condition=="HC"
  pb_hc      <- pb_sex[,hc_keep,drop=FALSE]
  pb_meta_hc <- pb_meta_sex[hc_keep,]
  for (comp in HC_age_comparisons) {
    comp_name <- comp$name
    all_f     <- file.path(deseq_path,paste0(current_sex,"_ALL_",comp_name,".csv"))
    cat("\n===",comp_name,"===\n")
    if (file.exists(all_f)) {
      cat("  ALL file exists -- skip\n")
      all_results[[comp_name]] <- split(
        read.csv(all_f,stringsAsFactors=FALSE),
        read.csv(all_f,stringsAsFactors=FALSE)$celltype)
      next
    }
    comp_results <- list()
    for (ct in celltypes) {
      cat("  Celltype:",ct,"\n")
      res <- run_deseq2(pb_hc,pb_meta_hc,"age_norm",
                         tolower(comp$num),tolower(comp$den),ct,comp_name)
      if (!is.null(res)) comp_results[[ct]] <- res
    }
    if (length(comp_results)>0) {
      all_results[[comp_name]] <- comp_results
      write.csv(dplyr::bind_rows(comp_results),all_f,row.names=FALSE)
    }
  }

  # -- EE age comparisons ----------------------------------------------------
  ee_keep    <- pb_meta_sex$condition=="EE"
  pb_ee      <- pb_sex[,ee_keep,drop=FALSE]
  pb_meta_ee <- pb_meta_sex[ee_keep,]
  for (comp in EE_age_comparisons) {
    comp_name <- comp$name
    all_f     <- file.path(deseq_path,paste0(current_sex,"_ALL_",comp_name,".csv"))
    cat("\n===",comp_name,"===\n")
    if (file.exists(all_f)) {
      cat("  ALL file exists -- skip\n")
      all_results[[comp_name]] <- split(
        read.csv(all_f,stringsAsFactors=FALSE),
        read.csv(all_f,stringsAsFactors=FALSE)$celltype)
      next
    }
    comp_results <- list()
    for (ct in celltypes) {
      cat("  Celltype:",ct,"\n")
      res <- run_deseq2(pb_ee,pb_meta_ee,"age_norm",
                         tolower(comp$num),tolower(comp$den),ct,comp_name)
      if (!is.null(res)) comp_results[[ct]] <- res
    }
    if (length(comp_results)>0) {
      all_results[[comp_name]] <- comp_results
      write.csv(dplyr::bind_rows(comp_results),all_f,row.names=FALSE)
    }
  }

  # -- OEEYHC: Old EE vs Young HC (cross-age, cross-condition) -------------
  cat("\n=== OEEYHC: Old EE vs Young HC ===\n")
  oeeyhc_mask <- (pb_meta_sex$age_norm==YOUNG_LABEL & pb_meta_sex$condition=="HC") |
                  (pb_meta_sex$age_norm==OLD_LABEL   & pb_meta_sex$condition=="EE")
  pb_meta_sex_oeeyhc <- pb_meta_sex
  pb_meta_sex_oeeyhc$grp_oeeyhc <- ifelse(
    pb_meta_sex_oeeyhc$age_norm==OLD_LABEL,"OEE","YHC")
  res_oeeyhc <- run_cross_age(
    pb_sex, pb_meta_sex_oeeyhc, oeeyhc_mask,
    "grp_oeeyhc","YHC","OEE","OEEYHC",current_sex,celltypes)
  if (length(res_oeeyhc)>0) all_results[["OEEYHC"]] <- res_oeeyhc

  # -- MEEYHC: Middle EE vs Young HC (cross-age, cross-condition) ----------
  # NEW: same structure as OEEYHC but using MIDDLE_LABEL instead of OLD_LABEL
  cat("\n=== MEEYHC: Middle EE vs Young HC ===\n")
  meeyhc_mask <- (pb_meta_sex$age_norm==YOUNG_LABEL   & pb_meta_sex$condition=="HC") |
                  (pb_meta_sex$age_norm==MIDDLE_LABEL & pb_meta_sex$condition=="EE")
  pb_meta_sex_meeyhc <- pb_meta_sex
  pb_meta_sex_meeyhc$grp_meeyhc <- ifelse(
    pb_meta_sex_meeyhc$age_norm==MIDDLE_LABEL,"MEE","YHC")
  res_meeyhc <- run_cross_age(
    pb_sex, pb_meta_sex_meeyhc, meeyhc_mask,
    "grp_meeyhc","YHC","MEE","MEEYHC",current_sex,celltypes)
  if (length(res_meeyhc)>0) all_results[["MEEYHC"]] <- res_meeyhc

  # -- VOLCANOS --------------------------------------------------------------
  cat("\n=== Generating volcano plots for",toupper(current_sex),"===\n")
  all_comp_names <- c(sapply(EEvHC_comparisons,"[[","name"),
                       sapply(HC_age_comparisons,"[[","name"),
                       sapply(EE_age_comparisons,"[[","name"))
  volcano_plots <- list()
  for (comp_name in all_comp_names) {
    if (!comp_name %in% names(all_results)) next
    for (ct in names(all_results[[comp_name]]))
      create_volcano(all_results[[comp_name]][[ct]],comp_name,ct,
                      volcano_path,current_sex)
    p <- create_combined_volcano(all_results[[comp_name]],comp_name,
                                  volcano_path,current_sex)
    if (!is.null(p)) volcano_plots[[comp_name]] <- p
  }
  if (length(volcano_plots)>0) {
    pdf_f <- file.path(volcano_path,
      paste0(current_sex,"_AllComparisons_Combined_Volcano.pdf"))
    if (!file.exists(pdf_f)) {
      pdf(pdf_f,width=max(10,length(celltypes)*0.5+4),
          height=7*length(volcano_plots))
      for (p in volcano_plots) print(p)
      dev.off()
    }
  }

  # -- DEG SUMMARY -----------------------------------------------------------
  all_sex_df <- dplyr::bind_rows(lapply(all_results,dplyr::bind_rows))
  if (nrow(all_sex_df)>0) {
    summary_wide <- all_sex_df %>%
      dplyr::filter(!is.na(padj)&padj<PADJ_THRESHOLD) %>%
      dplyr::group_by(celltype,comparison) %>%
      dplyr::summarise(n_total=n(),n_up=sum(log2FoldChange>0),
                        n_down=sum(log2FoldChange<0),.groups="drop") %>%
      tidyr::pivot_wider(names_from=comparison,values_from=n_total,values_fill=0)
    cat("\nDEG Summary --",toupper(current_sex),":\n")
    print(summary_wide)
    write.csv(summary_wide,
              file.path(deseq_path,paste0(current_sex,"_DEG_Summary.csv")),
              row.names=FALSE)
  }
  cat("\n--- Completed",toupper(current_sex),"---\n")
}

cat("\n=== DONE ===\n")
cat("Results:", deseq_path, "\n")
cat("\nComparisons now available:\n")
print(list.files(deseq_path,pattern="^male_ALL_",full.names=FALSE))
invisible(NULL)
}

################################################################################
### STEP 2a  Fig. 2a UMAP (cell type; age, stimulation state and sex insets)
################################################################################

SECTIONS[["2a"]] <- function() {
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
suppressPackageStartupMessages({
  library(conflicted)
  conflicts_prefer(dplyr::select, .quiet = TRUE)
  conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicts_prefer(base::setdiff, .quiet = TRUE)
})

fig_date <- format(Sys.Date(), "%Y-%m-%d")
out_dir  <- file.path(FIG_ROOT, "2a_UMAP")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
cat("Output:", out_dir, "\n\n")

################################################################################
### CELL TYPE COLOR PALETTE (used in every figure)
################################################################################

ct_levels <- sort(unique(df_harmony@meta.data$celltype))

CT_COLORS <- c(
  "DG ExNs"            ="#E53935",
  "CA1 ExNs"           ="#CC0000",
  "CA2 ExNs"           ="#8E24AA",
  "CA3 ExNs"           ="#5E35B1",
  "Sub ExNs"           ="#80D8FF",
  "IN"                 ="#00897B",
  "Astrocytes"         ="#43A047",
  "Oligodendrocytes"   ="#1E88E5",
  "OPC"                ="#00ACC1",
  "Microglia"          ="#7CB342",  
  "Endothelial"        ="#FB8C00",
  "Pericyte"           ="#6D4C41",
  "APC"                ="#546E7A",
  "Stromal"            ="#27408B", 
  "Macrophage"         ="#C62828",
  "Cajal-Retzius Cell" ="#FDD835",
  "T Cell"             ="#757575",
  "Ambiguous"          ="#BDBDBD")

# Fill fallback for any unexpected cell types
missing_cts <- setdiff(ct_levels, names(CT_COLORS))
if (length(missing_cts)>0) {
  CT_COLORS[missing_cts] <- "#999999"
  cat("Fallback color assigned to:", paste(missing_cts,collapse=", "), "\n")
}
CT_COLORS <- CT_COLORS[ct_levels]

################################################################################
### UMAPs
################################################################################

cat("\n=== Generating UMAPs ===\n")

umap_theme <- theme_classic(base_size=11) +
  theme(
    plot.title      = element_text(face="bold", size=13, hjust=0.5),
    legend.text     = element_text(size=7.5),
    legend.key.size = unit(0.4,"cm"),
    axis.line       = element_blank(),
    axis.ticks      = element_blank(),
    axis.text       = element_blank(),
    axis.title      = element_blank()
  )

# Extract UMAP coordinates
umap_df <- as.data.frame(
  Embeddings(df_harmony, reduction="full.umap.rpca.sketch")
)
names(umap_df) <- c("UMAP1","UMAP2")
umap_df <- cbind(umap_df, df_harmony@meta.data[
  rownames(umap_df),
  c("celltype","age","sex","experimental.condition")])

# Shuffle so no group is systematically on top
set.seed(42)
umap_df <- umap_df[sample(nrow(umap_df)), ]

# -- A: Cell type ------------------------------------------------------
p_ct <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=celltype)) +
  geom_point(size=0.3, alpha=1.0) +
  scale_color_manual(values=CT_COLORS, name=NULL) +
  guides(color=guide_legend(
    override.aes=list(size=3, alpha=1), ncol=1)) +
  umap_theme +
  labs(title="Cell type")

# -- B: Age -----------------------------------------------------------
# Age: Young=blue (Oligo #1565C0), Middle=teal, Old=red (CA1 #CC0000)
AGE_COLS <- c(young="#1565C0", middleage="#00897B", old="#FF1744")
# Normalize age labels to lowercase
umap_df$age_norm <- tolower(umap_df$age)

p_age <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=age_norm)) +
  geom_point(size=0.3, alpha=1.0) +
  scale_color_manual(
    values = AGE_COLS,
    labels = c(young="Young",middle="Middle",old="Old"),
    name   = NULL) +
  guides(color=guide_legend(override.aes=list(size=3, alpha=1))) +
  umap_theme +
  labs(title="Age",
       subtitle="Old = red | Young = blue | Middle = teal")

# -- C: Condition -----------------------------------------------------
# Condition: EE=red (CA1 ExNs red), HC=blue (Oligo blue)
COND_COLS <- c(EE="#FF1744", HC="#1565C0")
umap_df$cond_norm <- toupper(umap_df$experimental.condition)

p_cond <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=cond_norm)) +
  geom_point(size=0.3, alpha=1.0) +
  scale_color_manual(
    values = COND_COLS,
    labels = c(EE="Enriched Environment",HC="Home Cage"),
    name   = NULL) +
  guides(color=guide_legend(override.aes=list(size=3, alpha=1))) +
  umap_theme +
  labs(title="Condition",
       subtitle="EE = red | HC = blue")

# -- D: Sex ------------------------------------------------------------
# Sex: female=red (CA1 ExNs #CC0000), male=blue (Oligo #1565C0)
SEX_COLS <- c(male="#1565C0", female="#FF1744",
               Male="#1565C0", Female="#FF1744")

p_sex <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=sex)) +
  geom_point(size=0.3, alpha=1.0) +
  scale_color_manual(
    values = SEX_COLS,
    name   = NULL) +
  guides(color=guide_legend(override.aes=list(size=3, alpha=1))) +
  umap_theme +
  labs(title="Sex",
       subtitle="Female = red | Male = blue")

# -- Save each panel separately (PDF + PNG) ----------------------------
save_both <- function(p, stem, w, h) {
  ggsave(paste0(stem,".pdf"), p, width=w, height=h, device="pdf")
  ggsave(paste0(stem,".png"), p, width=w, height=h, device="png",
         dpi=300, bg="white")
  cat("  Saved:", basename(stem), ".pdf + .png\n")
}

save_both(p_ct,   file.path(out_dir,paste0(fig_date,"_UMAP_CellType")),   w=9, h=7)
save_both(p_age,  file.path(out_dir,paste0(fig_date,"_UMAP_Age")),        w=7, h=6)
save_both(p_cond, file.path(out_dir,paste0(fig_date,"_UMAP_Condition")),  w=7, h=6)
save_both(p_sex,  file.path(out_dir,paste0(fig_date,"_UMAP_Sex")),        w=7, h=6)

# -- Combined 2x2 panel ------------------------------------------------
p_combined <- (p_ct | p_age) / (p_cond | p_sex) +
  patchwork::plot_annotation(
    title    = "EE dataset - hippocampal single-cell atlas",
    subtitle = paste0(
      format(ncol(df_harmony), big.mark=","), " cells  |  ",
      length(unique(umap_df$celltype)), " cell types  |  ",
      "IN = merged PV IN + SST IN + VIP IN + Lamp5 IN"),
    theme    = theme(
      plot.title    = element_text(face="bold", size=14, hjust=0.5),
      plot.subtitle = element_text(size=9, color="grey40", hjust=0.5)
    )
  )

save_both(p_combined,
          file.path(out_dir,paste0(fig_date,"_UMAP_Combined_2x2")),
          w=14, h=12)

# -- Cell type UMAP with labels on the plot ----------------------------
# Compute centroid per cell type for labeling
centroids <- umap_df %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(UMAP1=median(UMAP1), UMAP2=median(UMAP2), .groups="drop")

p_labeled <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=celltype)) +
  geom_point(size=0.12, alpha=0.5) +
  ggrepel::geom_text_repel(
    data       = centroids,
    aes(x=UMAP1, y=UMAP2, label=celltype),
    color      = "grey10",
    size       = 2.8,
    fontface   = "bold",
    box.padding = 0.4,
    max.overlaps = 30,
    seed       = 42,
    segment.color = "grey60",
    segment.size  = 0.3
  ) +
  scale_color_manual(values=CT_COLORS, guide="none") +
  umap_theme +
  labs(title="Cell types - EE dataset (IN merged)",
       subtitle=paste0(format(ncol(df_harmony),big.mark=","),
                        " cells  |  ",
                        length(ct_levels), " cell types"))

save_both(p_labeled,
          file.path(out_dir,paste0(fig_date,"_UMAP_CellType_Labeled")),
          w=9, h=7)

cat("\n=== DONE ===\nOutputs:", out_dir, "\n")
invisible(NULL)
}

################################################################################
### STEP 2b  Fig. 2b IEG dotplot
################################################################################

SECTIONS[["2b"]] <- function() {
################################################################################
### IEG_Dotplot_EE_Ages
###
### Dotplot of canonical IEGs across Young / Middle / Old EE responses.
### Male and female separately.
###
### Design:
###   Y-axis:  IEG gene
###   X-axis:  cell type (all cell types including merged IN)
###   Color:   log2FC from per-cell-type DESeq2 (EE vs HC)
###   Size:    % of cells in EE condition expressing that gene
###   Facets:  age (Young EE / Middle EE / Old EE)
###
### Requires df_harmony in session for percent expressed.
### Reads the per-cell-type DESeq2 tables (step deseq2).
################################################################################

library(Seurat)
library(BPCells)
library(dplyr)
library(tidyr)
library(ggplot2)

suppressPackageStartupMessages({
  library(conflicted)
  conflicts_prefer(dplyr::select,   .quiet = TRUE)
  conflicts_prefer(dplyr::filter,   .quiet = TRUE)
  conflicts_prefer(base::intersect, .quiet = TRUE)
  conflicts_prefer(base::union,     .quiet = TRUE)
})

################################################################################
### GENE LIST
################################################################################

# Focused IEG list for main figure
IEGS <- c(
  # Clearly compromised with age -- tell the aging story
  "Fosl1",    # most dramatic delta
  "Fos",      # canonical, clearly reduced
  "Fosb",     # canonical AP-1, reduced
  "Egr3",     # moderate reduction, 5 sig cell types young

  # Middle ground -- partial compromise
  "Arc",      # THE canonical ARG, everyone knows it
  "Homer1",   # synaptic plasticity, modest reduction

  # Maintained -- balance the story
  "Nr4a1",    # flat across ages, very canonical
  "Junb",     # nearly flat, canonical AP-1
  "Nptx2",    # maintained/slight increase, well known
  "Npas4"     # increases with age, canonical ARG
)

################################################################################
### PATHS
################################################################################

deseq_path <- paste0(DESEQ_DIR, "/")
stopifnot("Run step deseq2 first (no MEEMHC per-cell-type files in DESeq2_Results)" =
  length(list.files(deseq_path, pattern="^(male|female)_.*_MEEMHC\\.csv$")) > 0)

fig_date   <- format(Sys.Date(),"%Y-%m-%d")
out_dir    <- file.path(FIG_ROOT, "2b_IEG_Dotplot")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
cat("Output:", out_dir, "\n\n")

PADJ_THR <- 0.05
COMPS    <- c("YEEYHC","MEEMHC","OEEOHC")
AGE_LABELS <- c(YEEYHC="Young EE", MEEMHC="Middle EE", OEEOHC="Old EE")

# Cell type order: excitatory -> IN -> glia -> other
CT_ORDER <- c(
  "DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs",
  "IN",
  "Astrocytes","Oligodendrocytes","OPC","Microglia",
  "Endothelial","Pericyte","APC","Stromal",
  "Cajal-Retzius Cell","Macrophage","Epithelial","T Cell"
)

################################################################################
### VERIFY df_harmony AND DETECT AGE / SEX LABELS
################################################################################

if (!exists("df_harmony")) stop("df_harmony not in session.")

# Detect middle age label
age_vals     <- unique(tolower(df_harmony@meta.data$age))
MIDDLE_NORM  <- age_vals[!age_vals %in% c("young","old")][1]
cat("Age labels detected:", paste(sort(age_vals), collapse=", "), "\n")
cat("Middle label:", MIDDLE_NORM, "\n")

AGE_MAP <- c(YEEYHC="young", MEEMHC=MIDDLE_NORM, OEEOHC="old")

################################################################################
### STEP 1: PERCENT EXPRESSED FROM df_harmony
################################################################################

cat("\n=== Step 1: Computing % expressed ===\n")

iegs_found <- intersect(IEGS, rownames(df_harmony))
cat(sprintf("IEGs in object: %d / %d\n", length(iegs_found), length(IEGS)))
missing <- setdiff(IEGS, iegs_found)
if (length(missing)>0) cat("Missing:", paste(missing, collapse=", "), "\n")

cat("Extracting counts (BPCells -> dgCMatrix)...\n")
ieg_counts <- as(
  GetAssayData(df_harmony, assay="RNA", layer="counts")[iegs_found, ],
  "dgCMatrix"
)
cat(sprintf("  Matrix: %d genes x %d cells\n", nrow(ieg_counts), ncol(ieg_counts)))

meta <- df_harmony@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::select(barcode, celltype,
                age_raw = age,
                sex     = sex,
                cond    = experimental.condition) %>%
  dplyr::mutate(age_norm = tolower(age_raw),
                sex_norm = tolower(sex))

pct_rows <- list()

for (current_sex in c("male","female")) {
  for (comp in COMPS) {
    age_val  <- AGE_MAP[comp]
    cells_ee <- meta$barcode[
      meta$sex_norm == current_sex &
      meta$age_norm == age_val &
      meta$cond     == "EE"]

    if (length(cells_ee) == 0) {
      cat(sprintf("  No EE cells: %s %s (age='%s')\n",
                  current_sex, comp, age_val)); next
    }

    sub_bc     <- intersect(cells_ee, colnames(ieg_counts))
    sub_counts <- ieg_counts[, sub_bc, drop=FALSE]
    sub_meta   <- meta[meta$barcode %in% sub_bc, ]

    for (ct in sort(unique(sub_meta$celltype))) {
      ct_bc  <- sub_meta$barcode[sub_meta$celltype == ct]
      ct_mat <- sub_counts[, intersect(ct_bc, colnames(sub_counts)), drop=FALSE]
      if (ncol(ct_mat) == 0) next
      pct <- rowMeans(ct_mat > 0) * 100
      for (g in iegs_found) {
        pct_rows[[paste(current_sex, comp, ct, g)]] <- data.frame(
          sex=current_sex, comparison=comp, celltype=ct,
          gene=g, pct_expr=pct[[g]], n_cells=ncol(ct_mat),
          stringsAsFactors=FALSE)
      }
    }
    cat(sprintf("  %s %s: %d EE cells, %d cell types\n",
                current_sex, comp, length(sub_bc),
                length(unique(sub_meta$celltype))))
  }
}

pct_df <- dplyr::bind_rows(pct_rows)
write.csv(pct_df,
          file.path(out_dir, paste0(fig_date,"_IEG_PctExpressed.csv")),
          row.names=FALSE)
cat(sprintf("Pct expressed: %d rows\n\n", nrow(pct_df)))

################################################################################
### STEP 2: LOG2FC FROM DESEQ2 CSVs
################################################################################

cat("=== Step 2: Loading log2FC from DESeq2 ===\n")

lfc_rows <- list()
for (current_sex in c("male","female")) {
  for (comp in COMPS) {
    # Load per-cell-type files (not the ALL file, to get per-CT log2FC)
    files <- list.files(deseq_path,
                         pattern=paste0("^",current_sex,"_.*_",comp,"\\.csv$"),
                         full.names=TRUE)
    files <- files[!grepl("_ALL_", files)]
    for (f in files) {
      df <- tryCatch(read.csv(f, stringsAsFactors=FALSE),
                     error=function(e) NULL)
      if (is.null(df) || nrow(df)==0) next
      if (!all(c("gene","celltype","log2FoldChange","padj") %in% names(df))) next
      sub <- df %>%
        dplyr::filter(gene %in% IEGS) %>%
        dplyr::select(gene, celltype, log2FoldChange, padj) %>%
        dplyr::mutate(sex=current_sex, comparison=comp)
      if (nrow(sub)>0)
        lfc_rows[[paste(current_sex, comp, basename(f))]] <- sub
    }
    cat(sprintf("  %s %s: %d files\n", current_sex, comp, length(files)))
  }
}

lfc_df <- dplyr::bind_rows(lfc_rows) %>%
  dplyr::mutate(sig = !is.na(padj) & padj < PADJ_THR)
cat(sprintf("LFC data: %d rows\n\n", nrow(lfc_df)))

################################################################################
### STEP 3: JOIN
################################################################################

dot_df <- dplyr::inner_join(lfc_df, pct_df,
                              by=c("sex","comparison","celltype","gene")) %>%
  dplyr::mutate(
    # Cap at 3 -- most IEG LFCs fall within this range
    # Tighter cap = more contrast between low and high expressors
    lfc_cap = pmax(pmin(log2FoldChange, 3), -3),
    age_lbl = factor(AGE_LABELS[comparison],
                      levels=c("Young EE","Middle EE","Old EE"))
  )

write.csv(dot_df,
          file.path(out_dir, paste0(fig_date,"_IEG_Dotplot_data.csv")),
          row.names=FALSE)
cat(sprintf("Joined data: %d rows\n", nrow(dot_df)))

# Gene order: most strongly induced in Young EE -> top of plot
gene_order <- dot_df %>%
  dplyr::filter(comparison=="YEEYHC", sig, log2FoldChange>0) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(mean_lfc=mean(log2FoldChange, na.rm=TRUE), .groups="drop") %>%
  dplyr::arrange(dplyr::desc(mean_lfc)) %>%
  dplyr::pull(gene)
gene_order <- unique(c(gene_order, IEGS))
gene_order <- intersect(gene_order, unique(dot_df$gene))

################################################################################
### STEP 4: PLOT
################################################################################

make_dotplot <- function(df, title, out_file) {
  if (nrow(df)==0) { cat("  No data for:", title, "\n"); return(invisible(NULL)) }

  # Enforce cell type and gene ordering
  ct_present  <- intersect(CT_ORDER, unique(df$celltype))
  ct_extra    <- setdiff(unique(df$celltype), CT_ORDER)
  ct_ordered  <- c(ct_present, sort(ct_extra))

  df <- df %>%
    dplyr::mutate(
      gene     = factor(gene,     levels=rev(gene_order)),
      celltype = factor(celltype, levels=ct_ordered)
    )

  p <- ggplot(df, aes(x=celltype, y=gene, color=lfc_cap, size=pct_expr)) +
    # All dots (grey if NS)
    geom_point(data=df[!df$sig, ],
               aes(x=celltype, y=gene, size=pct_expr),
               color="grey88", alpha=0.7) +
    # Significant dots colored by LFC
    geom_point(data=df[df$sig, ],
               aes(x=celltype, y=gene, color=lfc_cap, size=pct_expr),
               alpha=0.92) +
    # Tighter range (-3 to 3) + stronger color endpoints for contrast
    scale_color_gradient2(
      low  = "#1E88E5",   # Bright electric blue -- low / repression
      mid  = "grey97",    # Near-white midpoint
      high = "#E53935",   # deep crimson -- strong induction.  "#1E88E5", mid="white", high="#CC0000"
      midpoint =0,
      limits   =c(-3,3),
      oob      =scales::squish,
      name     ="log2FC\n(EE vs HC)"
    ) +
    # Larger max size + finer breaks for better visual separation
    scale_size_continuous(
      name   ="% EE cells\nexpressing",
      range  =c(0.3, 10),
      limits =c(0,100),
      breaks =c(5,10,25,50,75)
    ) +
    facet_wrap(~age_lbl, nrow=1) +
    theme_classic(base_size=9) +
    theme(
      plot.title       = element_text(face="bold", size=12, hjust=0.5),
      plot.subtitle    = element_text(size=7.5, color="grey40", hjust=0.5),
      axis.text.x      = element_text(angle=45, hjust=1, size=7.5),
      axis.text.y      = element_text(size=8, face="italic"),
      axis.title       = element_blank(),
      strip.text       = element_text(face="bold", size=9),
      legend.position  = "right",
      legend.text      = element_text(size=7.5),
      panel.grid.major = element_line(color="grey92", linewidth=0.3),
      panel.spacing    = unit(0.5, "lines")
    ) +
    guides(
      color = guide_colorbar(barwidth=0.8, barheight=5, title.position="top"),
      size  = guide_legend(title.position="top",
                            override.aes=list(color="grey50"))
    ) +
    labs(
      title    = title,
      subtitle = paste0(
        "Color = log2FC (EE vs HC, per cell type)  |  ",
        "Size = % EE cells expressing gene\n",
        "Grey = not significant (padj >= 0.05)  |  ",
        "Genes ordered by Young EE induction")
    )

  n_genes <- length(unique(df$gene))
  n_cts   <- length(unique(df$celltype))
  ggsave(out_file, p,
         width  = max(12, n_cts * 0.5 + 4),
         height = max(5,  n_genes * 0.38 + 3),
         device = "pdf")
  cat("  Saved:", basename(out_file), "\n")
}

for (current_sex in c("male","female")) {
  sub <- dot_df %>% dplyr::filter(sex == current_sex)
  if (nrow(sub)==0) next
  make_dotplot(
    sub,
    title    = paste0("IEG expression across the lifespan -- ",
                       tools::toTitleCase(current_sex)),
    out_file = file.path(out_dir,
                          paste0(fig_date,"_IEG_Dotplot_",
                                 tools::toTitleCase(current_sex),".pdf"))
  )
}

cat("\n=== DONE ===\nOutputs:", out_dir, "\n")

invisible(NULL)
}

################################################################################
### STEP 2e  Fig. 2b bars + Fig. 2e
################################################################################

SECTIONS[["2e"]] <- function() {
################################################################################
### Fig2E_StackedBar_FailedVsBasal
###
### Mirrored stacked bar (males left, females right) of what happens to the young NE program (YEEYHC) at
### middle age (MEEMHC) and old age (OEEOHC), with LOST genes split by why they were lost:
###
###   Maintained Up / Down   sig in YEEYHC and sig in the same direction at the older age
###   Flipped                sig in YEEYHC, sig in the opposite direction
###   Lost: failed induction           sig in YEEYHC, not sig at the older age, and NO basal change
###                                    in the same direction (aged HC vs young HC)          -> reds
###   Lost: basal shift (de-repression)  sig in YEEYHC, not sig at the older age, and the resting level
###                                    already moved that way (aged HC vs young HC, same direction)
###                                    Up = basally de-repressed, Down = basally co-suppressed -> purples
###   Emergent Up / Down     sig at the older age, not in YEEYHC (optional, INCLUDE_EMERGENT)
###
### Failed-induction segments share one hue family (red, dark = Up, light = Down); basal-shift segments
### share another (purple, dark = Up, light = Down), for every age and cell type.
### Same definitions as Figure 3 (basal = MHCYHC for middle, OHCYHC for old).
###
### Outputs (one shared x scale for middle and old):
###   _Counts.pdf / .png     absolute gene counts
###   _Pct.pdf / .png        % of each bar
###   _Simple_<Middle|Old>   Maintained / Flipped / Lost only (no split, no emergent): Fig 2b style
### Reads the DESeq2 tables (step deseq2). Skip-if-exists.
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
suppressPackageStartupMessages({
  library(conflicted)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(base::intersect, .quiet=TRUE); conflicts_prefer(base::setdiff, .quiet=TRUE)
  conflicts_prefer(base::union, .quiet=TRUE)
})

################################################################################
### PATHS AND SETTINGS
################################################################################

deseq_dir <- paste0(DESEQ_DIR, "/")
stopifnot("Run step deseq2 first" = file.exists(file.path(deseq_dir, "male_ALL_MEEMHC.csv")))
fig_date  <- format(Sys.Date(),"%Y-%m-%d")
out_dir   <- paste0(file.path(FIG_ROOT, "2e_StackedBar_FailedVsBasal"), "/")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
cat("Output:", out_dir, "\n\n")

PADJ_THR <- 0.05
INCLUDE_EMERGENT <- TRUE        # FALSE = YEEYHC fate only
COMPS   <- c(MEEMHC="MHCYHC", OEEOHC="OHCYHC")    # older-age NE comparison -> its basal (HC) aging comparison
COMP_LABELS <- c(MEEMHC="Middle age (MEEMHC)", OEEOHC="Old age (OEEOHC)")

CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN",
              "Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal",
              "Cajal-Retzius Cell","Macrophage","Epithelial","T Cell")

SEG_ORDER <- c("Maintained Up","Maintained Down","Flipped",
               "Failed induction Up","Failed induction Down",
               "Basal shift Up","Basal shift Down",
               "Emergent Up","Emergent Down")
SEG_COLS <- c("Maintained Up"         = "#238B45",
              "Maintained Down"       = "#74C476",
              "Flipped"               = "#969696",
              "Failed induction Up"   = "#B2182B",   # failed induction: red family
              "Failed induction Down" = "#F4A582",
              "Basal shift Up"        = "#54278F",   # basal shift: purple family
              "Basal shift Down"      = "#BCBDDC",
              "Emergent Up"           = "#2171B5",
              "Emergent Down"         = "#9ECAE1")
SEG_LABELS <- c("Maintained Up"         = "Maintained (Up)",
                "Maintained Down"       = "Maintained (Down)",
                "Flipped"               = "Flipped direction",
                "Failed induction Up"   = "Lost: failed induction (Up)",
                "Failed induction Down" = "Lost: failed induction (Down)",
                "Basal shift Up"        = "Lost: basal de-repression (Up, rising at rest)",
                "Basal shift Down"      = "Lost: basal co-suppression (Down, falling at rest)",
                "Emergent Up"           = "Emergent (Up, new at this age)",
                "Emergent Down"         = "Emergent (Down, new at this age)")
if (!INCLUDE_EMERGENT) SEG_ORDER <- setdiff(SEG_ORDER, c("Emergent Up","Emergent Down"))

should_skip <- function(f) { if (file.exists(f)) { cat("  Skipping:", basename(f), "\n"); return(TRUE) }; FALSE }
save_both <- function(p, stem, w, h) {
  ggsave(paste0(stem, ".pdf"), p, width=w, height=h, device="pdf")
  ggsave(paste0(stem, ".png"), p, width=w, height=h, device="png", dpi=300, bg="white")
  cat("Saved:", basename(stem), "[pdf+png]\n") }

################################################################################
### CLASSIFY
################################################################################

sig_dir <- function(d, dir) d$gene[!is.na(d$padj) & d$padj < PADJ_THR & (if (dir > 0) d$log2FoldChange > 0 else d$log2FoldChange < 0)]

classify_ct <- function(yee_df, comp_df, hc_df, ct) {
  y <- yee_df  %>% dplyr::filter(celltype == ct)
  d <- comp_df %>% dplyr::filter(celltype == ct)
  h <- hc_df   %>% dplyr::filter(celltype == ct)
  y_up <- sig_dir(y, 1); y_dn <- sig_dir(y, -1); y_sig <- union(y_up, y_dn)
  d_up <- sig_dir(d, 1); d_dn <- sig_dir(d, -1); d_sig <- union(d_up, d_dn)
  h_up <- sig_dir(h, 1); h_dn <- sig_dir(h, -1)
  lost_up <- setdiff(y_up, d_sig); lost_dn <- setdiff(y_dn, d_sig)
  data.frame(segment=c("Maintained Up","Maintained Down","Flipped",
                       "Failed induction Up","Failed induction Down","Basal shift Up","Basal shift Down",
                       "Emergent Up","Emergent Down"),
             n_genes=c(length(intersect(y_up, d_up)), length(intersect(y_dn, d_dn)),
                       length(intersect(y_up, d_dn)) + length(intersect(y_dn, d_up)),
                       length(setdiff(lost_up, h_up)), length(setdiff(lost_dn, h_dn)),
                       length(intersect(lost_up, h_up)), length(intersect(lost_dn, h_dn)),
                       length(setdiff(d_up, y_sig)), length(setdiff(d_dn, y_sig))),
             stringsAsFactors=FALSE) %>% dplyr::filter(segment %in% SEG_ORDER)
}

################################################################################
### BUILD DATA
################################################################################

all_rows <- list()
for (sx in c("male","female")) {
  cat(sprintf("\n### %s ###\n", toupper(sx)))
  load_df <- function(comp) { f <- file.path(deseq_dir, paste0(sx, "_ALL_", comp, ".csv"))
    if (!file.exists(f)) { cat("  Not found:", comp, "\n"); return(NULL) }; read.csv(f, stringsAsFactors=FALSE) }
  yee <- load_df("YEEYHC"); if (is.null(yee)) next
  for (comp_nm in names(COMPS)) {
    comp_df <- load_df(comp_nm); hc_df <- load_df(COMPS[[comp_nm]])
    if (is.null(comp_df) || is.null(hc_df)) next
    for (ct in sort(unique(yee$celltype))) {
      cl <- classify_ct(yee, comp_df, hc_df, ct)
      cl$celltype <- ct; cl$comparison <- comp_nm; cl$sex <- sx
      all_rows[[paste(sx, ct, comp_nm)]] <- cl
      g <- setNames(cl$n_genes, cl$segment)
      cat(sprintf("  %-20s %-7s maintained=%d flipped=%d failed=%d basal=%d%s\n", ct, comp_nm,
                  sum(g[c("Maintained Up","Maintained Down")]), g[["Flipped"]],
                  sum(g[c("Failed induction Up","Failed induction Down")]), sum(g[c("Basal shift Up","Basal shift Down")]),
                  if (INCLUDE_EMERGENT) sprintf(" emergent=%d", sum(g[c("Emergent Up","Emergent Down")])) else ""))
    }
  }
}
plot_df <- dplyr::bind_rows(all_rows) %>%
  dplyr::mutate(segment=factor(segment, levels=SEG_ORDER),
                comparison=factor(COMP_LABELS[comparison], levels=COMP_LABELS))
write.csv(plot_df, file.path(out_dir, paste0(fig_date, "_Fig2E_Counts_BothSexes.csv")), row.names=FALSE)
ct_levels <- rev(c(intersect(CT_ORDER, unique(plot_df$celltype)), sort(setdiff(unique(plot_df$celltype), CT_ORDER))))

################################################################################
### PLOT
################################################################################

build_rect <- function(df, segs) {
  df %>% dplyr::filter(segment %in% segs, n_genes > 0) %>%
    dplyr::mutate(segment=factor(as.character(segment), levels=segs)) %>%
    dplyr::group_by(sex, celltype, comparison) %>% dplyr::arrange(segment, .by_group=TRUE) %>%
    dplyr::mutate(abs_end=cumsum(n_genes), abs_start=abs_end - n_genes,
                  xmin=ifelse(sex == "male", -abs_end, abs_start), xmax=ifelse(sex == "male", -abs_start, abs_end)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(celltype=factor(celltype, levels=ct_levels), y_ctr=as.integer(celltype))
}

make_plot <- function(rect_df, segs, cols, labs_map, x_max, x_breaks, x_labels, title, x_lab, total_fmt) {
  BAR_H <- 0.38
  totals <- rect_df %>% dplyr::group_by(sex, celltype, comparison, y_ctr) %>%
    dplyr::summarise(total=sum(n_genes), tip=ifelse(sex[1] == "male", -max(abs_end), max(abs_end)), .groups="drop")
  x_lim <- x_max * 1.22
  ggplot() +
    geom_rect(data=rect_df, aes(xmin=xmin, xmax=xmax, ymin=y_ctr - BAR_H, ymax=y_ctr + BAR_H, fill=segment),
              color="white", linewidth=0.2) +
    geom_text(data=totals, aes(x=tip + ifelse(sex == "male", -x_max * 0.02, x_max * 0.02), y=y_ctr,
                               label=total_fmt(total), hjust=ifelse(sex == "male", 1, 0)), size=2.4, fontface="bold") +
    geom_vline(xintercept=0, linewidth=0.5, color="grey20") +
    annotate("text", x=-x_lim * 0.55, y=length(ct_levels) + 0.7, label="Males", fontface="bold", size=3.5) +
    annotate("text", x= x_lim * 0.55, y=length(ct_levels) + 0.7, label="Females", fontface="bold", size=3.5) +
    scale_fill_manual(values=cols[segs], breaks=segs, labels=labs_map[segs], name=NULL, drop=FALSE,
                      guide=guide_legend(ncol=2)) +
    scale_x_continuous(limits=c(-x_lim, x_lim), breaks=x_breaks, labels=x_labels, expand=c(0, 0)) +
    scale_y_continuous(breaks=seq_along(ct_levels), labels=ct_levels, expand=expansion(add=c(0.6, 1.0))) +
    facet_wrap(~comparison, ncol=1) +
    theme_classic(base_size=10) +
    theme(plot.title=element_text(face="bold", size=12, hjust=0.5), axis.text.y=element_text(size=9),
          axis.line.y=element_blank(), axis.ticks.y=element_blank(),
          strip.text=element_text(face="bold", size=10, color="white"), strip.background=element_rect(fill="grey35", color=NA),
          legend.position="bottom", legend.text=element_text(size=8),
          panel.grid.major.x=element_line(color="grey92", linewidth=0.3), panel.spacing=unit(0.6, "lines")) +
    labs(title=title, x=x_lab, y=NULL)
}
H <- max(8, length(ct_levels) * 0.45 + 5)

# Counts (one x scale for both ages)
stem <- file.path(out_dir, paste0(fig_date, "_Fig2E_FailedVsBasal_Counts"))
if (!should_skip(paste0(stem, ".pdf"))) {
  r <- build_rect(plot_df, SEG_ORDER)
  x_max <- max(r$abs_end)
  save_both(make_plot(r, SEG_ORDER, SEG_COLS, SEG_LABELS, x_max, pretty(c(-x_max, x_max), 6), function(x) abs(x),
                      "Young NE program at middle and old age: why genes are lost", "Number of genes", function(x) x),
            stem, 12, H)
}

# Percent of each bar
stem <- file.path(out_dir, paste0(fig_date, "_Fig2E_FailedVsBasal_Pct"))
if (!should_skip(paste0(stem, ".pdf"))) {
  pct <- plot_df %>% dplyr::group_by(sex, celltype, comparison) %>%
    dplyr::mutate(n_genes=ifelse(sum(n_genes) > 0, 100 * n_genes / sum(n_genes), 0)) %>% dplyr::ungroup()
  r <- build_rect(pct, SEG_ORDER)
  save_both(make_plot(r, SEG_ORDER, SEG_COLS, SEG_LABELS, 100, seq(-100, 100, 25), function(x) paste0(abs(x), "%"),
                      "Young NE program at middle and old age (% of each bar)", "% of genes",
                      function(x) sprintf("%.0f%%", x)), stem, 12, H)
}

# Simple version (Fig 2b style): Maintained / Flipped / Lost, one file per age, shared x scale
SIMPLE <- c("Maintained Up","Maintained Down","Lost Up","Lost Down","Flipped")
SIMPLE_COLS <- c("Maintained Up"="#238B45", "Maintained Down"="#74C476", "Lost Up"="#B2182B", "Lost Down"="#F4A582", "Flipped"="#4FA3C7")
simple_df <- plot_df %>% dplyr::filter(!grepl("^Emergent", segment)) %>%
  dplyr::mutate(segment=dplyr::case_when(segment %in% c("Failed induction Up","Basal shift Up") ~ "Lost Up",
                                         segment %in% c("Failed induction Down","Basal shift Down") ~ "Lost Down",
                                         TRUE ~ as.character(segment))) %>%
  dplyr::group_by(sex, celltype, comparison, segment) %>% dplyr::summarise(n_genes=sum(n_genes), .groups="drop")
r_all <- build_rect(simple_df, SIMPLE); x_max <- max(r_all$abs_end)
for (cp in levels(plot_df$comparison)) {
  stem <- file.path(out_dir, paste0(fig_date, "_Fig2b_Simple_", sub(" .*", "", cp)))
  if (should_skip(paste0(stem, ".pdf"))) next
  r <- r_all %>% dplyr::filter(comparison == cp)
  save_both(make_plot(r, SIMPLE, SIMPLE_COLS, setNames(SIMPLE, SIMPLE), x_max, pretty(c(-x_max, x_max), 6), function(x) abs(x),
                      paste0("Changes to the young ARG program: ", cp), "Number of genes", function(x) x),
            stem, 10, max(6, length(ct_levels) * 0.32 + 3))
}

cat("\n=== DONE ===\nOutputs:", out_dir, "\n")
cat("Lost = sig in YEEYHC, not sig at the older age. Failed induction = no same-direction basal change (aged HC vs\n")
cat("young HC); basal shift = the resting level already moved that way (de-repressed if Up, co-suppressed if Down).\n")

invisible(NULL)
}

################################################################################
### STEP 2c  Fig. 2c NEUROeSTIMator
################################################################################

SECTIONS[["2c"]] <- function() {
################################################################################
### NEUROeSTIMator PIPELINE -- EE/HC DATASET (36 Samples)
### Male + Female, Young + Middle-aged + Old, EE vs HC
### All violins, trajectories and heatmaps on one shared scale.
### Uses the merged IN cluster (PV IN + SST IN + VIP IN + Lamp5 IN -> IN)
################################################################################

library(Seurat)
library(BPCells)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(patchwork)
library(conflicted)
conflicts_prefer(dplyr::select,   .quiet=TRUE)
conflicts_prefer(dplyr::filter,   .quiet=TRUE)
conflicts_prefer(dplyr::summarise,.quiet=TRUE)
conflicts_prefer(base::intersect, .quiet=TRUE)
conflicts_prefer(BPCells::as.data.frame)
conflicts_prefer(SummarizedExperiment::as.data.frame)
conflicts_prefer(base::as.data.frame)

################################################################################
### PARAMETERS
################################################################################

MINICONDA_PATH <- "PATH_TO_MINICONDA"     # conda installation with the NEUROeSTIMator environment
CONDA_ENV      <- "neuroestimator"

OUT_DIR <- file.path(FIG_ROOT, "2c_NEUROeSTIMator")
dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
RAW_SCORES <- file.path(OUT_DIR, "NEUROeSTIMator_RawScores_AllCells.csv")    # per-cell scores; reused if present

# -- Neuronal cell types (merged IN) ----------------------------------------
NEURON_CELLTYPES <- c(
  "CA1 ExNs","CA2 ExNs","CA3 ExNs",
  "DG ExNs","Sub ExNs","IN")  # IN = merged PV/SST/VIP/Lamp5

# -- Cell type colors (used in every figure) ---------------------------------
CT_COLORS <- c(
  "DG ExNs"            ="#E53935",
  "CA1 ExNs"           ="#CC0000",
  "CA2 ExNs"           ="#8E24AA",
  "CA3 ExNs"           ="#5E35B1",
  "Sub ExNs"           ="#80D8FF",
  "IN"                 ="#00897B",
  "Astrocytes"         ="#43A047",
  "Oligodendrocytes"   ="#1E88E5",
  "OPC"                ="#00ACC1",
  "Microglia"          ="#7CB342",
  "Endothelial"        ="#FB8C00",
  "Pericyte"           ="#6D4C41",
  "APC"                ="#546E7A",
  "Stromal"            ="#27408B",
  "Macrophage"         ="#C62828",
  "Cajal-Retzius Cell" ="#FDD835",
  "T Cell"             ="#757575",
  "Ambiguous"          ="#BDBDBD")

# -- Age x condition factor levels + colors ----------------------------------
AGE_COND_LEVELS <- c(
  "Young_HC","Young_EE",
  "Middleage_HC","Middleage_EE",
  "Old_HC","Old_EE")

# Age x condition colors -- derived from standard age palette
# Young=#1565C0, Middleage=#00897B, Old=#FF1744
# HC = lighter/desaturated, EE = full saturation
AGE_COND_COLORS <- c(
  "Young_HC"    ="#90CAF9","Young_EE"    ="#1565C0",
  "Middleage_HC"="#80CBC4","Middleage_EE"="#00897B",
  "Old_HC"      ="#FF8A80","Old_EE"      ="#FF1744")

AGE_COND_LABELS <- c(
  "Young_HC"    ="Young HC",    "Young_EE"    ="Young EE",
  "Middleage_HC"="Mid-aged HC", "Middleage_EE"="Mid-aged EE",
  "Old_HC"      ="Old HC",      "Old_EE"      ="Old EE")

################################################################################
### STEP 0: VERIFY METADATA + MERGE IN + CONFIGURE PYTHON
################################################################################

cat("=== Step 0: Metadata verification ===\n")
if (!exists("df_harmony")) stop("df_harmony not in session.")

cat("Age values:\n");       print(table(df_harmony@meta.data$age))
cat("Condition values:\n"); print(table(df_harmony@meta.data$experimental.condition))
cat("Sex values:\n");       print(table(df_harmony@meta.data$sex))
cat("Cell types:\n");
print(sort(table(df_harmony@meta.data$celltype),decreasing=TRUE))

# Merge IN subtypes -> "IN" if not already done
OLD_IN <- c("PV IN","SST IN","VIP IN","Lamp5 IN")
in_present <- intersect(OLD_IN,unique(df_harmony@meta.data$celltype))
if (length(in_present)>0) {
  cat(sprintf("Merging %s -> IN\n",paste(in_present,collapse="+")))
  df_harmony@meta.data$celltype[
    df_harmony@meta.data$celltype %in% in_present] <- "IN"
  cat("After merge:\n")
  print(sort(table(df_harmony@meta.data$celltype),decreasing=TRUE))
} else {
  cat("IN already merged or IN subtypes not found -- using existing celltype\n")
}

# UMAP reduction
all_reds <- names(df_harmony@reductions)
n_cells  <- ncol(df_harmony)
umap_reduction <- NULL
for (r in c("full.umap.rpca.sketch","umap","umap.rpca","umap.harmony")) {
  if (r %in% all_reds) {
    emb <- df_harmony@reductions[[r]]@cell.embeddings
    if (nrow(emb)==n_cells) { umap_reduction <- r; break }
  }
}
if (is.null(umap_reduction)) stop("No UMAP reduction found")
cat("Using UMAP:",umap_reduction,"\n")



################################################################################
### STEP 1: CREATE age_condition COLUMN
################################################################################

cat("\n=== Step 1: Creating age_condition ===\n")

standardized_age <- dplyr::case_when(
  df_harmony@meta.data$age=="young"     ~ "Young",
  df_harmony@meta.data$age=="Middleage" ~ "Middleage",
  df_harmony@meta.data$age=="old"       ~ "Old",
  TRUE ~ as.character(df_harmony@meta.data$age))

df_harmony@meta.data$age_condition <- factor(
  paste0(standardized_age,"_",df_harmony@meta.data$experimental.condition),
  levels=AGE_COND_LEVELS)

cat("Age x condition:\n")
print(table(df_harmony@meta.data$age_condition))
cat("Age x condition x sex:\n")
print(table(df_harmony@meta.data$age_condition,df_harmony@meta.data$sex))

################################################################################
### STEP 2: EXTRACT RAW COUNTS
################################################################################

cat("\n=== Step 2: Extracting raw counts ===\n")

if (inherits(df_harmony@assays$RNA,"Assay5")) {
  cat("Assay5 -- extracting counts layer\n")
  expr_matrix <- LayerData(df_harmony,assay="RNA",layer="counts")
} else {
  cat("Standard Assay -- extracting counts slot\n")
  expr_matrix <- GetAssayData(df_harmony,assay="RNA",slot="counts")
}
cat("Matrix:",nrow(expr_matrix),"genes x",ncol(expr_matrix),"cells\n")
gene_names <- rownames(df_harmony)

################################################################################
### STEPS 3-4: SCORES (reused if already computed)
################################################################################

if (file.exists(RAW_SCORES)) {
  cat("\nReusing NEUROeSTIMator scores:", RAW_SCORES, "\n")
  activity_scores_df <- as.matrix(read.csv(RAW_SCORES, row.names=1))
} else {
cat("\n=== Configuring Python environment ===\n")
Sys.setenv(RETICULATE_MINICONDA_PATH=MINICONDA_PATH)
reticulate::use_condaenv(CONDA_ENV,required=TRUE)
reticulate::py_run_string(
  "import h5py; print('h5py version:', h5py.__version__)")
library(neuroestimator)
cat("NEUROeSTIMator loaded\n")

################################################################################
### STEP 3: TEST ON 1000 CELL SUBSET
################################################################################

cat("\n=== Step 3: Test run (1000 cells) ===\n")
set.seed(42)
test_cells  <- sample(colnames(df_harmony),1000)
test_matrix <- as.matrix(expr_matrix[,test_cells])
if (is.null(rownames(test_matrix))) rownames(test_matrix) <- gene_names

cat("Running NEUROeSTIMator on 1000 test cells...\n")
test_scores <- neuroestimator(test_matrix,species="mmusculus")
cat("Score range:",round(min(test_scores[,1]),3),"to",
    round(max(test_scores[,1]),3),"\n")
if (all(is.na(test_scores[,1]))) stop("Test FAILED -- check conda environment")
cat("Test PASSED\n")

################################################################################
### STEP 4: FULL DATASET IN BATCHES
################################################################################

cat("\n=== Step 4: Running on full dataset ===\n")
cat("Total cells:",ncol(expr_matrix),"\n")

batch_size   <- 50000
cell_indices <- seq_len(ncol(expr_matrix))
batches      <- split(cell_indices,ceiling(cell_indices/batch_size))
cat("Batches:",length(batches),"\n\n")

all_scores <- vector("list",length(batches))
for (i in seq_along(batches)) {
  cat("Batch",i,"of",length(batches),
      "| cells",min(batches[[i]]),"-",max(batches[[i]]),"...\n")
  batch_scores <- tryCatch({
    batch_dense <- as.matrix(expr_matrix[,batches[[i]]])
    if (is.null(rownames(batch_dense))) rownames(batch_dense) <- gene_names
    neuroestimator(batch_dense,species="mmusculus")
  },error=function(e) { cat("  ERROR batch",i,":",e$message,"\n"); NULL })
  if (!is.null(batch_scores)) {
    all_scores[[i]] <- batch_scores
    cat("  Done | range:",
        round(min(batch_scores[,1]),3),"to",
        round(max(batch_scores[,1]),3),"\n")
  }
  gc()
}

all_scores         <- all_scores[!sapply(all_scores,is.null)]
activity_scores_df <- do.call(rbind,all_scores)
cat("\nTotal scores:",nrow(activity_scores_df),"\n")
cat("Range:",round(min(activity_scores_df[,1]),3),"to",
    round(max(activity_scores_df[,1]),3),"\n")

write.csv(activity_scores_df, RAW_SCORES)
cat("Saved raw scores\n")
}

################################################################################
### STEP 5: ADD SCORES TO METADATA
################################################################################

cat("\n=== Step 5: Adding scores to metadata ===\n")
matched_scores <- activity_scores_df[colnames(df_harmony),1]
df_harmony@meta.data$neuro_activity <- matched_scores
cat("Scores added\n")
print(summary(df_harmony@meta.data$neuro_activity))

################################################################################
### STEP 6: NEURON-ONLY DATAFRAME
################################################################################

cat("\n=== Step 6: Building neuron dataframe ===\n")
available_neurons <- intersect(NEURON_CELLTYPES,
                                unique(df_harmony@meta.data$celltype))
cat("Neuron types found:",paste(available_neurons,collapse=", "),"\n")

neuron_df <- df_harmony@meta.data %>%
  tibble::rownames_to_column("cell_barcode") %>%
  dplyr::filter(celltype %in% available_neurons) %>%
  dplyr::mutate(
    standardized_age=dplyr::case_when(
      age=="young"     ~ "Young",
      age=="Middleage" ~ "Middleage",
      age=="old"       ~ "Old",
      TRUE             ~ as.character(age))) %>%
  dplyr::select(cell_barcode,celltype,age,experimental.condition,
                sex,age_condition,neuro_activity,standardized_age) %>%
  dplyr::mutate(
    celltype     =factor(celltype,
                          levels=intersect(NEURON_CELLTYPES,available_neurons)),
    age_condition=factor(age_condition,levels=AGE_COND_LEVELS),
    age          =factor(standardized_age,
                          levels=c("Young","Middleage","Old")),
    sex          =factor(tolower(sex),levels=c("female","male"))) %>%
  dplyr::select(-standardized_age)

cat("Neuron cells:",nrow(neuron_df),"\n")
cat("IN cells:",sum(neuron_df$celltype=="IN"),"\n")
cat("Cells per age_condition x sex:\n")
print(table(neuron_df$age_condition,neuron_df$sex))

# Neuronal colors
NEURON_COLORS <- CT_COLORS[intersect(NEURON_CELLTYPES,names(CT_COLORS))]

# One y range for every violin (both sexes, all cell types, all groups)
Y_LIM <- stats::quantile(neuron_df$neuro_activity, c(0.001, 0.999), na.rm=TRUE)
cat("Shared violin y range:", round(Y_LIM, 3), "\n")

################################################################################
### STEP 7: UMAP COLORED BY ACTIVITY SCORE (ggplot, no DimPlot fog)
################################################################################

cat("\n=== Step 7: UMAP plots ===\n")

# Build plot dataframe
umap_coords <- as.data.frame(
  df_harmony@reductions[[umap_reduction]]@cell.embeddings)
colnames(umap_coords) <- c("UMAP_1","UMAP_2")
umap_df_plot <- data.frame(
  UMAP_1        =umap_coords$UMAP_1,
  UMAP_2        =umap_coords$UMAP_2,
  neuro_activity=df_harmony@meta.data$neuro_activity,
  celltype      =df_harmony@meta.data$celltype,
  condition     =df_harmony@meta.data$experimental.condition,
  age           =standardized_age,
  sex           =tolower(df_harmony@meta.data$sex),
  age_condition =df_harmony@meta.data$age_condition,
  stringsAsFactors=FALSE) %>%
  dplyr::arrange(neuro_activity)  # high activity plotted on top

umap_theme <- theme_classic(base_size=11) +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.text=element_blank(),axis.title=element_text(size=9,color="grey60"),
        plot.title=element_text(face="bold",size=12,hjust=0.5),
        plot.subtitle=element_text(size=8.5,color="grey40",hjust=0.5),
        plot.background=element_rect(fill="white",color=NA),
        panel.background=element_rect(fill="white",color=NA))

save_both <- function(p,stem,w,h,dpi=600) {
  ggsave(paste0(stem,".pdf"),p,width=w,height=h,device="pdf")
  ggsave(paste0(stem,".png"),p,width=w,height=h,device="png",
         dpi=dpi,bg="white")
  cat("  Saved:",basename(stem),".pdf + .png\n")
}

# 7a: All cells activity score
p_all <- ggplot(umap_df_plot,
                aes(x=UMAP_1,y=UMAP_2,color=neuro_activity)) +
  geom_point(size=0.05,shape=16,alpha=1.0,stroke=0) +
  scale_color_gradient2(low="lightgrey",mid="#1565C0",high="#FF1744",
                         midpoint=median(umap_df_plot$neuro_activity,
                                          na.rm=TRUE),
                         name="Activity\nscore",na.value="grey90") +
  umap_theme +
  labs(title="NEUROeSTIMator Activity Score -- All Cells",
       subtitle="Higher = more transcriptionally active",
       x="UMAP 1",y="UMAP 2")
save_both(p_all,file.path(OUT_DIR,"UMAP_ActivityScore_AllCells"),w=9,h=7)

# 7b: Split by condition
for (cond in c("EE","HC")) {
  pd <- umap_df_plot %>% dplyr::mutate(
    act_plot=ifelse(condition==cond,neuro_activity,NA_real_)) %>%
    dplyr::arrange(is.na(act_plot),act_plot)
  p <- ggplot() +
    geom_point(data=pd%>%dplyr::filter(is.na(act_plot)),
               aes(x=UMAP_1,y=UMAP_2),
               color="lightgrey",size=0.04,shape=16,alpha=0.5,stroke=0) +
    geom_point(data=pd%>%dplyr::filter(!is.na(act_plot)),
               aes(x=UMAP_1,y=UMAP_2,color=act_plot),
               size=0.08,shape=16,alpha=1.0,stroke=0) +
    scale_color_gradient2(low="lightgrey",mid="#1565C0",high="#FF1744",
                           midpoint=median(umap_df_plot$neuro_activity,
                                            na.rm=TRUE),
                           name="Activity\nscore") +
    umap_theme +
    labs(title=paste0("NEUROeSTIMator -- ",cond),
         subtitle=paste0(cond," cells highlighted"),
         x="UMAP 1",y="UMAP 2")
  save_both(p,file.path(OUT_DIR,paste0("UMAP_ActivityScore_",cond)),w=9,h=7)
}

# 7c: Split by age (faceted)
p_age <- ggplot(umap_df_plot%>%
                  dplyr::mutate(age=factor(age,
                                            levels=c("Young","Middleage","Old"))),
                aes(x=UMAP_1,y=UMAP_2,color=neuro_activity)) +
  geom_point(size=0.04,shape=16,alpha=1.0,stroke=0) +
  scale_color_gradient2(low="lightgrey",mid="#1565C0",high="#FF1744",
                         midpoint=median(umap_df_plot$neuro_activity,
                                          na.rm=TRUE),
                         name="Activity") +
  facet_wrap(~age,nrow=1,
             labeller=labeller(age=c(Young="Young",
                                      Middleage="Middle-aged",
                                      Old="Old"))) +
  umap_theme +
  theme(strip.text=element_text(face="bold",size=9),
        strip.background=element_rect(fill="grey92")) +
  labs(title="NEUROeSTIMator -- split by age",x="UMAP 1",y="UMAP 2")
save_both(p_age,file.path(OUT_DIR,"UMAP_ActivityScore_ByAge"),w=18,h=7)

# 7d: Split by sex
p_sex <- ggplot(umap_df_plot%>%
                  dplyr::mutate(sex=factor(sex,levels=c("female","male"))),
                aes(x=UMAP_1,y=UMAP_2,color=neuro_activity)) +
  geom_point(size=0.04,shape=16,alpha=1.0,stroke=0) +
  scale_color_gradient2(low="lightgrey",mid="#1565C0",high="#FF1744",
                         midpoint=median(umap_df_plot$neuro_activity,
                                          na.rm=TRUE),
                         name="Activity") +
  facet_wrap(~sex,nrow=1,
             labeller=labeller(sex=c(female="Female",male="Male"))) +
  umap_theme +
  theme(strip.text=element_text(face="bold",size=9),
        strip.background=element_rect(fill="grey92")) +
  labs(title="NEUROeSTIMator -- split by sex",x="UMAP 1",y="UMAP 2")
save_both(p_sex,file.path(OUT_DIR,"UMAP_ActivityScore_BySex"),w=14,h=7)

cat("Saved UMAP plots\n")

################################################################################
### STEP 8: VIOLIN PLOTS
################################################################################

cat("\n=== Step 8: Violin plots ===\n")
for (sex_k in c("female","male")) {
  neuro_sex <- neuron_df %>% dplyr::filter(sex==sex_k)
  if (nrow(neuro_sex)==0) next
  p_vln <- ggplot(neuro_sex,
                  aes(x=age_condition,y=neuro_activity,fill=age_condition)) +
    geom_violin(trim=FALSE,scale="width",alpha=0.75) +
    geom_boxplot(width=0.1,outlier.size=0.3,alpha=0.5) +
    scale_fill_manual(values=AGE_COND_COLORS,labels=AGE_COND_LABELS,
                       name="Condition",drop=FALSE) +
    scale_x_discrete(labels=AGE_COND_LABELS,drop=FALSE) +
    facet_wrap(~celltype,ncol=3) + coord_cartesian(ylim=Y_LIM) +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
          strip.text=element_text(face="bold",size=9),
          strip.background=element_rect(fill="grey92"),
          plot.title=element_text(face="bold",size=13,hjust=0.5),
          plot.subtitle=element_text(hjust=0.5,size=9,color="grey40"),
          legend.position="bottom") +
    labs(title=paste0("NEUROeSTIMator Activity Score -- ",
                       tools::toTitleCase(sex_k)),
         subtitle="6 age x condition groups | IN = merged interneurons",
         x=NULL,y="Neuronal Activity Score")
  save_both(p_vln,
            file.path(OUT_DIR,paste0(sex_k,"_ViolinPlot_ActivityScore")),
            w=14,h=max(10,ceiling(length(available_neurons)/3)*4+2))
  cat("Saved violin:",sex_k,"\n")
}

################################################################################
### STEP 9: AGE TRAJECTORY
################################################################################

cat("\n=== Step 9: Age trajectory ===\n")
trajectory_df <- neuron_df %>%
  dplyr::group_by(celltype,age,experimental.condition,sex) %>%
  dplyr::summarise(
    mean_activity=mean(neuro_activity,na.rm=TRUE),
    se_activity  =sd(neuro_activity,na.rm=TRUE)/sqrt(dplyr::n()),
    n_cells      =dplyr::n(),
    .groups="drop") %>%
  dplyr::mutate(age=factor(age,levels=c("Young","Middleage","Old")))

write.csv(trajectory_df,
          file.path(OUT_DIR,"ActivityScore_AgeTrajectory_Summary.csv"),
          row.names=FALSE)
TRAJ_LIM <- range(c(trajectory_df$mean_activity - trajectory_df$se_activity,
                    trajectory_df$mean_activity + trajectory_df$se_activity), na.rm=TRUE)   # shared across sexes

for (sex_k in c("female","male")) {
  traj_sex <- trajectory_df %>% dplyr::filter(sex==sex_k)
  if (nrow(traj_sex)==0) next
  p_traj <- ggplot(traj_sex,
                   aes(x=age,y=mean_activity,
                       color=experimental.condition,
                       group=experimental.condition)) +
    geom_line(linewidth=1.2) +
    geom_point(size=3) +
    geom_errorbar(aes(ymin=mean_activity-se_activity,
                       ymax=mean_activity+se_activity),width=0.15) +
    scale_color_manual(values=c(EE="#FF1744",HC="#1565C0"),name="Condition") +
    scale_x_discrete(labels=c(Young="Young",Middleage="Middle-aged",Old="Old")) +
    facet_wrap(~celltype,ncol=3) + coord_cartesian(ylim=TRAJ_LIM) +
    theme_bw(base_size=10) +
    theme(strip.text=element_text(face="bold",size=9),
          strip.background=element_rect(fill="grey92"),
          plot.title=element_text(face="bold",size=13,hjust=0.5),
          plot.subtitle=element_text(hjust=0.5,size=9,color="grey40"),
          legend.position="bottom") +
    labs(title=paste0("Age Trajectory of Neuronal Activity -- ",
                       tools::toTitleCase(sex_k)),
         subtitle="Mean +/- SE | IN = merged interneurons",
         x="Age",y="Mean Neuronal Activity Score")
  save_both(p_traj,
            file.path(OUT_DIR,paste0(sex_k,"_AgeTrajectory_ActivityScore")),
            w=14,h=max(10,ceiling(length(available_neurons)/3)*4+2))
  cat("Saved trajectory:",sex_k,"\n")
}

################################################################################
### STEP 10: EE vs HC WILCOXON
################################################################################

cat("\n=== Step 10: EE vs HC stats ===\n")
stats_results <- neuron_df %>%
  dplyr::group_by(celltype,age,sex) %>%
  dplyr::summarise(
    n_EE   =sum(experimental.condition=="EE"),
    n_HC   =sum(experimental.condition=="HC"),
    mean_EE=mean(neuro_activity[experimental.condition=="EE"],na.rm=TRUE),
    mean_HC=mean(neuro_activity[experimental.condition=="HC"],na.rm=TRUE),
    delta  =mean_EE-mean_HC,
    pval   =tryCatch(
      wilcox.test(
        neuro_activity[experimental.condition=="EE"],
        neuro_activity[experimental.condition=="HC"])$p.value,
      error=function(e) NA_real_),
    .groups="drop") %>%
  dplyr::mutate(
    padj=p.adjust(pval,method="BH"),
    sig =dplyr::case_when(
      is.na(padj) ~"NA",padj<0.001~"***",
      padj<0.01~"**",padj<0.05~"*",TRUE~"ns"),
    age=factor(age,levels=c("Young","Middleage","Old"))) %>%
  dplyr::arrange(padj)

write.csv(stats_results,
          file.path(OUT_DIR,"ActivityScore_EEvsHC_Stats.csv"),
          row.names=FALSE)

cat("Significant EE vs HC (padj<0.05):\n")
sig <- stats_results %>% dplyr::filter(!is.na(padj)&padj<0.05)
if (nrow(sig)>0) {
  print(sig %>% dplyr::select(celltype,age,sex,mean_EE,mean_HC,delta,padj,sig))
} else cat("No significant comparisons\n")

################################################################################
### STEP 11: DELTA ACTIVITY HEATMAP
################################################################################

cat("\n=== Step 11: Delta heatmaps ===\n")
DELTA_LIM <- max(abs(stats_results$delta), na.rm=TRUE)   # shared colour scale across sexes
for (sex_k in c("female","male")) {
  stats_sex <- stats_results %>%
    dplyr::filter(sex==sex_k) %>%
    dplyr::mutate(
      age_label=factor(age,levels=c("Young","Middleage","Old"),
                        labels=c("Young","Middle-aged","Old")),
      sig_label=dplyr::case_when(
        is.na(padj)~"",padj<0.001~"***",
        padj<0.01~"**",padj<0.05~"*",TRUE~""),
      celltype=factor(celltype,
                       levels=rev(intersect(NEURON_CELLTYPES,
                                             unique(celltype)))))
  p_heat <- ggplot(stats_sex,
                   aes(x=age_label,y=celltype,fill=delta)) +
    geom_tile(color="white",linewidth=0.5) +
    geom_text(aes(label=paste0(round(delta,3),"\n",sig_label)),
              size=3,color="black") +
    scale_fill_gradient2(low="#1565C0",mid="white",high="#FF1744",
                          midpoint=0,limits=c(-DELTA_LIM,DELTA_LIM),name="EE-HC\nActivity") +
    theme_minimal(base_size=11) +
    theme(axis.text.x=element_text(size=10,face="bold"),
          axis.text.y=element_text(size=10),
          plot.title=element_text(face="bold",size=13,hjust=0.5),
          plot.subtitle=element_text(hjust=0.5,size=9,color="grey40"),
          panel.grid=element_blank()) +
    labs(title=paste0("EE Effect on Neuronal Activity -- ",
                       tools::toTitleCase(sex_k)),
         subtitle="Delta = mean(EE) - mean(HC) | Red = EE more active",
         caption="* padj<0.05  ** padj<0.01  *** padj<0.001",
         x=NULL,y=NULL)
  save_both(p_heat,
            file.path(OUT_DIR,paste0(sex_k,"_Heatmap_DeltaActivity_EEvsHC")),
            w=9,h=max(5,length(available_neurons)*0.6+3))
  cat("Saved delta heatmap:",sex_k,"\n")
}

################################################################################
### STEP 12: EE ONLY -- AGE COMPARISON
################################################################################

cat("\n=== Step 12: EE age comparison ===\n")
ee_only <- neuron_df %>% dplyr::filter(experimental.condition=="EE")

for (sex_k in c("female","male")) {
  ee_sex <- ee_only %>% dplyr::filter(sex==sex_k)
  if (nrow(ee_sex)==0) next
  p_ee <- ggplot(ee_sex,aes(x=age,y=neuro_activity,fill=age)) +
    geom_violin(trim=FALSE,scale="width",alpha=0.75) +
    geom_boxplot(width=0.1,outlier.size=0.3,alpha=0.5) +
    scale_fill_manual(
      values=c(Young="#1565C0",Middleage="#00897B",Old="#FF1744"),
      name="Age") +
    scale_x_discrete(labels=c(Young="Young",Middleage="Middle-aged",Old="Old")) +
    facet_wrap(~celltype,ncol=3) + coord_cartesian(ylim=Y_LIM) +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
          strip.text=element_text(face="bold",size=9),
          strip.background=element_rect(fill="grey92"),
          plot.title=element_text(face="bold",size=13,hjust=0.5),
          legend.position="bottom") +
    labs(title=paste0("EE Activation Across Ages -- ",
                       tools::toTitleCase(sex_k)),
         x=NULL,y="Neuronal Activity Score")
  save_both(p_ee,
            file.path(OUT_DIR,paste0(sex_k,"_ViolinPlot_EE_AgeComparison")),
            w=14,h=max(10,ceiling(length(available_neurons)/3)*4+2))
}

################################################################################
### STEP 13: HC ONLY -- BASAL AGING
################################################################################

cat("\n=== Step 13: HC basal aging ===\n")
hc_only <- neuron_df %>% dplyr::filter(experimental.condition=="HC")

for (sex_k in c("female","male")) {
  hc_sex <- hc_only %>% dplyr::filter(sex==sex_k)
  if (nrow(hc_sex)==0) next
  p_hc <- ggplot(hc_sex,aes(x=age,y=neuro_activity,fill=age)) +
    geom_violin(trim=FALSE,scale="width",alpha=0.75) +
    geom_boxplot(width=0.1,outlier.size=0.3,alpha=0.5) +
    scale_fill_manual(
      values=c(Young="#90CAF9",Middleage="#80CBC4",Old="#FF8A80"),
      name="Age") +
    scale_x_discrete(labels=c(Young="Young",Middleage="Middle-aged",Old="Old")) +
    facet_wrap(~celltype,ncol=3) + coord_cartesian(ylim=Y_LIM) +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
          strip.text=element_text(face="bold",size=9),
          strip.background=element_rect(fill="grey92"),
          plot.title=element_text(face="bold",size=13,hjust=0.5),
          legend.position="bottom") +
    labs(title=paste0("Basal Neuronal Activity Across Ages -- ",
                       tools::toTitleCase(sex_k)),
         x=NULL,y="Neuronal Activity Score")
  save_both(p_hc,
            file.path(OUT_DIR,paste0(sex_k,"_ViolinPlot_HC_AgeComparison")),
            w=14,h=max(10,ceiling(length(available_neurons)/3)*4+2))
}

cat("\n================================================\n")
cat("  NEUROeSTIMator PIPELINE COMPLETE\n")
cat("  IN = merged PV/SST/VIP/Lamp5 interneurons\n")
cat("================================================\n")
cat("Output:",OUT_DIR,"\n")

invisible(NULL)
}

################################################################################
### STEP 2d  Fig. 2d Venn
################################################################################

SECTIONS[["2d"]] <- function() {
################################################################################
### Fig2D_Venn  —  area-proportional 3-way Venn of the NE programs at each age
###
### Young NE (YEEYHC), Middle NE (MEEMHC), Old NE (OEEOHC): genes with padj < 0.05 in any cell type, per sex.
### Every region's AREA is proportional to its gene count (circle area = set size; centres fitted to the 7 regions),
### and each count sits inside its own region (thin regions get a short leader line). No extra packages.
### Per sex the run prints genes vs area share of each region and the overall fit error (share of area that
### circles cannot place exactly; < 0.05 is a close fit). Sexes are drawn side by side, each on its own scale.
### Reads the DESeq2 tables (step deseq2). Writes PDF + PNG per sex and both sexes together, region counts, gene lists.
################################################################################

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
HAS_PATCH <- requireNamespace("patchwork", quietly=TRUE)

deseq_path <- paste0(DESEQ_DIR, "/")
fig_date   <- format(Sys.Date(), "%Y-%m-%d")
out_dir    <- file.path(FIG_ROOT, "2d_Venn")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
cat("Output:", out_dir, "\n\n")

PADJ_THR  <- 0.05
SET_NAMES <- c("Young", "Old", "Middle")                 # A, B, C (Middle drawn below Young and Old)
SET_COLS  <- c(Young="#4A6FE3", Old="#E8743B", Middle="#7CC46B")

load_sig <- function(sex, comp) {
  f <- file.path(deseq_path, paste0(sex, "_ALL_", comp, ".csv"))
  if (!file.exists(f)) { cat("Not found:", basename(f), "\n"); return(character(0)) }
  df <- read.csv(f, stringsAsFactors=FALSE)
  unique(df$gene[!is.na(df$padj) & df$padj < PADJ_THR])
}

################################################################################
### Area-proportional 3-set Venn (Euler) drawing, no extra packages.
###   - circle area = set size; centres fitted so every one of the 7 regions has (as close as circles allow) an area
###     proportional to its gene count (soft-membership grid + optim, several starts)
###   - each count is placed at the most interior point of its own region; regions too thin to hold a label get a
###     short leader line from outside
###   - prints fitted vs actual share of each region and an overall fit error, so a poor fit is visible
################################################################################

euler3_fit <- function(n, n_grid=200, starts=8, seed=1) {
  # n: named counts of the 7 disjoint regions: A, B, C, AB, AC, BC, ABC (AB = in A and B only, etc.)
  tot <- c(A=n[["A"]] + n[["AB"]] + n[["AC"]] + n[["ABC"]], B=n[["B"]] + n[["AB"]] + n[["BC"]] + n[["ABC"]],
           C=n[["C"]] + n[["AC"]] + n[["BC"]] + n[["ABC"]])
  r <- sqrt(tot / pi); s <- max(r); r <- r / s                       # radii, scaled so the largest = 1
  target <- n[c("A","B","C","AB","AC","BC","ABC")] / s^2             # target region areas in the same units
  g <- seq(-3.2, 3.2, length.out=n_grid); G <- expand.grid(x=g, y=g); cell <- diff(g[1:2])^2
  eps <- 0.03
  regions <- function(cx, cy, soft=TRUE) {
    m <- sapply(1:3, function(i) { d <- sqrt((G$x - cx[i])^2 + (G$y - cy[i])^2)
      if (soft) stats::plogis((r[i] - d) / eps) else as.numeric(d <= r[i]) })
    a <- m[, 1]; b <- m[, 2]; c <- m[, 3]
    cbind(A=a*(1-b)*(1-c), B=(1-a)*b*(1-c), C=(1-a)*(1-b)*c, AB=a*b*(1-c), AC=a*(1-b)*c, BC=(1-a)*b*c, ABC=a*b*c)
  }
  pos <- function(p) list(cx=c(0, p[1], p[2]), cy=c(0, 0, p[3]))     # A at origin, B on the x-axis
  loss <- function(p) { q <- pos(p); a <- colSums(regions(q$cx, q$cy)) * cell; sum((a - target)^2) / sum(target^2) }
  set.seed(seed); best <- NULL
  for (k in seq_len(starts)) {
    p0 <- c(stats::runif(1, 0.2, 1.8), stats::runif(1, -0.5, 1.5), stats::runif(1, -1.6, 1.6))
    o <- stats::optim(p0, loss, method="Nelder-Mead", control=list(maxit=600, reltol=1e-7))
    if (is.null(best) || o$value < best$value) best <- o
  }
  q <- pos(best$par)
  cx <- q$cx - mean(q$cx); cy <- q$cy - mean(q$cy)                    # centre the layout
  if (cy[3] > mean(cy[1:2])) cy <- -cy                                # third set below the first two (canonical look)
  hard <- regions(cx, cy, soft=FALSE)
  fitted <- colSums(hard) * cell
  list(cx=cx, cy=cy, r=r, n=n, tot=tot, target=target, fitted=fitted, G=G, hard=hard,
       error=sum(abs(fitted / sum(fitted) - target / sum(target))) / 2)   # share of total area misplaced (0 = exact)
}

euler3_labels <- function(fit, min_inside=0.10) {
  # most interior grid point of each region = max over its points of the distance to the nearest circle edge
  edge <- sapply(1:3, function(i) abs(fit$r[i] - sqrt((fit$G$x - fit$cx[i])^2 + (fit$G$y - fit$cy[i])^2)))
  dmin <- apply(edge, 1, min)
  ctr <- c(mean(fit$cx), mean(fit$cy))
  do.call(rbind, lapply(colnames(fit$hard), function(rg) {
    nn <- fit$n[[rg]]; if (nn == 0) return(NULL)
    idx <- which(fit$hard[, rg] > 0)
    if (!length(idx)) {   # region the circles could not create: label outside, next to the sets it belongs to
      sets <- strsplit(rg, "")[[1]]; k <- match(sets, c("A","B","C"))
      px <- mean(fit$cx[k]); py <- mean(fit$cy[k]); v <- c(px, py) - ctr; v <- v / max(sqrt(sum(v^2)), 1e-6)
      return(data.frame(region=rg, n=nn, x=px + v[1] * 1.4, y=py + v[2] * 1.4, ax=px, ay=py, leader=TRUE, missing=TRUE))
    }
    j <- idx[which.max(dmin[idx])]; x <- fit$G$x[j]; y <- fit$G$y[j]
    if (dmin[j] >= min_inside) return(data.frame(region=rg, n=nn, x=x, y=y, ax=x, ay=y, leader=FALSE, missing=FALSE))
    v <- c(x, y) - ctr; v <- v / max(sqrt(sum(v^2)), 1e-6)            # thin region: label outside with a leader line
    data.frame(region=rg, n=nn, x=x + v[1] * 0.45, y=y + v[2] * 0.45, ax=x, ay=y, leader=TRUE, missing=FALSE)
  }))
}

euler3_plot <- function(fit, set_names, cols, title=NULL, alpha=0.35, count_size=3.3, name_size=4.2) {
  th <- seq(0, 2 * pi, length.out=361)
  circ <- do.call(rbind, lapply(1:3, function(i) data.frame(set=set_names[i], x=fit$cx[i] + fit$r[i] * cos(th),
                                                            y=fit$cy[i] + fit$r[i] * sin(th))))
  circ$set <- factor(circ$set, levels=set_names)
  lab <- euler3_labels(fit)
  ctr <- c(mean(fit$cx), mean(fit$cy))
  nm <- do.call(rbind, lapply(1:3, function(i) { v <- c(fit$cx[i], fit$cy[i]) - ctr; v <- v / max(sqrt(sum(v^2)), 1e-6)
    data.frame(set=set_names[i], x=fit$cx[i] + v[1] * fit$r[i] * 1.18, y=fit$cy[i] + v[2] * fit$r[i] * 1.18,
               hj=ifelse(v[1] > 0.3, 0, ifelse(v[1] < -0.3, 1, 0.5))) }))
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data=circ, ggplot2::aes(x, y, group=set, fill=set), alpha=alpha) +
    ggplot2::geom_path(data=circ, ggplot2::aes(x, y, group=set, color=set), linewidth=0.9) +
    ggplot2::geom_segment(data=lab[lab$leader, , drop=FALSE], ggplot2::aes(x=x, y=y, xend=ax, yend=ay), linewidth=0.3, color="grey30") +
    ggplot2::geom_label(data=lab, ggplot2::aes(x, y, label=format(n, big.mark=",")), size=count_size,
                        label.size=0, fill=grDevices::adjustcolor("white", 0.65), label.padding=ggplot2::unit(0.08, "lines"),
                        color=ifelse(lab$missing, "grey45", "grey10")) +
    ggplot2::geom_text(data=nm, ggplot2::aes(x, y, label=set, color=set, hjust=hj), size=name_size, fontface="bold") +
    ggplot2::scale_fill_manual(values=stats::setNames(cols, set_names), guide="none") +
    ggplot2::scale_color_manual(values=stats::setNames(cols, set_names), guide="none") +
    ggplot2::coord_equal(clip="off") + ggplot2::theme_void(base_size=11) +
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold", hjust=0.5, size=13), plot.margin=ggplot2::margin(20, 30, 20, 30))
  if (!is.null(title)) p <- p + ggplot2::labs(title=title)
  p
}

euler3_report <- function(fit, set_names) {
  rn <- c(A=set_names[1], B=set_names[2], C=set_names[3])
  nice <- sapply(names(fit$target), function(k) paste(rn[strsplit(k, "")[[1]]], collapse=" & "))
  data.frame(region=nice, genes=as.integer(fit$n[names(fit$target)]),
             pct_genes=round(100 * fit$target / sum(fit$target), 1), pct_area=round(100 * fit$fitted / sum(fit$fitted), 1),
             row.names=NULL)
}

################################################################################
### MAIN
################################################################################

plots <- list(); all_genesets <- list(); reports <- list()
for (current_sex in c("male", "female")) {
  cat("\n###", toupper(current_sex), "###\n")
  young  <- load_sig(current_sex, "YEEYHC"); middle <- load_sig(current_sex, "MEEMHC"); old <- load_sig(current_sex, "OEEOHC")
  if (any(lengths(list(young, middle, old)) == 0)) { cat("Skipping\n"); next }
  A <- young; B <- old; C <- middle
  n <- c(A=length(setdiff(A, union(B, C))), B=length(setdiff(B, union(A, C))), C=length(setdiff(C, union(A, B))),
         AB=length(setdiff(intersect(A, B), C)), AC=length(setdiff(intersect(A, C), B)), BC=length(setdiff(intersect(B, C), A)),
         ABC=length(Reduce(intersect, list(A, B, C))))
  cat(sprintf("  Young=%d  Middle=%d  Old=%d genes\n", length(young), length(middle), length(old)))
  cat("  fitting circles ...\n")
  fit <- euler3_fit(n)
  rep <- euler3_report(fit, SET_NAMES); rep$sex <- current_sex; reports[[current_sex]] <- rep
  print(rep, row.names=FALSE)
  cat(sprintf("  fit error = %.3f (share of area circles cannot place exactly; < 0.05 = close)\n", fit$error))
  sx <- tools::toTitleCase(current_sex)
  p <- euler3_plot(fit, SET_NAMES, unname(SET_COLS[SET_NAMES]), title=paste0(sx, "s"))
  plots[[current_sex]] <- p
  for (ext in c("pdf", "png")) ggsave(file.path(out_dir, paste0(fig_date, "_Fig2D_Venn_", sx, ".", ext)), p,
                                       width=6, height=5.6, dpi=600, bg="white")
  cat("  Saved:", paste0(fig_date, "_Fig2D_Venn_", sx, ".pdf/.png"), "\n")

  all_genes <- unique(c(young, middle, old))
  all_genesets[[current_sex]] <- data.frame(gene=all_genes, sex=current_sex, in_young=all_genes %in% young,
                                            in_middle=all_genes %in% middle, in_old=all_genes %in% old) %>%
    dplyr::mutate(set_label=dplyr::case_when(
      in_young & in_middle & in_old ~ "All three (core)", in_young & in_middle & !in_old ~ "Young + Middle",
      in_young & !in_middle & in_old ~ "Young + Old", !in_young & in_middle & in_old ~ "Middle + Old",
      in_young & !in_middle & !in_old ~ "Young only (lost)", !in_young & in_middle & !in_old ~ "Middle only",
      !in_young & !in_middle & in_old ~ "Old only (age-emergent)", TRUE ~ "other"))
}

if (length(plots) == 2 && HAS_PATCH) {
  pb <- patchwork::wrap_plots(plots[["male"]], plots[["female"]], nrow=1)
  for (ext in c("pdf", "png")) ggsave(file.path(out_dir, paste0(fig_date, "_Fig2D_Venn_BothSexes.", ext)), pb,
                                       width=11.5, height=5.6, dpi=600, bg="white")
  cat("\nSaved:", paste0(fig_date, "_Fig2D_Venn_BothSexes.pdf/.png"), "(each sex on its own scale)\n")
}
if (length(reports)) write.csv(dplyr::bind_rows(reports), file.path(out_dir, paste0(fig_date, "_Fig2D_Venn_RegionCounts.csv")), row.names=FALSE)
if (length(all_genesets)) write.csv(dplyr::bind_rows(all_genesets), file.path(out_dir, paste0(fig_date, "_Fig2D_GeneSets_All.csv")), row.names=FALSE)
cat("\n=== DONE ===\nOutputs:", out_dir, "\n")

invisible(NULL)
}

################################################################################
### STEP 2f  Fig. 2f pathways of lost and maintained genes
################################################################################

SECTIONS[["2f"]] <- function() {
################################################################################
### Fig2F_Pathway_LostMaintained
###
### GO/KEGG enrichment for Lost and Maintained gene sets from YEEYHC anchor.
###
### GENE SETS:
###   Maintained Up:   sig up in YEEYHC + sig up in comparison
###   Maintained Down: sig down in YEEYHC + sig down in comparison
###   Lost Up:         sig up in YEEYHC + NOT sig in comparison
###   Lost Down:       sig down in YEEYHC + NOT sig in comparison
###
### VISUALIZATION: Gene-level Sankey + Dot plot (SRplot style)
###   Left:   individual gene symbols -> flow -> pathway/GO term
###   Right:  Dot (size=Count, color=-log10 p.adjust)
###   Style:  bioinformatics.com.cn SRplot gene-sankey-dot
###
### Comparisons: MEEMHC (Middle EE) and OEEOHC (Old EE)
### Output: Figures/Figure2/2f_Pathway_LostMaintained/
################################################################################

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
suppressPackageStartupMessages({
  library(conflicted)
  library(ggalluvial)
  library(patchwork)
  conflicts_prefer(dplyr::select,           .quiet=TRUE)
  conflicts_prefer(dplyr::filter,           .quiet=TRUE)
  conflicts_prefer(base::intersect,         .quiet=TRUE)
  conflicts_prefer(base::setdiff,           .quiet=TRUE)
  conflicts_prefer(base::union,             .quiet=TRUE)
  conflicts_prefer(clusterProfiler::select, .quiet=TRUE)
  conflicts_prefer(AnnotationDbi::select,   .quiet=TRUE)
  conflicts_prefer(base::as.data.frame,     .quiet=TRUE)
})

################################################################################
### PATHS
################################################################################

deseq_dir <- paste0(DESEQ_DIR, "/")
fig_date  <- format(Sys.Date(),"%Y-%m-%d")
out_dir   <- paste0(file.path(FIG_ROOT, "2f_Pathway_LostMaintained"), "/")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
cat("Output:",out_dir,"\n\n")

PADJ_THR  <- 0.05
Q_THR     <- 0.2
MIN_GENES <- 5
TOP_N     <- 10   # top pathway terms shown in sankey+dot
COMPS     <- c("MEEMHC","OEEOHC")

should_skip <- function(f) {
  if (file.exists(f)) { cat("  Skipping:",basename(f),"\n"); return(TRUE) }
  FALSE
}

################################################################################
### PALETTE HELPERS
################################################################################

# Distinct palette for up to ~15 pathways; falls back to hue_pal for more
pathway_palette <- function(n) {
  base_cols <- c(
    "#E07B54","#5B8DB8","#6BAF6B","#8F7BB5","#C0B45A",
    "#7BBFBF","#C97BAB","#A07850","#7090B8","#90B870",
    "#B07070","#70A0A0","#A07090","#90A070","#A09060")
  if (n <= length(base_cols)) return(base_cols[seq_len(n)])
  scales::hue_pal()(n)
}

################################################################################
### GENE-LEVEL SANKEY + DOT PLOT
### Matches the bioinformatics.com.cn SRplot gene-sankey-dot style:
###   Left axis  : individual gene symbols (colored by their pathway)
###   Flows      : gene -> pathway, colored by pathway
###   Right panel: dot plot (x=GeneRatio, y=pathway, size=Count,
###                           color=-log10 p.adjust)
################################################################################

sankey_dot_plot <- function(enrich_res, title, out_file,
                             top_n    = TOP_N,
                             width    = 14,
                             height   = NULL,   # auto-sized from gene count
                             low_col  = "#3182BD",
                             high_col = "#E6550D") {

  if (is.null(enrich_res) || nrow(enrich_res) == 0) return(invisible(NULL))
  if (should_skip(out_file)) return(invisible(NULL))

  tryCatch({
    # ── 1. Prepare pathway-level data ───────────────────────────────────────
    df <- as.data.frame(enrich_res) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::mutate(
        GeneRatio_num = sapply(GeneRatio, function(x) {
          pts <- strsplit(x, "/")[[1]]
          as.numeric(pts[1]) / as.numeric(pts[2])
        }),
        neg_log10_padj = -log10(p.adjust),
        Desc_short = ifelse(
          nchar(Description) > 45,
          paste0(substr(Description, 1, 42), "..."),
          Description))

    # Pathway factor: top term (most significant) at top of dot plot
    pathway_levels <- rev(df$Desc_short)   # bottom=worst, top=best
    df$Desc_short  <- factor(df$Desc_short, levels = pathway_levels)

    # ── 2. Build gene-to-pathway edge table ─────────────────────────────────
    # geneID holds "/" delimited gene symbols (readable=TRUE in enrichGO/setReadable)
    edge_df <- df %>%
      dplyr::select(Desc_short, geneID, neg_log10_padj) %>%
      dplyr::mutate(genes = strsplit(as.character(geneID), "/")) %>%
      tidyr::unnest(genes) %>%
      dplyr::rename(Gene = genes, Pathway = Desc_short) %>%
      dplyr::mutate(Pathway = factor(Pathway, levels = levels(df$Desc_short)))

    # ── 3. Assign per-pathway colors ────────────────────────────────────────
    n_paths     <- nlevels(df$Desc_short)
    path_colors <- setNames(pathway_palette(n_paths), levels(df$Desc_short))

    # ── 4. Order genes: by rank of their most significant pathway ────────────
    gene_rank <- edge_df %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(top_rank = min(as.integer(Pathway)), .groups = "drop") %>%
      dplyr::arrange(top_rank, Gene)

    gene_levels  <- unique(gene_rank$Gene)
    edge_df$Gene <- factor(edge_df$Gene, levels = gene_levels)
    n_genes      <- length(gene_levels)

    # ── 5. LEFT PANEL: Gene-level Sankey ────────────────────────────────────
    # Using fill = Pathway in BOTH geom_alluvium and geom_stratum:
    #   * flows are colored by their destination pathway
    #   * gene strata (axis1) are split-colored by pathway proportions
    #     → produces the colored mini-bars on gene labels (matches reference)
    #   * pathway strata (axis2) are solid in their pathway color
    p_sankey <- ggplot(
      edge_df,
      aes(axis1 = Gene, axis2 = Pathway, y = 1)) +
      ggalluvial::geom_alluvium(
        aes(fill = Pathway),
        width      = 0.25,
        alpha      = 0.5,
        curve_type = "sigmoid",
        knot.pos   = 0.4) +
      ggalluvial::geom_stratum(
        aes(fill = Pathway),        # split-fills gene strata by pathway mix
        width     = 0.25,
        color     = "grey70",
        linewidth = 0.25) +
      # Gene labels on left stratum
      geom_text(
        stat     = "stratum",
        aes(label = after_stat(stratum)),
        size     = 2.6,
        hjust    = 1,
        nudge_x  = -0.15,
        color    = "grey15",
        fontface = "italic") +
      scale_fill_manual(values = path_colors, guide = "none") +
      scale_x_discrete(
        limits = c("Gene", "Pathway"),
        expand = c(0.35, 0.05)) +
      theme_void(base_size = 9) +
      theme(
        axis.text.x = element_text(size = 9, face = "bold", color = "grey30",
                                   margin = margin(t = 4)),
        plot.margin = margin(10, 0, 10, 65))

    # ── 6. RIGHT PANEL: Dot plot ─────────────────────────────────────────────
    # Thin colored strip on left edge mirrors pathway colors from Sankey
    dot_xmin <- min(df$GeneRatio_num) * 0.75
    strip_w  <- (max(df$GeneRatio_num) - min(df$GeneRatio_num)) * 0.04

    p_dot <- ggplot(df,
                    aes(x = GeneRatio_num, y = Desc_short)) +
      # Color strip
      geom_tile(
        aes(x = dot_xmin - strip_w / 2, fill = Desc_short),
        width  = strip_w,
        height = 0.75,
        color  = NA) +
      # Dots
      geom_point(
        aes(color = neg_log10_padj, size = Count),
        alpha = 0.9) +
      scale_color_gradient(
        low  = low_col, high = high_col,
        name = bquote(-log[10](italic(p)[adj]))) +
      scale_fill_manual(values = path_colors, guide = "none") +
      scale_size_continuous(range = c(2, 9), name = "Count") +
      scale_x_continuous(
        expand = expansion(mult = c(0.15, 0.18)),
        labels = scales::number_format(accuracy = 0.01)) +
      theme_classic(base_size = 9) +
      theme(
        axis.text.y     = element_text(size = 8, color = "grey15"),
        axis.ticks.y    = element_blank(),
        axis.line.y     = element_blank(),
        axis.text.x     = element_text(size = 7.5),
        axis.title.x    = element_text(size = 8.5),
        plot.margin     = margin(10, 12, 10, 4),
        legend.position = "right",
        legend.text     = element_text(size = 7),
        legend.title    = element_text(size = 7.5)) +
      labs(x = "Gene Ratio", y = NULL)

    # ── 7. Combine ───────────────────────────────────────────────────────────
    auto_h <- if (is.null(height)) max(5, n_genes * 0.22 + 2) else height

    combined <- p_sankey + p_dot +
      patchwork::plot_layout(widths = c(2.2, 1)) +
      patchwork::plot_annotation(
        title    = title,
        subtitle = paste0("Top ", nrow(df), " terms | ordered by p.adjust"),
        theme    = theme(
          plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle   = element_text(size = 8.5, color = "grey40", hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA)))

    ggsave(out_file, combined,
           width     = width,
           height    = auto_h,
           device    = "pdf",
           limitsize = FALSE)
    cat("  Saved:", basename(out_file), "\n")

  }, error = function(e) cat("  sankey_dot error:", e$message, "\n"))
}

################################################################################
### COMPARECLUSTER SANKEY+DOT (gene-level, faceted by cluster)
################################################################################

compare_sankey_dot <- function(gene_lists, universe, ont = "BP",
                                title, out_file, top_n = 8,
                                width = 16, height = 10) {
  gene_lists <- gene_lists[lengths(gene_lists) >= MIN_GENES]
  if (length(gene_lists) < 2) { cat("  Too few sets\n"); return(invisible(NULL)) }
  if (should_skip(out_file)) return(invisible(NULL))

  tryCatch({
    res <- compareCluster(
      geneClusters  = gene_lists, fun = "enrichGO",
      OrgDb         = org.Mm.eg.db, keyType = "SYMBOL", ont = ont,
      pAdjustMethod = "BH", pvalueCutoff = 0.05,
      universe      = universe, readable = TRUE)

    if (is.null(res) || nrow(res) == 0) {
      cat("  No significant terms\n"); return(invisible(NULL))
    }

    write.csv(as.data.frame(res),
              gsub(".pdf", ".csv", out_file, fixed = TRUE),
              row.names = FALSE)

    df <- as.data.frame(res) %>%
      dplyr::group_by(Cluster) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        GeneRatio_num = sapply(GeneRatio, function(x) {
          pts <- strsplit(x, "/")[[1]]
          as.numeric(pts[1]) / as.numeric(pts[2])
        }),
        neg_log10_padj = -log10(p.adjust),
        Desc_short = ifelse(nchar(Description) > 45,
                             paste0(substr(Description, 1, 42), "..."),
                             Description))

    top_terms <- df %>%
      dplyr::group_by(Description) %>%
      dplyr::summarise(min_p = min(p.adjust), .groups = "drop") %>%
      dplyr::arrange(min_p) %>%
      dplyr::pull(Description)

    df <- df %>%
      dplyr::mutate(Description = factor(Description, levels = rev(top_terms)))

    n_clusters <- length(unique(df$Cluster))

    p_dot <- ggplot(df,
                    aes(x = Cluster, y = Description,
                        color = neg_log10_padj, size = GeneRatio_num)) +
      geom_point(alpha = 0.9) +
      scale_color_gradient(
        low  = "#3182BD", high = "#E6550D",
        name = bquote(-log[10](italic(p)[adj]))) +
      scale_size_continuous(range = c(1.5, 7), name = "Gene ratio") +
      theme_classic(base_size = 8) +
      theme(
        axis.text.x   = element_text(angle = 35, hjust = 1, size = 7.5),
        axis.text.y   = element_text(size = 7.5),
        axis.title    = element_blank(),
        legend.position = "right",
        legend.text   = element_text(size = 7),
        legend.title  = element_text(size = 7.5),
        plot.margin   = margin(5, 10, 5, 5)) +
      patchwork::plot_annotation(
        title    = title,
        subtitle = paste0("Top ", top_n, " terms per gene set | GO ", ont),
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 8.5, color = "grey40", hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA)))

    ggsave(out_file, p_dot,
           width  = max(width, n_clusters * 1.0 + 3.5),
           height = height,
           device = "pdf")
    cat("  Saved compareCluster Dot:", basename(out_file), "\n")
    return(res)

  }, error = function(e) {
    cat("  compare_sankey_dot error:", e$message, "\n"); NULL
  })
}

################################################################################
### ENRICHMENT HELPERS
################################################################################

run_go <- function(genes, universe, ont = "BP") {
  genes <- intersect(genes, universe)
  if (length(genes) < MIN_GENES) return(NULL)
  tryCatch({
    res <- enrichGO(
      gene          = genes,  universe = universe,
      OrgDb         = org.Mm.eg.db, keyType = "SYMBOL",
      ont           = ont,    pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,   qvalueCutoff  = Q_THR,
      readable      = TRUE)
    if (is.null(res) || nrow(res) == 0) return(NULL)
    clusterProfiler::simplify(res, cutoff = 0.7, by = "p.adjust")
  }, error = function(e) NULL)
}

sym_to_entrez <- function(genes) {
  res <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = genes,
                                column   = "ENTREZID",
                                keytype  = "SYMBOL",
                                multiVals = "first")
  as.character(res[!is.na(res)])
}

run_kegg <- function(genes, universe) {
  genes <- intersect(genes, universe)
  if (length(genes) < MIN_GENES) return(NULL)
  eg <- sym_to_entrez(genes)
  eu <- sym_to_entrez(universe)
  if (length(eg) < MIN_GENES) return(NULL)
  tryCatch({
    res <- enrichKEGG(
      gene          = eg,  universe     = eu,
      organism      = "mmu", pAdjustMethod = "BH",
      pvalueCutoff  = 0.05, qvalueCutoff  = Q_THR)
    if (is.null(res) || nrow(res) == 0) return(NULL)
    # setReadable so geneID column contains symbols, not Entrez IDs
    setReadable(res, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  }, error = function(e) NULL)
}

################################################################################
### GENE SET BUILDER
################################################################################

build_genesets <- function(yee_df, comp_df, ct) {
  y <- yee_df  %>% dplyr::filter(celltype == ct)
  d <- comp_df %>% dplyr::filter(celltype == ct)

  y_up  <- y$gene[!is.na(y$padj) & y$padj < PADJ_THR & y$log2FoldChange > 0]
  y_dn  <- y$gene[!is.na(y$padj) & y$padj < PADJ_THR & y$log2FoldChange < 0]
  d_up  <- d$gene[!is.na(d$padj) & d$padj < PADJ_THR & d$log2FoldChange > 0]
  d_dn  <- d$gene[!is.na(d$padj) & d$padj < PADJ_THR & d$log2FoldChange < 0]
  d_sig <- union(d_up, d_dn)

  list(
    maintained_up   = intersect(y_up, d_up),
    maintained_down = intersect(y_dn, d_dn),
    lost_up         = setdiff(y_up, d_sig),
    lost_down       = setdiff(y_dn, d_sig))
}

################################################################################
### MAIN LOOP
################################################################################

for (current_sex in c("male", "female")) {
  cat("\n########################################\n")
  cat("### SEX:", toupper(current_sex), "\n")
  cat("########################################\n\n")

  sex_out <- file.path(out_dir, current_sex)
  dir.create(sex_out, recursive = TRUE, showWarnings = FALSE)

  yee_file <- file.path(deseq_dir,
                         paste0(current_sex, "_ALL_YEEYHC.csv"))
  if (!file.exists(yee_file)) { cat("  YEEYHC not found\n"); next }

  yee_all   <- read.csv(yee_file, stringsAsFactors = FALSE)
  universe  <- unique(yee_all$gene[!is.na(yee_all$gene)])
  celltypes <- sort(unique(yee_all$celltype))

  cat(sprintf("Universe: %d genes | Cell types: %d\n",
              length(universe), length(celltypes)))

  for (comp in COMPS) {
    cat(sprintf("\n--- Comparison: %s ---\n", comp))

    comp_file <- file.path(deseq_dir,
                            paste0(current_sex, "_ALL_", comp, ".csv"))
    if (!file.exists(comp_file)) {
      cat("  Not found:", basename(comp_file), "\n"); next
    }
    comp_all <- read.csv(comp_file, stringsAsFactors = FALSE)
    comp_out <- file.path(sex_out, comp)
    dir.create(comp_out, recursive = TRUE, showWarnings = FALSE)

    pool_up   <- list()
    pool_dn   <- list()
    pool_comb <- list()

    ##########################################################
    ### PER CELL TYPE
    ##########################################################
    cat("\n  Per cell type:\n")

    for (ct in celltypes) {
      gs      <- build_genesets(yee_all, comp_all, ct)
      ct_slug <- gsub("[^A-Za-z0-9]+", "_", ct)
      ct_out  <- file.path(comp_out, ct_slug)
      dir.create(ct_out, recursive = TRUE, showWarnings = FALSE)

      cat(sprintf("  %-22s  mU=%d mD=%d lU=%d lD=%d\n",
                  ct,
                  length(gs$maintained_up),   length(gs$maintained_down),
                  length(gs$lost_up),          length(gs$lost_down)))

      # Accumulate for pooled compareCluster
      if (length(gs$maintained_up)   >= MIN_GENES)
        pool_up[[paste0(ct, "\nMaintained Up")]]   <- gs$maintained_up
      if (length(gs$lost_up)         >= MIN_GENES)
        pool_up[[paste0(ct, "\nLost Up")]]         <- gs$lost_up
      if (length(gs$maintained_down) >= MIN_GENES)
        pool_dn[[paste0(ct, "\nMaintained Down")]] <- gs$maintained_down
      if (length(gs$lost_down)       >= MIN_GENES)
        pool_dn[[paste0(ct, "\nLost Down")]]       <- gs$lost_down

      comb_maintained <- unique(c(gs$maintained_up,  gs$maintained_down))
      comb_lost       <- unique(c(gs$lost_up,        gs$lost_down))
      if (length(comb_maintained) >= MIN_GENES)
        pool_comb[[paste0(ct, "\nMaintained")]] <- comb_maintained
      if (length(comb_lost)       >= MIN_GENES)
        pool_comb[[paste0(ct, "\nLost")]]       <- comb_lost

      # Per cell type sankey+dot plots
      sets_to_run <- list(
        list(genes = gs$maintained_up,   slug = "Maintained_Up",
             label = "Maintained Up"),
        list(genes = gs$maintained_down, slug = "Maintained_Down",
             label = "Maintained Down"),
        list(genes = gs$lost_up,         slug = "Lost_Up",
             label = "Lost Up"),
        list(genes = gs$lost_down,       slug = "Lost_Down",
             label = "Lost Down"),
        list(genes = comb_maintained,    slug = "Maintained_Combined",
             label = "Maintained (Up+Down)"),
        list(genes = comb_lost,          slug = "Lost_Combined",
             label = "Lost (Up+Down)"))

      for (s in sets_to_run) {
        if (length(s$genes) < MIN_GENES) next

        # GO BP sankey+dot
        out_go <- file.path(ct_out,
          paste0(fig_date, "_", ct_slug, "_", s$slug, "_GO_BP_SankeyDot.pdf"))
        if (!should_skip(out_go)) {
          r <- run_go(s$genes, universe, "BP")
          if (!is.null(r)) {
            write.csv(as.data.frame(r),
                      gsub(".pdf", ".csv", out_go, fixed = TRUE),
                      row.names = FALSE)
            sankey_dot_plot(r,
              paste0(ct, " | ", s$label, " | GO BP | ", comp),
              out_go, top_n = TOP_N)
          }
        }

        # KEGG sankey+dot
        out_kk <- file.path(ct_out,
          paste0(fig_date, "_", ct_slug, "_", s$slug, "_KEGG_SankeyDot.pdf"))
        if (!should_skip(out_kk)) {
          rk <- run_kegg(s$genes, universe)
          if (!is.null(rk)) {
            write.csv(as.data.frame(rk),
                      gsub(".pdf", ".csv", out_kk, fixed = TRUE),
                      row.names = FALSE)
            sankey_dot_plot(rk,
              paste0(ct, " | ", s$label, " | KEGG | ", comp),
              out_kk, top_n = TOP_N)
          }
        }
      }
    }

    ##########################################################
    ### POOLED COMPARECLUSTER DOT PLOT
    ##########################################################
    cat(sprintf("\n  Pooled compareCluster Dot (%s):\n", comp))

    compare_sankey_dot(
      pool_up, universe, ont = "BP",
      title    = paste0("GO BP -- Up genes | ",
                         tools::toTitleCase(current_sex), " | ", comp,
                         "\nMaintained Up vs Lost Up per cell type"),
      out_file = file.path(sex_out,
        paste0(fig_date, "_", current_sex, "_", comp,
               "_UP_compareCluster_Dot_GO_BP.pdf")),
      top_n = 8, width = 16, height = 10)

    compare_sankey_dot(
      pool_dn, universe, ont = "BP",
      title    = paste0("GO BP -- Down genes | ",
                         tools::toTitleCase(current_sex), " | ", comp,
                         "\nMaintained Down vs Lost Down per cell type"),
      out_file = file.path(sex_out,
        paste0(fig_date, "_", current_sex, "_", comp,
               "_DOWN_compareCluster_Dot_GO_BP.pdf")),
      top_n = 8, width = 16, height = 10)

    compare_sankey_dot(
      pool_comb, universe, ont = "BP",
      title    = paste0("GO BP -- Combined | ",
                         tools::toTitleCase(current_sex), " | ", comp,
                         "\nMaintained vs Lost (Up+Down) per cell type"),
      out_file = file.path(sex_out,
        paste0(fig_date, "_", current_sex, "_", comp,
               "_Combined_compareCluster_Dot_GO_BP.pdf")),
      top_n = 8, width = 16, height = 10)

    ##########################################################
    ### SUPER-POOLED (all cell types merged)
    ##########################################################
    cat(sprintf("\n  Super-pooled (%s):\n", comp))

    super <- list(
      "Maintained Up"   = unique(unlist(lapply(celltypes, function(ct)
        build_genesets(yee_all, comp_all, ct)$maintained_up))),
      "Maintained Down" = unique(unlist(lapply(celltypes, function(ct)
        build_genesets(yee_all, comp_all, ct)$maintained_down))),
      "Lost Up"         = unique(unlist(lapply(celltypes, function(ct)
        build_genesets(yee_all, comp_all, ct)$lost_up))),
      "Lost Down"       = unique(unlist(lapply(celltypes, function(ct)
        build_genesets(yee_all, comp_all, ct)$lost_down))))

    for (nm in names(super)) {
      genes   <- super[[nm]]
      nm_slug <- gsub("[^A-Za-z0-9]+", "_", nm)
      if (length(genes) < MIN_GENES) next
      cat(sprintf("    %s: %d genes\n", nm, length(genes)))

      out_sp_go <- file.path(sex_out,
        paste0(fig_date, "_", current_sex, "_", comp,
               "_SuperPooled_", nm_slug, "_GO_BP_SankeyDot.pdf"))
      if (!should_skip(out_sp_go)) {
        r <- run_go(genes, universe, "BP")
        if (!is.null(r)) {
          write.csv(as.data.frame(r),
                    gsub(".pdf", ".csv", out_sp_go, fixed = TRUE),
                    row.names = FALSE)
          sankey_dot_plot(r,
            paste0("Super-pooled: ", nm, " | ", current_sex, " | ", comp),
            out_sp_go, top_n = TOP_N, width = 13)
        }
      }

      out_sp_kk <- file.path(sex_out,
        paste0(fig_date, "_", current_sex, "_", comp,
               "_SuperPooled_", nm_slug, "_KEGG_SankeyDot.pdf"))
      if (!should_skip(out_sp_kk)) {
        rk <- run_kegg(genes, universe)
        if (!is.null(rk)) {
          write.csv(as.data.frame(rk),
                    gsub(".pdf", ".csv", out_sp_kk, fixed = TRUE),
                    row.names = FALSE)
          sankey_dot_plot(rk,
            paste0("Super-pooled: ", nm, " | ", current_sex, " | ", comp, " | KEGG"),
            out_sp_kk, top_n = TOP_N, width = 13)
        }
      }
    }

    # Super-pooled compareCluster dot (all 4 sets together)
    super_filtered <- super[lengths(super) >= MIN_GENES]
    if (length(super_filtered) >= 2) {
      compare_sankey_dot(
        super_filtered, universe, ont = "BP",
        title    = paste0("GO BP -- Super-pooled | ",
                           tools::toTitleCase(current_sex), " | ", comp),
        out_file = file.path(sex_out,
          paste0(fig_date, "_", current_sex, "_", comp,
                 "_SuperPooled_compareCluster_Dot_GO_BP.pdf")),
        top_n = 5, width = 14, height = 11)
    }
  }
}

cat("\n=== DONE ===\nOutput:", out_dir, "\n\n")
cat("GENE SET DEFINITIONS:\n")
cat("  Maintained Up:   sig up in YEEYHC + sig up in comparison\n")
cat("  Maintained Down: sig down in YEEYHC + sig down in comparison\n")
cat("  Lost Up:         sig up in YEEYHC + NOT sig in comparison\n")
cat("  Lost Down:       sig down in YEEYHC + NOT sig in comparison\n\n")
cat("VISUALIZATION: Gene-level Sankey+Dot (SRplot style)\n")
cat("  Left:  Gene symbols -> flow (colored by pathway) -> Pathway terms\n")
cat("  Right: Dot (size=Count, color=-log10 p.adjust, strip=pathway color)\n\n")
cat("REQUIRES: ggalluvial, patchwork, scales, clusterProfiler, enrichplot\n")

invisible(NULL)
}

################################################################################
### STEP sv_all  Strip volcanoes, every comparison
################################################################################

SECTIONS[["sv_all"]] <- function() {

################################################################################
### StripVolcano_Violin_AllComparisons
### No log2FC ceiling (axis fits the data, or set Y_FIXED).
###
### Violin-shaped strip volcano for every EE comparison, both sexes.
### Each comparison is self-contained: sig (padj < 0.05) coloured,
### NS grey.  No cross-comparison grey-out.
###
### Requires: ggbeeswarm (geom_quasirandom); falls back to manual spread.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(conflicted)
  conflicts_prefer(dplyr::select,   .quiet = TRUE)
  conflicts_prefer(dplyr::filter,   .quiet = TRUE)
  conflicts_prefer(base::intersect, .quiet = TRUE)
  conflicts_prefer(base::setdiff,   .quiet = TRUE)
  conflicts_prefer(base::union,     .quiet = TRUE)
})

USE_BEESWARM <- tryCatch({ library(ggbeeswarm); TRUE },
                          error = function(e) {
                            message("ggbeeswarm unavailable — manual spread")
                            FALSE
                          })

################################################################################
### PATHS & CONSTANTS
################################################################################

deseq_dir <- DESEQ_DIR
fig_date  <- format(Sys.Date(), "%Y-%m-%d")
out_dir   <- file.path(FIG_ROOT, "StripVolcano_AllComparisons")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output dir:", out_dir, "\n\n")

PADJ_THR <- 0.05
Y_FIXED  <- NULL       # NULL = each plot's axis fits its own genes (see SV_Y_FIT at the top);
                       #   a number (e.g. 12) = same symmetric axis on every plot
PT_SIZE_NS  <- 0.20    # point size for NS
PT_SIZE_SIG <- 0.55    # point size for significant
STRIP_W  <- 0.35       # half-width of quasirandom spread

# Comparisons to loop over (csv prefix after "{sex}_ALL_")
COMPARISONS <- c("YEEYHC", "MEEMHC", "OEEOHC",
                  "OEEYEE", "MEEYEE", "OEEMEE",
                  "OHCYHC", "MHCYHC")

CT_ORDER <- c(
  "DG ExNs", "CA3 ExNs", "CA1 ExNs", "CA2 ExNs", "Sub ExNs",
  "IN",
  "Astrocytes", "Oligodendrocytes", "OPC", "Microglia",
  "Endothelial", "Pericyte", "APC", "Stromal",
  "Cajal-Retzius Cell", "Macrophage", "Epithelial", "T Cell"
)

col_vals <- c(
  "Up"   = "#D6604D",
  "Down" = "#4393C3",
  "NS"   = "grey82"
)

################################################################################
### HELPERS
################################################################################

# save_both(): saves PDF + PNG at 600 dpi
save_both <- function(p, stem, w, h) {
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h,
         device = "pdf", limitsize = FALSE)
  ggsave(paste0(stem, ".png"), p, width = w, height = h,
         device = "png", dpi = 600, limitsize = FALSE)
  cat("  Saved:", basename(stem), "[pdf + png]\n")
}

# manual_xoff(): quasirandom x-offsets without ggbeeswarm
# Scales spread by sqrt(bin density) so dense bins spread more — violin-like.
manual_xoff <- function(lfc_vec, max_width = STRIP_W, n_bins = 80) {
  bins   <- cut(lfc_vec, breaks = n_bins, labels = FALSE,
                include.lowest = TRUE)
  bins[is.na(bins)] <- 1L
  bsize  <- ave(rep(1, length(bins)), bins, FUN = length)
  idx    <- ave(seq_along(bins), bins, FUN = seq_along)
  # centre-offset scaled by sqrt(bsize) relative to global max
  scale  <- sqrt(bsize) / sqrt(max(bsize))
  xoff   <- ((idx - 1) - (bsize - 1) / 2) / pmax(1, bsize - 1) *
              max_width * scale
  pmax(pmin(xoff, max_width), -max_width)
}

################################################################################
### BUILD PLOT FUNCTION
################################################################################

make_violin_strip <- function(df, comp, sex, ct_levels) {

  y_max <- if (!is.null(Y_FIXED)) Y_FIXED else sv_ymax(df$lfc_plot, df$color_cat != "NS")
  y_brk <- sv_break_at(df$lfc_plot, df$color_cat != "NS", y_max)

  n_cts <- length(ct_levels)
  df$celltype <- factor(df$celltype, levels = ct_levels)

  # ── Violin background ──────────────────────────────────────────────────────
  # Draw violin for ALL points to show the density shape; fill by sig status
  # but keep it light so points are readable.
  violin_layer <- geom_violin(
    data        = df,
    aes(x = celltype, y = lfc_plot),
    fill        = "grey92",
    color       = "grey70",
    linewidth   = 0.25,
    scale       = "width",     # all violins same max width — comparable shape
    trim        = TRUE,
    bw          = "nrd0",
    alpha       = 0.6
  )

  # ── Point layers ───────────────────────────────────────────────────────────
  if (USE_BEESWARM) {

    p <- ggplot(mapping = aes(x = celltype, y = lfc_plot)) +
      violin_layer +
      geom_quasirandom(
        data    = df %>% dplyr::filter(color_cat == "NS"),
        color   = col_vals["NS"], size = PT_SIZE_NS, alpha = 0.30,
        width   = STRIP_W, method = "quasirandom") +
      geom_quasirandom(
        data    = df %>% dplyr::filter(color_cat != "NS"),
        aes(color = color_cat),
        size    = PT_SIZE_SIG, alpha = 0.85,
        width   = STRIP_W, method = "quasirandom")

  } else {

    # Add manual offsets per cell type
    df <- df %>%
      dplyr::group_by(celltype) %>%
      dplyr::mutate(.xoff = manual_xoff(lfc_plot)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(.xnum = as.integer(celltype) + .xoff)

    p <- ggplot() +
      # Re-draw violin using standard aes (celltype still factor)
      geom_violin(
        data      = df,
        aes(x = celltype, y = lfc_plot),
        fill      = "grey92", color = "grey70", linewidth = 0.25,
        scale     = "width", trim = TRUE, alpha = 0.6
      ) +
      geom_point(
        data  = df %>% dplyr::filter(color_cat == "NS"),
        aes(x = .xnum, y = lfc_plot),
        color = col_vals["NS"], size = PT_SIZE_NS, alpha = 0.30
      ) +
      geom_point(
        data  = df %>% dplyr::filter(color_cat != "NS"),
        aes(x = .xnum, y = lfc_plot, color = color_cat),
        size = PT_SIZE_SIG, alpha = 0.85
      ) +
      scale_x_continuous(
        breaks = seq_along(ct_levels),
        labels = ct_levels,
        expand = expansion(add = 0.7)
      )
  }

  # ── Shared aesthetics ──────────────────────────────────────────────────────
  p +
    geom_hline(yintercept = 0, linewidth = 0.45, color = "grey10") +
    scale_color_manual(
      values = col_vals,
      breaks = c("Up", "Down"),
      labels = c("up-regulated", "down-regulated"),
      name   = NULL
    ) +
    sv_break_marks(y_brk) +
    scale_y_continuous(
      limits = c(-y_max, y_max),
      oob    = scales::squish,          # genes beyond the axis sit on its edge (only if SV_Y_QUANTILE < 1)
      trans  = sv_trans(y_brk),         # see SV_Y_TRANS at the top
      breaks = sv_breaks(y_max, y_brk),
      labels = function(x) as.character(x),
      expand = expansion(mult = c(0, 0))
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 2.5, alpha = 1),
      keywidth = 0.4, keyheight = 0.4, nrow = 1
    )) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1,
                                       size = 8.5, color = "grey10"),
      axis.text.y      = element_text(size = 8.5),
      axis.title.y     = element_text(size = 10),
      axis.title.x     = element_blank(),
      axis.line        = element_line(color = "grey20", linewidth = 0.35),
      axis.ticks       = element_line(color = "grey20", linewidth = 0.30),
      panel.grid       = element_blank(),
      plot.title       = element_text(size = 10, hjust = 0.5,
                                       margin = margin(b = 6)),
      legend.position      = "top",
      legend.justification = "right",
      legend.text          = element_text(size = 8),
      legend.spacing.x     = unit(0.15, "cm"),
      plot.margin          = margin(6, 8, 4, 4)
    ) +
    labs(
      title = paste0("snRNA-seq  –  ", comp, "  –  ",
                      tools::toTitleCase(sex)),
      y = bquote(log[2](FC))
    )
}

################################################################################
### MAIN LOOP
################################################################################

for (current_sex in c("male", "female")) {
  cat(sprintf("\n========== %s ==========\n", toupper(current_sex)))

  for (comp in COMPARISONS) {
    cat(sprintf("\n--- %s | %s ---\n", current_sex, comp))

    # ── Load ────────────────────────────────────────────────────────────────
    f <- file.path(deseq_dir,
                    paste0(current_sex, "_ALL_", comp, ".csv"))
    if (!file.exists(f)) {
      cat("  File not found — skipping:", basename(f), "\n")
      next
    }
    df <- read.csv(f, stringsAsFactors = FALSE)
    cat(sprintf("  Loaded: %d rows\n", nrow(df)))

    # ── Classify ─────────────────────────────────────────────────────────────
    df <- df %>%
      dplyr::mutate(
        sig       = !is.na(padj) & padj < PADJ_THR,
        color_cat = dplyr::case_when(
          sig & log2FoldChange > 0 ~ "Up",
          sig & log2FoldChange < 0 ~ "Down",
          TRUE                      ~ "NS"
        ),
        lfc_plot  = log2FoldChange          # no ceiling: every gene at its real log2FC
      ) %>%
      dplyr::filter(!is.na(lfc_plot))
    cat(sprintf("  log2FC range: %.2f to %.2f\n", min(df$lfc_plot), max(df$lfc_plot)))

    # DEG summary
    counts_tbl <- df %>%
      dplyr::filter(color_cat != "NS") %>%
      dplyr::count(celltype, color_cat) %>%
      tidyr::pivot_wider(names_from  = color_cat,
                          values_from = n,
                          values_fill = 0L)
    cat("  DEG summary (sig only):\n")
    print(counts_tbl, n = Inf)

    # ── Cell-type ordering ───────────────────────────────────────────────────
    present   <- intersect(CT_ORDER, unique(df$celltype))
    extra     <- sort(setdiff(unique(df$celltype), CT_ORDER))
    ct_levels <- c(present, extra)
    n_cts     <- length(ct_levels)

    # ── Plot ─────────────────────────────────────────────────────────────────
    p <- make_violin_strip(df, comp, current_sex, ct_levels)

    w <- max(7, n_cts * 0.55 + 2)
    h <- 5

    stem <- file.path(out_dir,
                       paste0(fig_date, "_", current_sex, "_", comp,
                               "_ViolinStrip"))
    save_both(p, stem, w, h)

    # ── Classified CSV ───────────────────────────────────────────────────────
    out_csv <- file.path(out_dir,
                          paste0(fig_date, "_", current_sex, "_", comp,
                                  "_classified.csv"))
    write.csv(
      df %>% dplyr::select(-any_of(c("lfc_plot", "sig",
                                      ".xoff", ".xnum"))),
      out_csv, row.names = FALSE
    )
  }
}

cat("\n=== DONE ===\nOutputs:", out_dir, "\n")

invisible(NULL)
}

################################################################################
### STEP sv_mee  Strip volcano MEEYEE
################################################################################

SECTIONS[["sv_mee"]] <- function() {

################################################################################
### StripVolcano_MEEYEE_ExcludeMHCYHC  (no log2FC ceiling)
###
### Strip volcano: Middle EE vs Young EE (MEEYEE)
### Genes sig in MHCYHC (basal middle aging) shown in grey — not coloured.
### Violin-shaped quasirandom spread (density-proportional x offset).
### Male and female separately.
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

suppressPackageStartupMessages({
  library(conflicted)
  conflicts_prefer(dplyr::select,   .quiet=TRUE)
  conflicts_prefer(dplyr::filter,   .quiet=TRUE)
  conflicts_prefer(base::intersect, .quiet=TRUE)
  conflicts_prefer(base::setdiff,   .quiet=TRUE)
  conflicts_prefer(base::union,     .quiet=TRUE)
})

USE_BEESWARM <- tryCatch({
  library(ggbeeswarm); TRUE
}, error=function(e) {
  cat("ggbeeswarm unavailable -- using manual quasirandom spread\n"); FALSE
})

################################################################################
### VIOLIN-SHAPED QUASIRANDOM HELPER
### density-proportional x spread: densest y-bin gets full MAX_W,
### sparser bins taper proportionally. FIX: .bsize/max(.bsize) not /20.
################################################################################

violin_offsets <- function(df, y_col, group_col, max_width=0.35, n_bins=80) {
  df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::mutate(
      .bin  = as.integer(cut(.data[[y_col]], breaks=n_bins,
                              include.lowest=TRUE)),
      .bin  = ifelse(is.na(.bin), 1L, .bin)) %>%
    dplyr::group_by(.data[[group_col]], .bin) %>%
    dplyr::mutate(
      .bsize = dplyr::n(),
      .idx   = seq_along(.bin)) %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::mutate(
      .max_bsize = max(.bsize),
      .xoff = ((.idx-1) - (.bsize-1)/2) /
               max(1, .bsize-1) * max_width * (.bsize / .max_bsize)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(.xoff = pmax(pmin(.xoff, max_width), -max_width))
}

################################################################################
### PATHS
################################################################################

deseq_dir <- paste0(DESEQ_DIR, "/")
fig_date  <- format(Sys.Date(),"%Y-%m-%d")
out_dir   <- paste0(file.path(FIG_ROOT, "StripVolcano_MEEYEE_ExclMHCYHC"), "/")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
cat("Output:", out_dir, "\n\n")

PADJ_THR <- 0.05
Y_FIXED  <- NULL     # NULL = axis fits this plot's own genes (see SV_Y_FIT at the top); or a number (e.g. 12) for all
MAX_W    <- 0.35

CT_ORDER <- c(
  "DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs",
  "IN",
  "Astrocytes","Oligodendrocytes","OPC","Microglia",
  "Endothelial","Pericyte","APC","Stromal",
  "Cajal-Retzius Cell","Macrophage","Epithelial","T Cell"
)

col_vals <- c(
  "EE-specific Up"   = "#D6604D",
  "EE-specific Down" = "#4393C3",
  "Also in MHCYHC"   = "grey60",
  "NS"               = "grey82"
)

################################################################################
### MAIN
################################################################################

for (current_sex in c("male","female")) {

  cat(sprintf("\n### %s ###\n", toupper(current_sex)))

  load_df <- function(comp) {
    f <- file.path(deseq_dir, paste0(current_sex,"_ALL_",comp,".csv"))
    if (!file.exists(f)) { cat("  Not found:", comp, "\n"); return(NULL) }
    df <- read.csv(f, stringsAsFactors=FALSE)
    cat(sprintf("  Loaded %-10s : %d rows\n", comp, nrow(df)))
    df
  }

  meeyee <- load_df("MEEYEE")
  mhcyhc <- load_df("MHCYHC")
  if (is.null(meeyee)) { cat("  MEEYEE missing -- skip\n"); next }

  mhcyhc_sig <- if (!is.null(mhcyhc)) {
    mhcyhc %>%
      dplyr::filter(!is.na(padj), padj < PADJ_THR) %>%
      dplyr::select(gene, celltype)
  } else {
    data.frame(gene=character(0), celltype=character(0))
  }

  df <- meeyee %>%
    dplyr::mutate(
      sig_meeyee = !is.na(padj) & padj < PADJ_THR,
      in_mhcyhc  = paste(gene, celltype) %in%
                   paste(mhcyhc_sig$gene, mhcyhc_sig$celltype),
      color_cat  = dplyr::case_when(
        sig_meeyee & !in_mhcyhc & log2FoldChange > 0 ~ "EE-specific Up",
        sig_meeyee & !in_mhcyhc & log2FoldChange < 0 ~ "EE-specific Down",
        sig_meeyee &  in_mhcyhc                       ~ "Also in MHCYHC",
        TRUE                                           ~ "NS"
      ),
      lfc_plot = log2FoldChange          # no ceiling: every gene at its real log2FC
    )
  df    <- df %>% dplyr::filter(!is.na(lfc_plot))
  Y_MAX <- if (!is.null(Y_FIXED)) Y_FIXED else sv_ymax(df$lfc_plot, df$sig_meeyee)
  Y_BRK <- sv_break_at(df$lfc_plot, df$sig_meeyee, Y_MAX)
  cat(sprintf("  log2FC range: %.2f to %.2f\n", min(df$lfc_plot), max(df$lfc_plot)))

  ct_p      <- intersect(CT_ORDER, unique(df$celltype))
  ct_e      <- sort(setdiff(unique(df$celltype), CT_ORDER))
  ct_levels <- c(ct_p, ct_e)
  df$celltype <- factor(df$celltype, levels=ct_levels)
  n_cts <- length(ct_levels)

  counts_df <- df %>%
    dplyr::filter(color_cat %in% c("EE-specific Up","EE-specific Down")) %>%
    dplyr::count(celltype, color_cat) %>%
    tidyr::pivot_wider(names_from=color_cat, values_from=n, values_fill=0L)
  cat("EE-specific DEG summary:\n"); print(counts_df)

  # ── Build plot ──────────────────────────────────────────────────────────────
  if (USE_BEESWARM) {
    p <- ggplot(mapping=aes(x=celltype, y=lfc_plot)) +
      geom_quasirandom(
        data  = df %>% dplyr::filter(color_cat=="NS"),
        color = col_vals["NS"], size=0.25, alpha=0.25, width=MAX_W) +
      geom_quasirandom(
        data  = df %>% dplyr::filter(color_cat=="Also in MHCYHC"),
        color = col_vals["Also in MHCYHC"], size=0.35, alpha=0.30, width=MAX_W) +
      geom_quasirandom(
        data  = df %>% dplyr::filter(color_cat %in%
                                      c("EE-specific Up","EE-specific Down")),
        aes(color=color_cat), size=0.65, alpha=0.80, width=MAX_W)

  } else {
    df <- violin_offsets(df, y_col="lfc_plot", group_col="celltype",
                          max_width=MAX_W, n_bins=80) %>%
      dplyr::mutate(.xnum = as.integer(celltype) + .xoff)

    p <- ggplot() +
      geom_point(
        data  = df %>% dplyr::filter(color_cat=="NS"),
        aes(x=.xnum, y=lfc_plot),
        color=col_vals["NS"], size=0.25, alpha=0.25) +
      geom_point(
        data  = df %>% dplyr::filter(color_cat=="Also in MHCYHC"),
        aes(x=.xnum, y=lfc_plot),
        color=col_vals["Also in MHCYHC"], size=0.35, alpha=0.30) +
      geom_point(
        data  = df %>% dplyr::filter(color_cat %in%
                                      c("EE-specific Up","EE-specific Down")),
        aes(x=.xnum, y=lfc_plot, color=color_cat),
        size=0.65, alpha=0.80) +
      scale_x_continuous(
        breaks = seq_along(ct_levels),
        labels = ct_levels,
        expand = expansion(add=0.6))
  }

  # ── Shared theme ────────────────────────────────────────────────────────────
  p <- p +
    geom_hline(yintercept=0, linewidth=0.5, color="grey10") +
    scale_color_manual(
      values = col_vals,
      breaks = c("EE-specific Up","EE-specific Down"),
      labels = c("up-regulated","down-regulated"),
      name   = NULL
    ) +
    sv_break_marks(Y_BRK) +
    scale_y_continuous(
      limits = c(-Y_MAX, Y_MAX),
      oob    = scales::squish,          # genes beyond the axis sit on its edge (only if SV_Y_QUANTILE < 1)
      trans  = sv_trans(Y_BRK),         # see SV_Y_TRANS at the top
      breaks = sv_breaks(Y_MAX, Y_BRK),
      labels = function(x) as.character(x),
      expand = expansion(mult=c(0,0))
    ) +
    guides(color=guide_legend(
      override.aes=list(size=2.5, alpha=1),
      keywidth=0.4, keyheight=0.4, nrow=1)) +
    theme_classic(base_size=10) +
    theme(
      axis.text.x          = element_text(angle=45, hjust=1, size=8.5,
                                           color="grey10"),
      axis.text.y          = element_text(size=8.5),
      axis.title.y         = element_text(size=10),
      axis.title.x         = element_blank(),
      axis.line            = element_line(color="grey20", linewidth=0.35),
      axis.ticks           = element_line(color="grey20", linewidth=0.3),
      panel.grid           = element_blank(),
      plot.title           = element_text(size=10, hjust=0.5,
                                           margin=margin(b=4)),
      plot.subtitle        = element_text(size=8, hjust=0.5, color="grey40",
                                           margin=margin(b=6)),
      legend.position      = "top",
      legend.justification = "right",
      legend.text          = element_text(size=8),
      legend.spacing.x     = unit(0.15,"cm"),
      plot.margin          = margin(6, 8, 4, 4)
    ) +
    labs(
      title    = paste0("snRNA-seq  --  MEEYEE  --  ",
                         tools::toTitleCase(current_sex)),
      subtitle = paste0("Genes sig in MHCYHC (basal middle aging) shown in grey",
                         " | padj<", PADJ_THR),
      y = bquote(log[2](FC))
    )

  w <- max(7, n_cts*0.52+2)
  h <- 5

  out_pdf <- file.path(out_dir,
                        paste0(fig_date,"_",current_sex,
                               "_MEEYEE_ExcludeMHCYHC_StripVolcano.pdf"))
  out_png <- file.path(out_dir,
                        paste0(fig_date,"_",current_sex,
                               "_MEEYEE_ExcludeMHCYHC_StripVolcano.png"))

  ggsave(out_pdf, p, width=w, height=h, device="pdf")
  cat("Saved:", basename(out_pdf), "\n")

  ggsave(out_png, p, width=w, height=h, device="png", dpi=300)
  cat("Saved:", basename(out_png), "\n")

  write.csv(
    df %>% dplyr::select(-any_of(c("lfc_plot",".bin",".max_bsize",
                                    ".bsize",".idx",".xoff",".xnum"))),
    file.path(out_dir, paste0(fig_date,"_",current_sex,
                               "_MEEYEE_classified.csv")),
    row.names=FALSE)
}

cat("\n=== DONE ===\nOutputs:", out_dir, "\n")

invisible(NULL)
}

################################################################################
### STEP sv_oee  Strip volcano OEEYEE
################################################################################

SECTIONS[["sv_oee"]] <- function() {

################################################################################
### StripVolcano_OEEYEE_ExcludeOHCYHC  (no log2FC ceiling: all genes shown at their real value)
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

suppressPackageStartupMessages({
  library(conflicted)
  conflicts_prefer(dplyr::select,   .quiet=TRUE)
  conflicts_prefer(dplyr::filter,   .quiet=TRUE)
  conflicts_prefer(base::intersect, .quiet=TRUE)
  conflicts_prefer(base::setdiff,   .quiet=TRUE)
  conflicts_prefer(base::union,     .quiet=TRUE)
})

# Try ggbeeswarm; fall back to manual quasirandom if corrupt/missing
USE_BEESWARM <- tryCatch({
  library(ggbeeswarm)
  TRUE
}, error = function(e) {
  cat("ggbeeswarm unavailable -- using manual quasirandom spread\n")
  FALSE
})

################################################################################
### MANUAL QUASIRANDOM (mimics geom_quasirandom density spread)
### Groups y into bins per cell type, then offsets x by rank within bin.
################################################################################

quasirandom_offsets <- function(df, y_col, group_col,
                                 max_width=0.35, n_bins=100) {
  df <- df %>% dplyr::mutate(.row = dplyr::row_number())
  df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::mutate(
      .bin   = cut(.data[[y_col]], breaks=n_bins, labels=FALSE,
                   include.lowest=TRUE),
      .bin   = ifelse(is.na(.bin), 1L, .bin),
      .rank  = dplyr::row_number(.bin),        # rank within each bin
      .n_bin = dplyr::n(),                      # total in this group
      # bin size = how many points share this bin
      .bsize = ave(rep(1, dplyr::n()), .bin, FUN=length),
      # center-offset within bin
      .idx   = ave(seq_along(.bin), .bin, FUN=seq_along),
      .xoff  = ((.idx - 1) - (.bsize - 1) / 2) /
               max(1, max(.bsize - 1)) * max_width * (.bsize / max(.bsize))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::pull(.xoff)
}

################################################################################
### PATHS
################################################################################

deseq_dir <- paste0(DESEQ_DIR, "/")
fig_date  <- format(Sys.Date(),"%Y-%m-%d")
out_dir   <- paste0(file.path(FIG_ROOT, "StripVolcano_OEEYEE_ExclOHCYHC"), "/")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
cat("Output:", out_dir, "\n\n")

PADJ_THR <- 0.05
Y_FIXED  <- NULL     # NULL = axis fits this plot's own genes (see SV_Y_FIT at the top); or a number (e.g. 12) for all

CT_ORDER <- c(
  "DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs",
  "IN",
  "Astrocytes","Oligodendrocytes","OPC","Microglia",
  "Endothelial","Pericyte","APC","Stromal",
  "Cajal-Retzius Cell","Macrophage","Epithelial","T Cell"
)

col_vals <- c(
  "EE-specific Up"   = "#D6604D",
  "EE-specific Down" = "#4393C3",
  "Also in OHCYHC"   = "grey60",
  "NS"               = "grey82"
)

################################################################################
### MAIN
################################################################################

for (current_sex in c("male","female")) {

  cat(sprintf("\n### %s ###\n", toupper(current_sex)))

  load_df <- function(comp) {
    f <- file.path(deseq_dir, paste0(current_sex,"_ALL_",comp,".csv"))
    if (!file.exists(f)) { cat("  Not found:", comp, "\n"); return(NULL) }
    df <- read.csv(f, stringsAsFactors=FALSE)
    cat(sprintf("  Loaded %-10s : %d rows\n", comp, nrow(df)))
    df
  }

  oeeyee <- load_df("OEEYEE")
  ohcyhc <- load_df("OHCYHC")
  if (is.null(oeeyee)) { cat("  OEEYEE missing -- skip\n"); next }

  ohcyhc_sig <- if (!is.null(ohcyhc)) {
    ohcyhc %>%
      dplyr::filter(!is.na(padj), padj < PADJ_THR) %>%
      dplyr::select(gene, celltype)
  } else {
    data.frame(gene=character(0), celltype=character(0))
  }

  df <- oeeyee %>%
    dplyr::mutate(
      sig_oeeyee = !is.na(padj) & padj < PADJ_THR,
      in_ohcyhc  = paste(gene, celltype) %in%
                   paste(ohcyhc_sig$gene, ohcyhc_sig$celltype),
      color_cat  = dplyr::case_when(
        sig_oeeyee & !in_ohcyhc & log2FoldChange > 0  ~ "EE-specific Up",
        sig_oeeyee & !in_ohcyhc & log2FoldChange < 0  ~ "EE-specific Down",
        sig_oeeyee &  in_ohcyhc                        ~ "Also in OHCYHC",
        TRUE                                            ~ "NS"
      ),
      lfc_plot = log2FoldChange          # no ceiling: every gene at its real log2FC
    )
  df    <- df %>% dplyr::filter(!is.na(lfc_plot))
  Y_MAX <- if (!is.null(Y_FIXED)) Y_FIXED else sv_ymax(df$lfc_plot, df$sig_oeeyee)
  Y_BRK <- sv_break_at(df$lfc_plot, df$sig_oeeyee, Y_MAX)
  cat(sprintf("  log2FC range: %.2f to %.2f\n", min(df$lfc_plot), max(df$lfc_plot)))

  ct_p      <- intersect(CT_ORDER, unique(df$celltype))
  ct_e      <- sort(setdiff(unique(df$celltype), CT_ORDER))
  ct_levels <- c(ct_p, ct_e)
  df$celltype <- factor(df$celltype, levels=ct_levels)
  n_cts <- length(ct_levels)

  # Summary
  counts_df <- df %>%
    dplyr::filter(color_cat %in% c("EE-specific Up","EE-specific Down")) %>%
    dplyr::count(celltype, color_cat) %>%
    tidyr::pivot_wider(names_from=color_cat, values_from=n, values_fill=0L)
  cat("EE-specific DEG summary:\n"); print(counts_df)

  # ── Build plot ────────────────────────────────────────────────────────────
  if (USE_BEESWARM) {
    # ── ggbeeswarm path ──────────────────────────────────────────────────────
    p <- ggplot(mapping = aes(x=celltype, y=lfc_plot)) +
      geom_quasirandom(
        data    = df %>% dplyr::filter(color_cat == "NS"),
        color=col_vals["NS"], size=0.25, alpha=0.25, width=0.35) +
      geom_quasirandom(
        data    = df %>% dplyr::filter(color_cat == "Also in OHCYHC"),
        color=col_vals["Also in OHCYHC"], size=0.35, alpha=0.30, width=0.35) +
      geom_quasirandom(
        data    = df %>% dplyr::filter(color_cat %in%
                                        c("EE-specific Up","EE-specific Down")),
        aes(color=color_cat), size=0.65, alpha=0.80, width=0.35)

  } else {
    # ── Manual quasirandom path ──────────────────────────────────────────────
    MAX_W <- 0.35
    df <- df %>%
      dplyr::group_by(celltype) %>%
      dplyr::mutate(
        .bin   = as.integer(cut(lfc_plot, breaks=80,
                                include.lowest=TRUE)),
        .bin   = ifelse(is.na(.bin), 1L, .bin)) %>%
      dplyr::group_by(celltype, .bin) %>%
      dplyr::mutate(
        .bsize = dplyr::n(),
        .idx   = seq_along(.bin),
        .xoff  = ((.idx-1) - (.bsize-1)/2) /
                  max(1,.bsize-1) * MAX_W * (.bsize / 20)) %>%
      dplyr::ungroup()

    # Clamp offsets
    df$.xoff <- pmax(pmin(df$.xoff, MAX_W), -MAX_W)

    # Convert celltype to numeric for manual x
    df$.xnum <- as.integer(df$celltype) + df$.xoff

    plot_layers <- function(sub_df, col, sz, al) {
      geom_point(data=sub_df,
                 aes(x=.xnum, y=lfc_plot),
                 color=col, size=sz, alpha=al)
    }

    p <- ggplot() +
      plot_layers(df %>% dplyr::filter(color_cat=="NS"),
                  col_vals["NS"], 0.25, 0.25) +
      plot_layers(df %>% dplyr::filter(color_cat=="Also in OHCYHC"),
                  col_vals["Also in OHCYHC"], 0.35, 0.30) +
      geom_point(
        data  = df %>% dplyr::filter(color_cat %in%
                                      c("EE-specific Up","EE-specific Down")),
        aes(x=.xnum, y=lfc_plot, color=color_cat),
        size=0.65, alpha=0.80) +
      scale_x_continuous(
        breaks = seq_along(ct_levels),
        labels = ct_levels,
        expand = expansion(add=0.6))
  }

  # ── Shared theme layers ───────────────────────────────────────────────────
  p <- p +
    geom_hline(yintercept=0, linewidth=0.5, color="grey10") +
    scale_color_manual(
      values = col_vals,
      breaks = c("EE-specific Up","EE-specific Down","NS"),
      labels = c("up-regulated","down-regulated","not significant"),
      name   = NULL
    ) +
    sv_break_marks(Y_BRK) +
    scale_y_continuous(
      limits = c(-Y_MAX, Y_MAX),
      oob    = scales::squish,          # genes beyond the axis sit on its edge (only if SV_Y_QUANTILE < 1)
      trans  = sv_trans(Y_BRK),         # see SV_Y_TRANS at the top
      breaks = sv_breaks(Y_MAX, Y_BRK),
      labels = function(x) as.character(x),
      expand = expansion(mult=c(0,0))
    ) +
    guides(color = guide_legend(
      override.aes = list(size=2.5, alpha=1),
      keywidth=0.4, keyheight=0.4, nrow=1)) +
    theme_classic(base_size=10) +
    theme(
      axis.text.x        = element_text(angle=45, hjust=1, size=8.5,
                                         color="grey10"),
      axis.text.y        = element_text(size=8.5),
      axis.title.y       = element_text(size=10),
      axis.title.x       = element_blank(),
      axis.line          = element_line(color="grey20", linewidth=0.35),
      axis.ticks         = element_line(color="grey20", linewidth=0.3),
      panel.grid         = element_blank(),
      plot.title         = element_text(size=10, hjust=0.5,
                                         margin=margin(b=6)),
      legend.position    = "top",
      legend.justification = "right",
      legend.text        = element_text(size=8),
      legend.spacing.x   = unit(0.15,"cm"),
      plot.margin        = margin(6, 8, 4, 4)
    ) +
    labs(
      title = paste0("snRNA-seq  –  OEEYEE  –  ",
                      tools::toTitleCase(current_sex)),
      y = bquote(log[2](FC))
    )

  w <- max(7, n_cts*0.52+2)
  h <- 5

  out_pdf <- file.path(out_dir,
                        paste0(fig_date,"_",current_sex,
                               "_OEEYEE_ExcludeOHCYHC_StripVolcano.pdf"))
  out_png <- file.path(out_dir,
                        paste0(fig_date,"_",current_sex,
                               "_OEEYEE_ExcludeOHCYHC_StripVolcano.png"))

  ggsave(out_pdf, p, width=w, height=h, device="pdf")
  cat("Saved:", basename(out_pdf), "\n")

  ggsave(out_png, p, width=w, height=h, device="png", dpi=300)
  cat("Saved:", basename(out_png), "\n")

  write.csv(df %>% dplyr::select(-any_of(c("lfc_plot",".bin",
                                             ".bsize",".idx",".xoff",".xnum"))),
            file.path(out_dir, paste0(fig_date,"_",current_sex,
                                       "_OEEYEE_classified.csv")),
            row.names=FALSE)
}

cat("\n=== DONE ===\nOutputs:", out_dir, "\n")

invisible(NULL)
}

################################################################################
### STEP ed4cd  ED Fig. 4c,d YEEYHC pathways
################################################################################

SECTIONS[["ed4cd"]] <- function() {
################################################################################
### Pathway_SankeyDot_YEEYHC
###
### GO BP + KEGG enrichment (Sankey+Dot / SR-plot style) for YEEYHC,
### both sexes.
###
### Per comparison × sex:
###   1. Per cell type   — GO BP + KEGG sankey+dot for Up and Down sig genes
###   2. Pooled          — compareCluster dot + cluster_dot_plot (f/h style)
###   3. Super-pooled    — all CTs merged, GO BP + KEGG sankey+dot
###
### No GSEA / ridgeplot section.
###
### Output: Figures/Figure2/ED4cd_Pathway_YEEYHC/
################################################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggalluvial)
  library(patchwork)
  library(conflicted)
  conflicts_prefer(dplyr::select,           .quiet = TRUE)
  conflicts_prefer(dplyr::filter,           .quiet = TRUE)
  conflicts_prefer(base::intersect,         .quiet = TRUE)
  conflicts_prefer(base::setdiff,           .quiet = TRUE)
  conflicts_prefer(base::union,             .quiet = TRUE)
  conflicts_prefer(AnnotationDbi::select,   .quiet = TRUE)
  conflicts_prefer(base::as.data.frame,     .quiet = TRUE)
})

################################################################################
### PATHS & CONSTANTS
################################################################################

deseq_dir <- DESEQ_DIR
fig_date  <- format(Sys.Date(), "%Y-%m-%d")
out_dir   <- file.path(FIG_ROOT, "ED4cd_Pathway_YEEYHC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output dir:", out_dir, "\n\n")

PADJ_THR  <- 0.05
Q_THR     <- 0.2
MIN_GENES <- 5
TOP_N     <- 10   # pathway terms shown per sankey+dot
TOP_N_CMP          <- 8    # terms shown in compareCluster dot
TOP_N_CLUSTER      <- 5    # terms shown per cell type in cluster_dot_plot
MAX_GENES_PER_PATH <- 6   # cap genes shown per pathway in Sankey (5=tight, 8=looser)

COMP <- "YEEYHC" ###Can change this to whichever comparison I want or all together

################################################################################
### HELPERS
################################################################################

should_skip <- function(f) {
  if (file.exists(f)) { cat("  Skipping:", basename(f), "\n"); return(TRUE) }
  FALSE
}

# Save PDF + 600 dpi PNG
save_both <- function(p, stem, w, h) {
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h,
         device = "pdf", limitsize = FALSE)
  ggsave(paste0(stem, ".png"), p, width = w, height = h,
         device = "png", dpi = 600, limitsize = FALSE)
  cat("  Saved:", basename(stem), "[pdf+png]\n")
}

pathway_palette <- function(n) {
  base_cols <- c(
    "#E07B54","#5B8DB8","#6BAF6B","#8F7BB5","#C0B45A",
    "#7BBFBF","#C97BAB","#A07850","#7090B8","#90B870",
    "#B07070","#70A0A0","#A07090","#90A070","#A09060"
  )
  if (n <= length(base_cols)) return(base_cols[seq_len(n)])
  scales::hue_pal()(n)
}

################################################################################
### ENRICHMENT FUNCTIONS
################################################################################

run_go <- function(genes, universe, ont = "BP") {
  genes <- intersect(genes, universe)
  if (length(genes) < MIN_GENES) return(NULL)
  tryCatch({
    res <- enrichGO(
      gene          = genes,
      universe      = universe,
      OrgDb         = org.Mm.eg.db,
      keyType       = "SYMBOL",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = Q_THR,
      readable      = TRUE
    )
    if (is.null(res) || nrow(res) == 0) return(NULL)
    clusterProfiler::simplify(res, cutoff = 0.7, by = "p.adjust")
  }, error = function(e) NULL)
}

sym_to_entrez <- function(genes) {
  res <- AnnotationDbi::mapIds(
    org.Mm.eg.db, keys = genes,
    column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"
  )
  as.character(res[!is.na(res)])
}

run_kegg <- function(genes, universe) {
  genes <- intersect(genes, universe)
  if (length(genes) < MIN_GENES) return(NULL)
  eg <- sym_to_entrez(genes)
  eu <- sym_to_entrez(universe)
  if (length(eg) < MIN_GENES) return(NULL)
  tryCatch({
    res <- enrichKEGG(
      gene          = eg,
      universe      = eu,
      organism      = "mmu",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = Q_THR
    )
    if (is.null(res) || nrow(res) == 0) return(NULL)
    setReadable(res, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  }, error = function(e) NULL)
}

################################################################################
### SANKEY + DOT PLOT  (SR-plot style)
###
### Left panel : gene symbols → flow (colored by pathway) → pathway terms
### Right panel: dot plot (x=GeneRatio, y=pathway, size=Count,
###              color=-log10 p.adjust, thin color strip on left edge)
################################################################################

sankey_dot_plot <- function(enrich_res, title, out_stem,
                             top_n              = TOP_N,
                             max_genes_per_path = MAX_GENES_PER_PATH,
                             width              = 14,
                             height             = NULL) {

  if (is.null(enrich_res) || nrow(enrich_res) == 0) return(invisible(NULL))
  if (should_skip(paste0(out_stem, ".pdf"))) return(invisible(NULL))

  tryCatch({
    # ── 1. Pathway-level data ─────────────────────────────────────────────────
    df <- as.data.frame(enrich_res) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::mutate(
        GeneRatio_num  = sapply(GeneRatio, function(x) {
          v <- strsplit(x, "/")[[1]]
          as.numeric(v[1]) / as.numeric(v[2])
        }),
        neg_log10_padj = -log10(p.adjust),
        Desc_short     = ifelse(nchar(Description) > 45,
                                paste0(substr(Description, 1, 42), "..."),
                                Description)
      ) %>%
      # Disambiguate any truncated descriptions that collide
      dplyr::group_by(Desc_short) %>%
      dplyr::mutate(Desc_short = if (dplyr::n() > 1)
                      paste0(Desc_short, " (", dplyr::row_number(), ")")
                    else Desc_short) %>%
      dplyr::ungroup()

    # Pathway factor: most significant term at top
    pathway_levels <- rev(df$Desc_short)
    df$Desc_short  <- factor(df$Desc_short, levels = pathway_levels)

    # ── 2. Gene-to-pathway edges ──────────────────────────────────────────────
    edge_df <- df %>%
      dplyr::select(Desc_short, geneID, neg_log10_padj) %>%
      dplyr::mutate(genes = strsplit(as.character(geneID), "/")) %>%
      tidyr::unnest(genes) %>%
      dplyr::rename(Gene = genes, Pathway = Desc_short) %>%
      dplyr::mutate(Pathway = factor(Pathway, levels = levels(df$Desc_short))) %>%
      # cap representative genes per pathway to reduce Sankey density
      dplyr::group_by(Pathway) %>%
      dplyr::slice_head(n = max_genes_per_path) %>%
      dplyr::ungroup()

    # ── 3. Colors ─────────────────────────────────────────────────────────────
    n_paths     <- nlevels(df$Desc_short)
    path_colors <- setNames(pathway_palette(n_paths), levels(df$Desc_short))

    # ── 4. Gene order by most-significant pathway ─────────────────────────────
    gene_rank    <- edge_df %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(top_rank = min(as.integer(Pathway)), .groups = "drop") %>%
      dplyr::arrange(top_rank, Gene)
    gene_levels  <- unique(gene_rank$Gene)
    edge_df$Gene <- factor(edge_df$Gene, levels = gene_levels)
    n_genes      <- length(gene_levels)

    # ── 5. Left panel: Sankey ─────────────────────────────────────────────────
    p_sankey <- ggplot(
        edge_df,
        aes(axis1 = Gene, axis2 = Pathway, y = 1)) +
      ggalluvial::geom_alluvium(
        aes(fill = Pathway),
        width      = 0.25,
        alpha      = 0.50,
        curve_type = "sigmoid",
        knot.pos   = 0.4
      ) +
      ggalluvial::geom_stratum(
        aes(fill = Pathway),
        width     = 0.25,
        color     = "grey70",
        linewidth = 0.25
      ) +
      geom_text(
        stat     = "stratum",
        aes(label = after_stat(stratum)),
        size     = 2.6,
        hjust    = 1,
        nudge_x  = -0.15,
        color    = "grey15",
        fontface = "italic"
      ) +
      scale_fill_manual(values = path_colors, guide = "none") +
      scale_x_discrete(
        limits = c("Gene", "Pathway"),
        expand = c(0.35, 0.05)
      ) +
      theme_void(base_size = 9) +
      theme(
        axis.text.x = element_text(size = 9, face = "bold", color = "grey30",
                                   margin = margin(t = 4)),
        plot.margin = margin(10, 0, 10, 65)
      )

    # ── 6. Right panel: Dot plot ──────────────────────────────────────────────
    dot_xmin <- min(df$GeneRatio_num) * 0.75
    strip_w  <- (max(df$GeneRatio_num) - min(df$GeneRatio_num)) * 0.04

    p_dot <- ggplot(df, aes(x = GeneRatio_num, y = Desc_short)) +
      geom_tile(
        aes(x = dot_xmin - strip_w / 2, fill = Desc_short),
        width = strip_w, height = 0.75, color = NA
      ) +
      geom_point(aes(color = neg_log10_padj, size = Count), alpha = 0.9) +
      scale_color_gradient(
        low  = "#3182BD",
        high = "#E6550D",
        name = bquote(-log[10](italic(p)[adj]))
      ) +
      scale_fill_manual(values = path_colors, guide = "none") +
      scale_size_continuous(range = c(2, 9), name = "Count") +
      scale_x_continuous(
        expand = expansion(mult = c(0.15, 0.18)),
        labels = scales::number_format(accuracy = 0.01)
      ) +
      theme_classic(base_size = 9) +
      theme(
        axis.text.y     = element_text(size = 8, color = "grey15"),
        axis.ticks.y    = element_blank(),
        axis.line.y     = element_blank(),
        axis.text.x     = element_text(size = 7.5),
        axis.title.x    = element_text(size = 8.5),
        plot.margin     = margin(10, 12, 10, 4),
        legend.position = "right",
        legend.text     = element_text(size = 7),
        legend.title    = element_text(size = 7.5)
      ) +
      labs(x = "Gene Ratio", y = NULL)

    # ── 7. Combine ────────────────────────────────────────────────────────────
    auto_h  <- if (is.null(height)) max(5, n_genes * 0.22 + 2) else height
    combined <- p_sankey + p_dot +
      patchwork::plot_layout(widths = c(2.2, 1)) +
      patchwork::plot_annotation(
        title    = title,
        subtitle = paste0("Top ", nrow(df), " terms | ordered by p.adjust",
                          " | max ", max_genes_per_path, " genes/pathway shown"),
        theme    = theme(
          plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle   = element_text(size = 8.5, color = "grey40",
                                         hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA)
        )
      )

    save_both(combined, out_stem, width, auto_h)
    write.csv(as.data.frame(enrich_res),
              paste0(out_stem, ".csv"), row.names = FALSE)

  }, error = function(e) cat("  sankey_dot error:", e$message, "\n"))
}

################################################################################
### COMPARE CLUSTER DOT PLOT (pooled across cell types — classic x=CT style)
################################################################################

compare_dot <- function(gene_lists, universe, ont = "BP",
                         title, out_stem, top_n = TOP_N_CMP,
                         width = 16, height = 10) {

  gene_lists <- gene_lists[lengths(gene_lists) >= MIN_GENES]
  if (length(gene_lists) < 2) {
    cat("  Too few gene sets for compareCluster\n")
    return(invisible(NULL))
  }
  if (should_skip(paste0(out_stem, ".pdf"))) return(invisible(NULL))

  tryCatch({
    res <- compareCluster(
      geneClusters  = gene_lists,
      fun           = "enrichGO",
      OrgDb         = org.Mm.eg.db,
      keyType       = "SYMBOL",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      universe      = universe,
      readable      = TRUE
    )
    if (is.null(res) || nrow(res) == 0) {
      cat("  No significant terms in compareCluster\n")
      return(invisible(NULL))
    }

    write.csv(as.data.frame(res), paste0(out_stem, ".csv"), row.names = FALSE)

    df <- as.data.frame(res) %>%
      dplyr::group_by(Cluster) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        GeneRatio_num  = sapply(GeneRatio, function(x) {
          v <- strsplit(x, "/")[[1]]
          as.numeric(v[1]) / as.numeric(v[2])
        }),
        neg_log10_padj = -log10(p.adjust),
        Description    = factor(
          Description,
          levels = rev(unique(
            (dplyr::arrange(., p.adjust))$Description
          ))
        )
      )

    n_cls <- length(unique(df$Cluster))

    p <- ggplot(df, aes(x = Cluster, y = Description,
                         color = neg_log10_padj, size = GeneRatio_num)) +
      geom_point(alpha = 0.9) +
      scale_color_gradient(
        low  = "#3182BD",
        high = "#E6550D",
        name = bquote(-log[10](italic(p)[adj]))
      ) +
      scale_size_continuous(range = c(1.5, 7), name = "Gene ratio") +
      theme_classic(base_size = 8) +
      theme(
        axis.text.x     = element_text(angle = 35, hjust = 1, size = 7.5),
        axis.text.y     = element_text(size = 7.5),
        axis.title      = element_blank(),
        legend.position = "right",
        legend.text     = element_text(size = 7),
        legend.title    = element_text(size = 7.5),
        plot.margin     = margin(5, 10, 5, 5)
      ) +
      patchwork::plot_annotation(
        title    = title,
        subtitle = paste0("Top ", top_n, " terms per gene set | GO ", ont),
        theme    = theme(
          plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle   = element_text(size = 8.5, color = "grey40",
                                         hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA)
        )
      )

    save_both(p, out_stem,
              w = max(width, n_cls * 1.0 + 3.5),
              h = height)
    return(res)

  }, error = function(e) {
    cat("  compare_dot error:", e$message, "\n"); NULL
  })
}

################################################################################
### CLUSTER DOT PLOT  (panels f/h style — y=term, x=-log10p, facet=cell type)
###
### Takes the compareCluster result object from compare_dot().
### Each cell type appears as a labelled facet on the LEFT; top TOP_N_CLUSTER
### terms per CT are shown.  Dots: size = Count, fill = -log10(p.adjust).
################################################################################

cluster_dot_plot <- function(compare_res, title, out_stem,
                              top_n  = TOP_N_CLUSTER,
                              width  = 11,
                              height = NULL) {

  if (is.null(compare_res) || nrow(compare_res) == 0) return(invisible(NULL))
  if (should_skip(paste0(out_stem, ".pdf"))) return(invisible(NULL))

  tryCatch({
    df <- as.data.frame(compare_res) %>%
      dplyr::group_by(Cluster) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        neg_log_padj = -log10(p.adjust),
        Count        = as.numeric(Count),
        term_short   = ifelse(
          nchar(Description) > 52,
          paste0(substr(Description, 1, 49), "..."),
          Description
        )
      )

    if (nrow(df) == 0) return(invisible(NULL))

    # Preserve cluster order as it came from compareCluster
    cluster_levels <- as.character(unique(df$Cluster))
    n_cls          <- length(cluster_levels)
    cls_cols       <- setNames(pathway_palette(n_cls), cluster_levels)

    df <- df %>%
      dplyr::mutate(Cluster = factor(Cluster, levels = cluster_levels)) %>%
      dplyr::arrange(Cluster, p.adjust) %>%
      # Within each cluster, most significant term on top → reverse for ggplot
      dplyr::group_by(Cluster) %>%
      dplyr::mutate(rank_in_ct = dplyr::n() - dplyr::row_number()) %>%
      dplyr::ungroup()

    # Create a unique y-key per cluster × term (avoids label collisions across CTs)
    df <- df %>%
      dplyr::mutate(
        y_key = paste0(Cluster, "___", rank_in_ct, "___", term_short),
        y_key = factor(y_key, levels = rev(unique(y_key)))
      )

    n_terms <- nrow(df)
    auto_h  <- if (is.null(height)) max(5, n_terms * 0.30 + 2) else height

    # ── Build plot ────────────────────────────────────────────────────────────
    p <- ggplot(df, aes(x = neg_log_padj, y = y_key)) +

      # Main enrichment dots
      geom_point(
        aes(size = Count, fill = neg_log_padj),
        shape = 21, color = "white", stroke = 0.25, alpha = 0.95
      ) +

      # Facet by cell type; strip on the LEFT, coloured per CT
      facet_grid(
        Cluster ~ .,
        scales = "free_y",
        space  = "free_y",
        switch = "y"
      ) +

      scale_fill_gradient(
        low  = "#F4B860",
        high = "#C0392B",
        name = bquote(-log[10](italic(p)[adj]))
      ) +
      scale_size_continuous(range = c(2, 8), name = "# genes") +

      # Strip y-axis labels to just the term name (drop the CT___rank___ prefix)
      scale_y_discrete(
        labels = function(x) sub("^[^_]+___[^_]+___", "", x)
      ) +

      # x-axis: start at 0
      scale_x_continuous(
        expand = expansion(mult = c(0, 0.15)),
        limits = c(0, NA)
      ) +

      theme_classic(base_size = 9) +
      theme(
        # Left strip = cell type label box
        strip.placement      = "outside",
        strip.background     = element_rect(fill = "grey88", color = NA),
        strip.text.y.left    = element_text(
          angle = 0, size = 7, face = "bold",
          hjust = 1, margin = margin(r = 4)
        ),
        # y-axis term labels
        axis.text.y          = element_text(size = 7.5, color = "grey15"),
        axis.ticks.y         = element_blank(),
        axis.line.y          = element_blank(),
        # x-axis
        axis.text.x          = element_text(size = 7.5),
        axis.title.x         = element_text(size = 8.5),
        # spacing between facet panels
        panel.spacing.y      = unit(3, "pt"),
        # legends
        legend.position      = "right",
        legend.text          = element_text(size = 7),
        legend.title         = element_text(size = 7.5),
        # title / subtitle
        plot.title           = element_text(face = "bold", size = 12,
                                             hjust = 0.5),
        plot.subtitle        = element_text(size = 8.5, color = "grey40",
                                             hjust = 0.5),
        plot.background      = element_rect(fill = "white", color = NA)
      ) +
      labs(
        title    = title,
        subtitle = paste0(
          "Top ", top_n, " GO BP terms per cell type | ",
          "ordered by p.adjust | dot size = # genes"
        ),
        x = bquote(-log[10](italic(p)[adj])),
        y = NULL
      )

    save_both(p, out_stem, width, auto_h)

  }, error = function(e) cat("  cluster_dot error:", e$message, "\n"))
}

################################################################################
### MAIN LOOP
################################################################################

for (current_sex in c("male", "female")) {
  cat(sprintf("\n========== %s ==========\n", toupper(current_sex)))

  for (comp in COMP) {
    cat(sprintf("\n--- %s | %s ---\n", current_sex, comp))

    # ── Load ──────────────────────────────────────────────────────────────────
    f <- file.path(deseq_dir,
                    paste0(current_sex, "_ALL_", comp, ".csv"))
    if (!file.exists(f)) {
      cat("  File not found — skipping:", basename(f), "\n")
      next
    }
    all_df    <- read.csv(f, stringsAsFactors = FALSE)
    universe  <- unique(all_df$gene[!is.na(all_df$gene)])
    celltypes <- sort(unique(all_df$celltype))
    cat(sprintf("  Universe: %d genes | %d cell types\n",
                length(universe), length(celltypes)))

    # Output dirs
    comp_out <- file.path(out_dir, current_sex, comp)
    dir.create(comp_out, recursive = TRUE, showWarnings = FALSE)

    # Accumulate for pooled compareCluster
    pool_up <- list()
    pool_dn <- list()

    ##########################################################
    ### 1. PER CELL TYPE
    ##########################################################
    cat("  Per cell type enrichment:\n")

    for (ct in celltypes) {
      ct_df  <- all_df %>% dplyr::filter(celltype == ct)
      sig_up <- ct_df$gene[!is.na(ct_df$padj) & ct_df$padj < PADJ_THR &
                              ct_df$log2FoldChange > 0]
      sig_dn <- ct_df$gene[!is.na(ct_df$padj) & ct_df$padj < PADJ_THR &
                              ct_df$log2FoldChange < 0]

      cat(sprintf("    %-25s up=%d  down=%d\n", ct,
                  length(sig_up), length(sig_dn)))

      ct_slug <- gsub("[^A-Za-z0-9]+", "_", ct)
      ct_out  <- file.path(comp_out, ct_slug)
      dir.create(ct_out, recursive = TRUE, showWarnings = FALSE)

      if (length(sig_up) >= MIN_GENES) pool_up[[ct]] <- sig_up
      if (length(sig_dn) >= MIN_GENES) pool_dn[[ct]] <- sig_dn

      # UP —————————————————————————————————————————————
      if (length(sig_up) >= MIN_GENES) {
        # GO BP
        r_go <- run_go(sig_up, universe, "BP")
        sankey_dot_plot(
          r_go,
          title    = paste0(ct, " | UP | GO BP | ", comp,
                             " | ", tools::toTitleCase(current_sex)),
          out_stem = file.path(ct_out,
                                paste0(fig_date, "_", ct_slug,
                                       "_UP_GO_BP_SankeyDot"))
        )
        # KEGG
        r_kk <- run_kegg(sig_up, universe)
        sankey_dot_plot(
          r_kk,
          title    = paste0(ct, " | UP | KEGG | ", comp,
                             " | ", tools::toTitleCase(current_sex)),
          out_stem = file.path(ct_out,
                                paste0(fig_date, "_", ct_slug,
                                       "_UP_KEGG_SankeyDot"))
        )
      }

      # DOWN ————————————————————————————————————————————
      if (length(sig_dn) >= MIN_GENES) {
        # GO BP
        r_go <- run_go(sig_dn, universe, "BP")
        sankey_dot_plot(
          r_go,
          title    = paste0(ct, " | DOWN | GO BP | ", comp,
                             " | ", tools::toTitleCase(current_sex)),
          out_stem = file.path(ct_out,
                                paste0(fig_date, "_", ct_slug,
                                       "_DOWN_GO_BP_SankeyDot"))
        )
        # KEGG
        r_kk <- run_kegg(sig_dn, universe)
        sankey_dot_plot(
          r_kk,
          title    = paste0(ct, " | DOWN | KEGG | ", comp,
                             " | ", tools::toTitleCase(current_sex)),
          out_stem = file.path(ct_out,
                                paste0(fig_date, "_", ct_slug,
                                       "_DOWN_KEGG_SankeyDot"))
        )
      }
    }

    ##########################################################
    ### 2. POOLED COMPARE CLUSTER (across cell types)
    ##########################################################
    cat("\n  Pooled compareCluster:\n")

    # UP — classic x=CT dot + new cluster_dot_plot style
    cmp_up <- compare_dot(
      gene_lists = pool_up,
      universe   = universe,
      ont        = "BP",
      title      = paste0("GO BP — UP genes per cell type | ",
                           tools::toTitleCase(current_sex), " | ", comp),
      out_stem   = file.path(comp_out,
                              paste0(fig_date, "_", current_sex, "_", comp,
                                     "_Pooled_UP_compareCluster_Dot"))
    )
    cluster_dot_plot(
      cmp_up,
      title    = paste0("GO BP — UP genes per cell type | ",
                         tools::toTitleCase(current_sex), " | ", comp),
      out_stem = file.path(comp_out,
                            paste0(fig_date, "_", current_sex, "_", comp,
                                   "_Pooled_UP_ClusterDot"))
    )

    # DOWN — classic x=CT dot + new cluster_dot_plot style
    cmp_dn <- compare_dot(
      gene_lists = pool_dn,
      universe   = universe,
      ont        = "BP",
      title      = paste0("GO BP — DOWN genes per cell type | ",
                           tools::toTitleCase(current_sex), " | ", comp),
      out_stem   = file.path(comp_out,
                              paste0(fig_date, "_", current_sex, "_", comp,
                                     "_Pooled_DOWN_compareCluster_Dot"))
    )
    cluster_dot_plot(
      cmp_dn,
      title    = paste0("GO BP — DOWN genes per cell type | ",
                         tools::toTitleCase(current_sex), " | ", comp),
      out_stem = file.path(comp_out,
                            paste0(fig_date, "_", current_sex, "_", comp,
                                   "_Pooled_DOWN_ClusterDot"))
    )

    ##########################################################
    ### 3. SUPER-POOLED (all CTs merged into one Up / one Down set)
    ##########################################################
    cat("\n  Super-pooled (all cell types merged):\n")

    super_up <- unique(unlist(pool_up))
    super_dn <- unique(unlist(pool_dn))
    cat(sprintf("    UP: %d genes | DOWN: %d genes\n",
                length(super_up), length(super_dn)))

    for (nm in c("UP", "DOWN")) {
      genes <- if (nm == "UP") super_up else super_dn
      if (length(genes) < MIN_GENES) next

      # GO BP
      r_go <- run_go(genes, universe, "BP")
      sankey_dot_plot(
        r_go,
        title    = paste0("Super-pooled | ", nm, " | GO BP | ", comp,
                           " | ", tools::toTitleCase(current_sex)),
        out_stem = file.path(comp_out,
                              paste0(fig_date, "_", current_sex, "_", comp,
                                     "_SuperPooled_", nm, "_GO_BP_SankeyDot")),
        width    = 13
      )

      # KEGG
      r_kk <- run_kegg(genes, universe)
      sankey_dot_plot(
        r_kk,
        title    = paste0("Super-pooled | ", nm, " | KEGG | ", comp,
                           " | ", tools::toTitleCase(current_sex)),
        out_stem = file.path(comp_out,
                              paste0(fig_date, "_", current_sex, "_", comp,
                                     "_SuperPooled_", nm, "_KEGG_SankeyDot")),
        width    = 13
      )
    }

    # compareCluster dot for super-pooled UP vs DOWN together
    super_sets <- list()
    if (length(super_up) >= MIN_GENES) super_sets[["UP (all CTs)"]]   <- super_up
    if (length(super_dn) >= MIN_GENES) super_sets[["DOWN (all CTs)"]] <- super_dn

    cmp_super <- compare_dot(
      gene_lists = super_sets,
      universe   = universe,
      ont        = "BP",
      title      = paste0("GO BP — Super-pooled UP vs DOWN | ",
                           tools::toTitleCase(current_sex), " | ", comp),
      out_stem   = file.path(comp_out,
                              paste0(fig_date, "_", current_sex, "_", comp,
                                     "_SuperPooled_UPvsDOWN_compareCluster_Dot")),
      width = 10, height = 10
    )
    cluster_dot_plot(
      cmp_super,
      title    = paste0("GO BP — Super-pooled UP vs DOWN | ",
                         tools::toTitleCase(current_sex), " | ", comp),
      out_stem = file.path(comp_out,
                            paste0(fig_date, "_", current_sex, "_", comp,
                                   "_SuperPooled_UPvsDOWN_ClusterDot")),
      top_n = 8
    )

    cat(sprintf("\n  Done: %s | %s\n", current_sex, comp))
  }
}

cat("\n=== DONE ===\nOutputs:", out_dir, "\n")
cat("\nOutput structure:\n")
cat("  <sex>/<comp>/                      — pooled + super-pooled files\n")
cat("  <sex>/<comp>/<cell_type>/          — per-CT sankey+dot PDFs + PNGs\n")
cat("\nFile naming:\n")
cat("  *_UP_GO_BP_SankeyDot.{pdf,png,csv}\n")
cat("  *_DOWN_GO_BP_SankeyDot.{pdf,png,csv}\n")
cat("  *_UP_KEGG_SankeyDot.{pdf,png,csv}\n")
cat("  *_DOWN_KEGG_SankeyDot.{pdf,png,csv}\n")
cat("  *_Pooled_{UP,DOWN}_compareCluster_Dot.{pdf,png,csv}   ← classic style\n")
cat("  *_Pooled_{UP,DOWN}_ClusterDot.{pdf,png}               ← panel f/h style\n")
cat("  *_SuperPooled_{UP,DOWN}_{GO_BP,KEGG}_SankeyDot.{pdf,png,csv}\n")
cat("  *_SuperPooled_UPvsDOWN_compareCluster_Dot.{pdf,png,csv}\n")
cat("  *_SuperPooled_UPvsDOWN_ClusterDot.{pdf,png}\n")

invisible(NULL)
}

################################################################################
### RUN
################################################################################

status <- data.frame(section=character(0), result=character(0), minutes=numeric(0), stringsAsFactors=FALSE)
for (id in intersect(names(SECTION_INFO), RUN)) {
  if (!id %in% CRITICAL && !fig_ready()) {
    say("DESeq2 tables are incomplete: run step deseq2 first (add \"deseq2\" to RUN)"); break }
  banner(id)
  t0 <- Sys.time()
  res <- tryCatch({ SECTIONS[[id]](); "done" },
                  error=function(e) { say("  ERROR in step ", id, ": ", conditionMessage(e)); paste("FAILED:", conditionMessage(e)) })
  mins <- round(as.numeric(difftime(Sys.time(), t0, units="mins")), 1)
  invisible(gc())
  status <- rbind(status, data.frame(section=id, result=res, minutes=mins))
  say(sprintf("STEP %s %s in %.1f min", id, if (res == "done") "finished" else "FAILED", mins))
  if (res != "done" && id %in% CRITICAL) { say("Step ", id, " is required by everything after it: stopping"); break }
}

cat("\n", strrep("=", 90), "\n", sep="")
say("Figure 2 pipeline finished. Summary:")
for (i in seq_len(nrow(status)))
  say(sprintf("  %-7s %-6s %5.1f min  %s", status$section[i], ifelse(status$result[i] == "done", "ok", "FAIL"),
              status$minutes[i], SECTION_INFO[[status$section[i]]]))
if (any(status$result != "done")) say("Failed steps can be rerun alone, e.g. RUN <- c(\"2c\").")
say("Outputs: ", FIG_ROOT, " (one subfolder per panel). Log: ", LOG)
