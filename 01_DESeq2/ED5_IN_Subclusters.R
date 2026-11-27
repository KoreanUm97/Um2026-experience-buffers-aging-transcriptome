################################################################################
### ED5_IN_Subclusters.R
###
### Extended Data Fig. 5: activity-regulated genes in IN subtypes.
###
### Complete pipeline for annotated IN subtypes (in_subtype):
###   Sec 0 — IN subclustering + manual subtype annotation (from the integrated eNE object);
###           ED Fig. 5a: annotated UMAP, age / sex / condition UMAPs, subtype marker dot plot
###   Sec 1 — DESeq2 pseudobulk (EE vs HC per subtype × age × sex)
###   Sec 2 — IEG DotViolin
###   Sec 3 — ARG Discovery
###              exclusion: ExNs + all other IN subtypes
###   Sec 4 — Strip Volcano + ARG Transcriptome Similarity
###
### Male / Female / Combined throughout.
### Requires: the integrated eNE object loaded as df_harmony (for Sec 0), or in_obj already in the session,
###           or the subcluster object saved by Sec 0 (SUBCLUSTER_RDS).
###           The cell-type-level DESeq2 tables from Fig2_ED4_Pipeline.R (step deseq2), for the ExN exclusion in Sec 3.
### Outputs: <PATH_TO_OUTPUT>/Figures/ED5_IN_Subtypes/ and DESeq2_Results/IN_Subtypes/ (existing files are skipped).
################################################################################

suppressPackageStartupMessages({
  library(conflicted)
  library(Seurat); library(BPCells)
  library(DESeq2)
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggrepel); library(patchwork)
  conflicts_prefer(dplyr::select, .quiet=TRUE)
  conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::rename, .quiet=TRUE)
  conflicts_prefer(base::intersect,  .quiet=TRUE)
  conflicts_prefer(base::setdiff,    .quiet=TRUE)
  conflicts_prefer(dplyr::union)
  conflicts_prefer(GenomicRanges::union)
  conflicts_prefer(base::union)
  conflicts_prefer(SummarizedExperiment::as.data.frame)
  conflicts_prefer(BPCells::as.data.frame)
  conflicts_prefer(base::as.data.frame)
})

USE_BEESWARM <- tryCatch({ library(ggbeeswarm); TRUE },
                          error=function(e){ cat("ggbeeswarm unavailable\n"); FALSE })

################################################################################
### PARAMETERS
################################################################################

BASE    <- "PATH_TO_OUTPUT"                # same folder as BASE in Fig2_ED4_Pipeline.R
if (grepl("^PATH_TO_", BASE)) stop("Set BASE at the top of this script")
SUBCLUSTER_RDS <- file.path(BASE, "IN_Subcluster_Annotated.rds")   # written by Sec 0, reused on reruns

# Cell-type-level DESeq2 results (for ExN exclusion in Sec 3)
ext_deseq_dir <- paste0(file.path(BASE, "DESeq2_Results"), "/")

fig_date  <- "E5"   # file-name prefix
out_dir   <- paste0(file.path(BASE, "Figures", "ED5_IN_Subtypes"), "/")
pb_dir    <- paste0(file.path(BASE, "DESeq2_Results", "IN_Subtypes"), "/")  # pseudobulk DESeq2 CSVs
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(pb_dir,  recursive=TRUE, showWarnings=FALSE)
cat("Output :", out_dir, "\n")
cat("DESeq2 :", pb_dir,  "\n\n")

# ── thresholds ────────────────────────────────────────────────────────────────
PADJ_THR  <- 0.05
LFC_THR   <- 0
TOP_N     <- 30
TOP_PER_CT <- 15
MIN_PB    <- 3      # min pseudobulk samples per condition to run DESeq2

# ── subtype constants ─────────────────────────────────────────────────────────
SUBTYPE_LEVELS <- c("PV_IN","SST_IN","VIP_IN","LAMP5_IN")
SUBTYPE_COLS   <- c(PV_IN="#E41A1C", SST_IN="#FF7F00",
                    VIP_IN="#4DAF4A", LAMP5_IN="#984EA3")
EXCIT_CTS      <- c("DG ExNs","CA1 ExNs","CA2 ExNs","CA3 ExNs","Sub ExNs")

# ── comparisons ───────────────────────────────────────────────────────────────
COMPS <- list(
  list(name="YEEYHC", age="young",     label="Young EE vs HC"),
  list(name="MEEMHC", age="middleage", label="Middle EE vs HC"),
  list(name="OEEOHC", age="old",       label="Old EE vs HC")
)
AGE_LABELS <- c(YEEYHC="Young EE", MEEMHC="Middle EE", OEEOHC="Old EE")

# ── IEG gene list ─────────────────────────────────────────────────────
IEG_GENES <- c(
  "Fos","Fosb","Fosl1","Fosl2","Jun","Junb","Jund","Arc","Npas4",
  "Egr1","Egr2","Egr3","Egr4","Nr4a1","Nr4a2","Nr4a3",
  "Bdnf","Homer1","Nptx2","Dusp1","Dusp5","Gadd45b","Ier2","Per1"
)

# ── colors ────────────────────────────────────────────────────────────────────
COND_COLS  <- c(HC="#1565C0", EE="#FF1744")
VOL_COLS   <- c(Up="#D6604D", Down="#4393C3", NS="grey82")
MAX_W      <- 0.35
LFC_CAP    <- 5

# ── helpers ───────────────────────────────────────────────────────────────────
SKIP_EXISTING <- TRUE   # TRUE = do not overwrite figures that already exist
save_both <- function(p, stem, w, h, dpi=300) {
  if (SKIP_EXISTING && file.exists(paste0(stem, ".pdf")) && file.exists(paste0(stem, ".png"))) {
    cat("  Skip (exists):", basename(stem), "\n"); return(invisible(NULL)) }
  ggsave(paste0(stem, ".pdf"), p, width=w, height=h,
         device="pdf", limitsize=FALSE)
  ggsave(paste0(stem, ".png"), p, width=w, height=h,
         device="png", dpi=dpi, bg="white", limitsize=FALSE)
  cat("  Saved:", basename(stem), "\n")
}

skip_csv <- function(path) file.exists(path)

skip_fig <- function(stem) {
  file.exists(paste0(stem, ".pdf")) && file.exists(paste0(stem, ".png"))
}

vln_theme <- function(base_size=9) {
  theme_classic(base_size=base_size) +
  theme(
    strip.text.y.left  = element_text(angle=0, hjust=1, size=7, face="italic",
                                       color="grey10", margin=margin(r=4,l=2)),
    strip.text.x       = element_text(size=9, face="bold", margin=margin(b=3)),
    strip.background   = element_blank(), strip.placement="outside",
    axis.text.x        = element_text(angle=45, hjust=1, vjust=1,
                                       size=8, color="grey15"),
    axis.ticks.x       = element_line(linewidth=0.3, color="grey40"),
    axis.line.x        = element_line(linewidth=0.35, color="grey40"),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.line.y        = element_blank(),
    axis.title         = element_blank(),
    panel.spacing.y    = unit(1, "pt"), panel.spacing.x = unit(4, "pt"),
    panel.grid         = element_blank(),
    panel.background   = element_rect(fill="grey98", color=NA),
    legend.position    = "none",
    plot.title         = element_text(size=10, face="bold", hjust=0.5),
    plot.subtitle      = element_text(size=8, color="grey40", hjust=0.5,
                                       margin=margin(b=6)),
    plot.margin        = margin(6,8,4,4)
  )
}

################################################################################
### SEC 0: IN SUBCLUSTERING (skipped if in_obj is in the session or SUBCLUSTER_RDS exists)
###   Subset IN from the integrated eNE object, Harmony by sample, cluster at
###   resolution 0.55 and apply the manual cluster -> subtype annotation.
################################################################################

build_in_subclusters <- function() {
  library(harmony)
  if (!exists("df_harmony", envir=.GlobalEnv)) stop("Load the integrated eNE object as df_harmony first")
  df_harmony <- get("df_harmony", envir=.GlobalEnv)

  # ── 1. Subset ──────────────────────────────────────────────────────────────
  in_obj <- subset(df_harmony, cells=colnames(df_harmony)[df_harmony$celltype=="IN"])
  in_obj <- JoinLayers(in_obj, assay="RNA")
  for (lyr in intersect(c("counts","data"), Layers(in_obj[["RNA"]]))) {   # deposited object has counts only
    mat <- GetAssayData(in_obj, assay="RNA", layer=lyr)
    if (!inherits(mat,"dgCMatrix")) {
      mat <- as(mat,"dgCMatrix")
      in_obj[["RNA"]] <- SetAssayData(in_obj[["RNA"]], layer=lyr, new.data=mat)
    }
  }
  cat("IN cells:", ncol(in_obj), "\n")

  # ── 2. Process ─────────────────────────────────────────────────────────────
  in_obj <- NormalizeData(in_obj, verbose=FALSE)
  in_obj <- FindVariableFeatures(in_obj, nfeatures=2000, verbose=FALSE)
  in_obj <- ScaleData(in_obj, verbose=FALSE)
  in_obj <- RunPCA(in_obj, npcs=30, verbose=FALSE)
  in_obj <- RunHarmony(in_obj, group.by.vars="orig.ident",
                        reduction="pca", reduction.save="harmony_in",
                        verbose=FALSE)
  in_obj <- FindNeighbors(in_obj, reduction="harmony_in",
                           dims=1:20, verbose=FALSE)
  in_obj <- RunUMAP(in_obj, reduction="harmony_in", dims=1:20,
                     reduction.name="umap_in", verbose=FALSE)

  # ── 3. Cluster at res=0.55 ─────────────────────────────────────────────────
  in_obj <- FindClusters(in_obj, resolution=0.55, verbose=FALSE,
                          cluster.name="in_cluster_res0.55")
  in_obj@meta.data$in_cluster_res0.55 <- paste0(
    "C", as.character(in_obj@meta.data$in_cluster_res0.55))
  cat("Clusters found:",
      paste(sort(unique(in_obj@meta.data$in_cluster_res0.55)), collapse=", "), "\n")

  # ── 4. Manual annotation ───────────────────────────────────────────────────
  CLUSTER_MAP <- c(
    "C0"="PV_IN",    "C1"="SST_IN",   "C2"="LAMP5_IN", "C3"="VIP_IN",
    "C4"="PV_IN",    "C5"="SST_IN",   "C6"="VIP_IN",   "C7"="SST_IN",
    "C8"="SST_IN",   "C9"="VIP_IN",   "C10"="VIP_IN",  "C11"="PV_IN",
    "C12"="LAMP5_IN","C13"="LAMP5_IN","C14"="LAMP5_IN","C15"="Ambiguous",
    "C16"="PV_IN",   "C17"="SST_IN",  "C18"="PV_IN",   "C19"="SST_IN"
  )

  # Check all clusters are covered
  found <- unique(in_obj@meta.data$in_cluster_res0.55)
  missing <- setdiff(found, names(CLUSTER_MAP))
  if (length(missing)>0) {
    cat("WARNING: unmapped clusters:", paste(missing, collapse=", "), "\n")
    CLUSTER_MAP[missing] <- "Ambiguous"
  }

  in_obj@meta.data$in_subtype <- factor(
    CLUSTER_MAP[as.character(in_obj@meta.data$in_cluster_res0.55)],
    levels=c("PV_IN","SST_IN","VIP_IN","LAMP5_IN","Ambiguous"))

  # Use the deposited labels if present: the integrated eNE object in GEO carries in_subtype
  # (the annotation used in the paper); cluster numbering can change between runs.
  if ("in_subtype" %in% colnames(df_harmony@meta.data)) {
    dep <- as.character(df_harmony@meta.data[colnames(in_obj), "in_subtype"])
    cat("\nAgreement of this run's annotation with the deposited in_subtype labels:\n")
    print(table(this_run=as.character(in_obj$in_subtype), deposited=dep, useNA="ifany"))
    in_obj@meta.data$in_subtype <- factor(dep, levels=c("PV_IN","SST_IN","VIP_IN","LAMP5_IN","Ambiguous"))
    cat("Using the deposited in_subtype labels\n")
  }

  cat("\nSubtype counts:\n")
  print(sort(table(in_obj@meta.data$in_subtype), decreasing=TRUE))

  # ── 5. Save ────────────────────────────────────────────────────────────────
  saveRDS(in_obj, SUBCLUSTER_RDS)
  cat("Saved:", SUBCLUSTER_RDS, "\n")
  in_obj
}

################################################################################
### LOAD in_obj
################################################################################

if (!exists("in_obj")) {
  if (file.exists(SUBCLUSTER_RDS)) {
    cat("Loading:", SUBCLUSTER_RDS, "\n")
    in_obj <- readRDS(SUBCLUSTER_RDS)
  } else {
    in_obj <- build_in_subclusters()
  }
}

################################################################################
### SEC 0b: ED Fig. 5a — annotated IN UMAP (+ age / sex / condition insets) and
###          subtype marker dot plot (Ambiguous cells still included here)
################################################################################

plot_in_subclusters <- function(in_obj) {
  sub_out <- file.path(out_dir, "Sec0_Subclustering")
  dir.create(sub_out, recursive=TRUE, showWarnings=FALSE)
  COL_NM <- "in_cluster_res0.55"

  SUBTYPE_COLS_A <- c(
    "PV_IN"     = "#E41A1C",
    "SST_IN"    = "#FF7F00",
    "VIP_IN"    = "#4DAF4A",
    "LAMP5_IN"  = "#984EA3",
    "Ambiguous" = "#999999"
  )
  SUBTYPE_ORDER <- c("PV_IN","SST_IN","VIP_IN","LAMP5_IN","Ambiguous")
  in_obj@meta.data$in_subtype <- factor(as.character(in_obj@meta.data$in_subtype),
    levels=intersect(SUBTYPE_ORDER, unique(as.character(in_obj@meta.data$in_subtype))))
  Idents(in_obj) <- in_obj@meta.data$in_subtype

  make_umap <- function(color_by, cols=NULL, title, out_file, label=TRUE) {
    Idents(in_obj) <- in_obj@meta.data[[color_by]]
    p <- DimPlot(in_obj,
                 reduction = "umap_in",
                 label     = FALSE,
                 pt.size   = 0.4) +
      labs(title=title) +
      theme(plot.title  = element_text(face="bold",size=13,hjust=0.5),
            legend.text = element_text(size=9))
    if (!is.null(cols)) p <- p + scale_color_manual(values=cols)
    if (label) {
      umap_coords <- Embeddings(in_obj, "umap_in")
      label_df <- data.frame(
        x       = umap_coords[,1],
        y       = umap_coords[,2],
        cluster = as.character(Idents(in_obj)),
        stringsAsFactors=FALSE
      ) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(x=median(x), y=median(y), .groups="drop")
      p <- p + geom_text(
        data        = label_df,
        mapping     = aes(x=x, y=y, label=cluster),
        size        = 4,
        fontface    = "bold",
        color       = "black",
        inherit.aes = FALSE
      )
    }
    ggsave(out_file, p, width=9, height=8, device="pdf")
    cat("Saved:", basename(out_file), "\n")
  }

  # Main annotated UMAP
  make_umap("in_subtype",
            cols     = SUBTYPE_COLS_A,
            title    = "IN subclusters -- manually annotated",
            out_file = file.path(sub_out, paste0(fig_date,"_IN_UMAP_Annotated.pdf")))

  # Split by subtype (one panel per subtype)
  out_split <- file.path(sub_out, paste0(fig_date,"_IN_UMAP_Split_bySubtype.pdf"))
  Idents(in_obj) <- in_obj@meta.data$in_subtype
  p_split <- DimPlot(in_obj,
                      reduction = "umap_in",
                      split.by  = "in_subtype",
                      ncol      = 3,
                      pt.size   = 0.3,
                      label     = FALSE,
                      cols      = SUBTYPE_COLS_A) +
    labs(title="IN subclusters split by subtype") +
    theme(plot.title=element_text(face="bold",size=12,hjust=0.5))
  n_types <- length(unique(in_obj@meta.data$in_subtype))
  ggsave(out_split, p_split, width=9, height=ceiling(n_types/3)*3.5, device="pdf")
  cat("Saved:", basename(out_split), "\n")

  # UMAP showing the res-0.55 cluster IDs (reference)
  if (COL_NM %in% names(in_obj@meta.data))
    make_umap(COL_NM,
              title    = "IN subclusters (res=0.55) -- original IDs",
              out_file = file.path(sub_out, paste0(fig_date,"_IN_UMAP_OriginalClusters.pdf")))

  # Context UMAPs (insets)
  for (col in c("age","sex","experimental.condition")) {
    if (col %in% names(in_obj@meta.data)) {
      make_umap(col,
                title    = paste0("IN cells -- ", col),
                out_file = file.path(sub_out, paste0(fig_date,"_IN_UMAP_",col,".pdf")),
                label    = FALSE)
    }
  }

  # Dot plot of canonical markers per annotated subtype
  IN_MARKERS_ORDERED <- c(
    "Pvalb","Kcnc1","Kcnmb4","Pthlh",          # PV
    "Sst","Npy","Cort","Chodl","Elfn1",          # SST
    "Vip","Calb2","Cck","Chrm2",                 # VIP
    "Lamp5","Nos1","Ndnf","Lhx6",               # LAMP5
    "Gad1","Gad2","Slc32a1"                      # pan-IN
  )
  markers_present <- intersect(IN_MARKERS_ORDERED, rownames(in_obj))
  out_dot <- file.path(sub_out, paste0(fig_date,"_IN_DotPlot_Annotated.pdf"))
  p_dot <- DotPlot(in_obj,
                    features = markers_present,
                    assay    = "RNA",
                    cols     = c("grey90","darkblue")) +
    scale_size(name="% expressing") +
    theme_classic(base_size=10) +
    theme(
      axis.text.x = element_text(angle=45,hjust=1,size=9,face="italic"),
      axis.text.y = element_text(size=10),
      plot.title  = element_text(face="bold",size=12,hjust=0.5)
    ) +
    labs(title="IN subtypes -- canonical marker expression", x=NULL, y=NULL)
  ggsave(out_dot, p_dot, width=max(10, length(markers_present)*0.45+3), height=5, device="pdf")
  cat("Saved:", basename(out_dot), "\n")

  # Annotation summary
  summ <- in_obj@meta.data %>%
    dplyr::group_by(in_subtype) %>%
    dplyr::summarise(
      n_cells    = dplyr::n(),
      n_clusters = if (COL_NM %in% names(in_obj@meta.data)) length(unique(.data[[COL_NM]])) else NA_integer_,
      clusters   = if (COL_NM %in% names(in_obj@meta.data)) paste(sort(unique(.data[[COL_NM]])), collapse=", ") else NA_character_,
      .groups="drop") %>%
    dplyr::arrange(dplyr::desc(n_cells))
  print(summ)
  utils::write.csv(summ, file.path(sub_out, paste0(fig_date,"_IN_Annotation_Summary.csv")), row.names=FALSE)
  invisible(NULL)
}
plot_in_subclusters(in_obj)

DefaultAssay(in_obj) <- "RNA"
in_obj <- JoinLayers(in_obj, assay="RNA")
for (lyr in c("counts","data")) {
  mat <- GetAssayData(in_obj, assay="RNA", layer=lyr)
  if (!inherits(mat,"dgCMatrix")) {
    mat <- as(mat,"dgCMatrix")
    in_obj[["RNA"]] <- SetAssayData(in_obj[["RNA"]], layer=lyr, new.data=mat)
  }
}
rm(mat); gc()

in_obj@meta.data$in_subtype <- factor(
  in_obj@meta.data$in_subtype, levels=SUBTYPE_LEVELS)

in_obj@meta.data$sex_norm <- tolower(in_obj@meta.data$sex)
in_obj@meta.data$experimental.condition <- factor(
  in_obj@meta.data$experimental.condition, levels=c("HC","EE"))

# ── recover age from sex×age condition column if needed ───────────────────────
# 'condition' encodes sex×age (e.g. "male_young", "female_middleage").
# When in_obj$age has NAs (common after subclustering), recover from there.
n_age_na <- sum(is.na(in_obj@meta.data$age))
if (n_age_na > 0) {
  cat(sprintf("NOTE: age has %d NAs — recovering from condition column\n", n_age_na))
  age_from_cond <- sub("^(male|female)_", "",
                        tolower(as.character(in_obj@meta.data$condition)))
  age_from_cond <- dplyr::recode(age_from_cond,
                                   "young"     = "Young",
                                   "middleage" = "Middleage",
                                   "old"       = "Old",
                                   .default    = NA_character_)
  na_idx <- is.na(in_obj@meta.data$age)
  in_obj@meta.data$age[na_idx] <- age_from_cond[na_idx]
  in_obj@meta.data$age <- factor(in_obj@meta.data$age,
                                   levels=c("Young","Middleage","Old"))
  cat("Recovered age:\n"); print(table(in_obj@meta.data$age, useNA="always")); cat("\n")
}

subtypes_present <- intersect(SUBTYPE_LEVELS,
                               levels(droplevels(in_obj@meta.data$in_subtype)))
genes_present    <- intersect(IEG_GENES, rownames(in_obj))
cat(sprintf("Cells: %d | Subtypes: %s\n", ncol(in_obj),
            paste(subtypes_present, collapse=", ")))
cat(sprintf("IEGs found: %d/%d\n\n", length(genes_present), length(IEG_GENES)))


################################################################################
### SEC 1: DESeq2 PSEUDOBULK PER SUBTYPE × AGE × SEX
################################################################################

cat("╔══════════════════════════════════════════════════╗\n")
cat("  Sec 1: DESeq2 pseudobulk\n")
cat("╚══════════════════════════════════════════════════╝\n\n")

counts_mat <- GetAssayData(in_obj, assay="RNA", layer="counts")
meta_full  <- in_obj@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::mutate(age_lower = tolower(as.character(age))) %>%
  dplyr::select(barcode, in_subtype, age_lower, sex_norm,
                 cond=experimental.condition, orig.ident)



run_deseq2 <- function(pb_mat, col_data, design_formula) {
  # require MIN_PB samples per condition
  cond_tab <- table(col_data$cond)
  if (any(cond_tab < MIN_PB)) return(NULL)

  dds <- tryCatch(
    DESeqDataSetFromMatrix(countData = round(pb_mat),
                            colData   = col_data,
                            design    = design_formula),
    error=function(e){ message("  DESeq2 setup error: ", e$message); NULL })
  if (is.null(dds)) return(NULL)

  # filter low-count genes
  keep <- rowSums(counts(dds) >= 5) >= MIN_PB
  dds  <- dds[keep, ]
  if (nrow(dds) < 100) return(NULL)

  dds <- tryCatch(DESeq(dds, quiet=TRUE),
                   error=function(e){ message("  DESeq error: ", e$message); NULL })
  if (is.null(dds)) return(NULL)

  res <- results(dds, contrast=c("cond","EE","HC"), independentFiltering=TRUE)
  as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::select(gene, log2FoldChange, padj, baseMean, lfcSE) %>%
    dplyr::filter(!is.na(log2FoldChange))
}

make_pseudobulk <- function(meta_sub) {
  uniq_oids <- unique(meta_sub$orig.ident)
  if (length(uniq_oids) < (MIN_PB * 2)) return(NULL)

  n_genes <- nrow(counts_mat)
  pb_mat  <- vapply(uniq_oids, function(oid) {
    bc <- intersect(meta_sub$barcode[meta_sub$orig.ident == oid],
                    colnames(counts_mat))
    if (length(bc) == 0) return(rep(0.0, n_genes))
    as.numeric(Matrix::rowSums(counts_mat[, bc, drop=FALSE]))
  }, numeric(n_genes))
  # vapply always returns a matrix (n_genes × n_oids)
  rownames(pb_mat) <- rownames(counts_mat)
  colnames(pb_mat) <- uniq_oids

  col_data <- data.frame(
    sample = uniq_oids,
    cond   = factor(
      sapply(uniq_oids, function(oid)
        as.character(meta_sub$cond[meta_sub$orig.ident == oid][1])),
      levels=c("HC","EE")),
    row.names = uniq_oids,
    stringsAsFactors=FALSE
  )
  list(mat=pb_mat, coldata=col_data)
}

all_deseq <- list()

for (sx in c("male","female")) {
  cat(sprintf("── %s ──\n", sx))
  for (comp in COMPS) {
    cat(sprintf("  %s\n", comp$label))
    for (st in subtypes_present) {
      csv_path <- file.path(pb_dir,
                             sprintf("%s_%s_%s.csv", sx, st, comp$name))
      if (skip_csv(csv_path)) {
        cat(sprintf("    %s — load CSV\n", st))
        res <- read.csv(csv_path, stringsAsFactors=FALSE)
      } else {
        meta_sub <- meta_full %>%
          dplyr::filter(sex_norm==sx,
                         age_lower==comp$age,
                         in_subtype==st)
        cat(sprintf("    %s: %d cells", st, nrow(meta_sub)))
        pb <- make_pseudobulk(meta_sub)
        if (is.null(pb)) { cat(" — skip (too few samples)\n"); next }
        res <- run_deseq2(pb$mat, pb$coldata, ~ cond)
        if (is.null(res)) { cat(" — DESeq2 failed\n"); next }
        cat(sprintf(" — %d genes | %d sig\n",
                    nrow(res),
                    sum(!is.na(res$padj) & res$padj < PADJ_THR &
                          abs(res$log2FoldChange) >= LFC_THR)))
        write.csv(res, csv_path, row.names=FALSE)
      }
      res$celltype   <- st
      res$sex        <- sx
      res$comparison <- comp$name
      all_deseq[[paste(sx, st, comp$name)]] <- res
    }
  }
}

# Combined: all sexes, sex as covariate
cat("── combined ──\n")
for (comp in COMPS) {
  cat(sprintf("  %s\n", comp$label))
  for (st in subtypes_present) {
    csv_path <- file.path(pb_dir, sprintf("combined_%s_%s.csv", st, comp$name))
    if (skip_csv(csv_path)) {
      cat(sprintf("    %s — load CSV\n", st))
      res <- read.csv(csv_path, stringsAsFactors=FALSE)
    } else {
      meta_sub <- meta_full %>%
        dplyr::filter(age_lower==comp$age, in_subtype==st)
      cat(sprintf("    %s: %d cells", st, nrow(meta_sub)))
      pb <- make_pseudobulk(meta_sub)
      if (is.null(pb)) { cat(" — skip\n"); next }
      # add sex to colData
      pb$coldata$sex <- factor(
        sapply(rownames(pb$coldata), function(oid)
          meta_sub$sex_norm[meta_sub$orig.ident == oid][1]))
      if (length(levels(droplevels(pb$coldata$sex))) > 1) {
        design_f <- ~ sex + cond
      } else {
        design_f <- ~ cond
      }
      res <- run_deseq2(pb$mat, pb$coldata, design_f)
      if (is.null(res)) { cat(" — DESeq2 failed\n"); next }
      cat(sprintf(" — %d genes | %d sig\n",
                  nrow(res),
                  sum(!is.na(res$padj) & res$padj < PADJ_THR &
                        abs(res$log2FoldChange) >= LFC_THR)))
      write.csv(res, csv_path, row.names=FALSE)
    }
    res$celltype   <- st
    res$sex        <- "combined"
    res$comparison <- comp$name
    all_deseq[[paste("combined", st, comp$name)]] <- res
  }
}

full_df <- dplyr::bind_rows(all_deseq) %>%
  dplyr::mutate(sig = !is.na(padj) & padj < PADJ_THR &
                        abs(log2FoldChange) >= LFC_THR)

cat(sprintf("\nSec 1 complete: %d result rows | %d sig\n\n",
            nrow(full_df), sum(full_df$sig, na.rm=TRUE)))

################################################################################
### SEC 2: IEG DOTVIOLIN
################################################################################

cat("╔══════════════════════════════════════════════════╗\n")
cat("  Sec 2: IEG DotViolin\n")
cat("╚══════════════════════════════════════════════════╝\n\n")

# % expressed: compute from in_obj for EE cells per subtype × age
ig_counts <- as(
  GetAssayData(in_obj, assay="RNA", layer="counts")[genes_present, ],
  "dgCMatrix")

compute_pct <- function(meta_rows) {
  bc <- intersect(meta_rows$barcode, colnames(ig_counts))
  if (length(bc) == 0) return(NULL)
  mat <- ig_counts[, bc, drop=FALSE]
  data.frame(gene    = genes_present,
             pct_expr = rowMeans(mat > 0) * 100,
             stringsAsFactors=FALSE)
}

make_ieg_dotplot <- function(dot_df, label) {
  if (is.null(dot_df) || nrow(dot_df)==0) return(invisible())
  dot_df$celltype <- factor(dot_df$celltype, levels=subtypes_present)   # PV, SST, VIP, LAMP5 left to right
  p <- ggplot(dot_df, aes(x=celltype, y=gene)) +
    geom_point(data=dot_df %>% dplyr::filter(!sig_dot),
               aes(size=pct_expr), color="grey88", alpha=0.7) +
    geom_point(data=dot_df %>% dplyr::filter(sig_dot),
               aes(color=lfc_cap, size=pct_expr), alpha=0.92) +
    scale_color_gradient2(low="#1E88E5", mid="grey97", high="#E53935",
                           midpoint=0, limits=c(-3,3), oob=scales::squish,
                           name="log2FC\n(EE vs HC)") +
    scale_size_continuous(name="% EE cells\nexpressing",
                           range=c(0.3,10), limits=c(0,100),
                           breaks=c(5,10,25,50,75)) +
    facet_wrap(~age_lbl, nrow=1) +
    scale_y_discrete(limits=rev(genes_present)) +
    theme_classic(base_size=9) +
    theme(plot.title   = element_text(face="bold", size=11, hjust=0.5),
          plot.subtitle= element_text(size=7.5, color="grey40", hjust=0.5),
          axis.text.x  = element_text(angle=45, hjust=1, size=8),
          axis.text.y  = element_text(size=8, face="italic"),
          axis.title   = element_blank(),
          strip.text   = element_text(face="bold", size=9),
          legend.position="right",
          panel.grid.major=element_line(color="grey92", linewidth=0.3),
          panel.spacing=unit(0.5,"lines")) +
    guides(color=guide_colorbar(barwidth=0.8, barheight=5,
                                 title.position="top"),
           size=guide_legend(title.position="top",
                              override.aes=list(color="grey50"))) +
    labs(title    = paste0("IEG/ARG expression — IN subtypes — ", label),
         subtitle = "Color=log2FC (DESeq2 EE vs HC) | Size=% EE expressing | Grey=NS")
  save_both(p, file.path(out_dir, paste0(fig_date, "_Sec2_IEG_Dotplot_", label)),
            w=max(7, length(subtypes_present)*0.9+3),
            h=max(6, length(genes_present)*0.35+3))
}

for (sx_label in c("Male","Female","Combined")) {
  sx_key <- tolower(sx_label)
  stem   <- file.path(out_dir, paste0(fig_date, "_Sec2_IEG_Dotplot_", sx_label))
  if (skip_fig(stem)) { cat("Skip:", basename(stem), "\n"); next }

  pct_rows <- list()
  for (comp in COMPS) {
    for (st in subtypes_present) {
      m <- meta_full %>%
        dplyr::filter(age_lower==comp$age,
                       in_subtype==st,
                       cond=="EE",
                       if (sx_key=="combined") TRUE else sex_norm==sx_key)
      pct <- compute_pct(m)
      if (!is.null(pct))
        pct_rows[[paste(comp$name,st)]] <- pct %>%
          dplyr::mutate(celltype=st, comparison=comp$name)
    }
  }
  pct_df <- dplyr::bind_rows(pct_rows)

  lfc_df <- full_df %>%
    dplyr::filter(sex==(if(sx_key=="combined") "combined" else sx_key),
                   gene %in% genes_present) %>%
    dplyr::mutate(sig_dot  = sig,
                  lfc_cap  = pmax(pmin(log2FoldChange,3),-3),
                  age_lbl  = factor(AGE_LABELS[comparison],
                                     levels=c("Young EE","Middle EE","Old EE")))

  dot_df <- dplyr::inner_join(lfc_df, pct_df,
                               by=c("celltype","comparison","gene")) %>%
    dplyr::mutate(gene=factor(gene, levels=genes_present))

  make_ieg_dotplot(dot_df, sx_label)
}

# ── Stacked violin: genes × condition × age, per subtype ─────────────────────
cat("\n  [Stacked violin]\n")

expr_df <- FetchData(in_obj,
                      vars=c("experimental.condition","sex_norm","age",
                             "in_subtype", genes_present),
                      layer="data") %>%
  dplyr::mutate(
    experimental.condition=factor(experimental.condition, levels=c("HC","EE")),
    in_subtype=factor(in_subtype, levels=subtypes_present),
    age=factor(age, levels=c("Young","Middleage","Old")))

expr_long <- expr_df %>%
  tidyr::pivot_longer(cols=all_of(genes_present),
                       names_to="gene", values_to="expression") %>%
  dplyr::mutate(gene=factor(gene, levels=genes_present))

make_stacked_violin <- function(df, label) {
  if (nrow(df)==0) return(invisible())
  # across ages, HC vs EE, facet gene × age, separate panel per subtype
  for (st in subtypes_present) {
    sub <- df %>% dplyr::filter(in_subtype==st)
    if (nrow(sub)<200) next
    stem <- file.path(out_dir,
                       paste0(fig_date,"_Sec2_IEG_Violin_",st,"_",label))
    if (skip_fig(stem)) { cat("  Skip:", basename(stem), "\n"); next }
    p <- ggplot(sub, aes(x=experimental.condition, y=expression,
                          fill=experimental.condition)) +
      geom_violin(scale="width", trim=TRUE, linewidth=0.12, color="grey30") +
      facet_grid(gene ~ age, scales="free_y", switch="y") +
      scale_fill_manual(values=COND_COLS) +
      scale_x_discrete(expand=expansion(add=0.4)) +
      scale_y_continuous(expand=expansion(mult=c(0.05,0.15))) +
      vln_theme() +
      labs(title    = paste0(st, " — IEG/ARG HC vs EE | ", label),
           subtitle = "Columns=Young|Middle|Old | Blue=HC | Red=EE")
    save_both(p, stem,
              w=7, h=max(8, length(genes_present)*0.45+3))
  }
}

for (sx_label in c("Male","Female","Combined")) {
  sub <- if (sx_label=="Combined") expr_long else
    expr_long %>% dplyr::filter(sex_norm==tolower(sx_label))
  make_stacked_violin(sub, sx_label)
}

cat("\nSec 2 complete.\n\n")

################################################################################
################################################################################
### SEC 3: ARG DISCOVERY — SUBTYPE-SPECIFIC PER AGE
### Exclusion: sig in any ExN OR sig in any other IN subtype
### Runs per age comparison (YEEYHC, MEEMHC, OEEOHC) independently.
### Per sex: Male, Female, Combined
### Per age: separate plots + union plots across ages
################################################################################

cat("╔══════════════════════════════════════════════════╗\n")
cat("  Sec 3: ARG Discovery (subtype-specific, per age)\n")
cat("╚══════════════════════════════════════════════════╝\n\n")

# Load existing ExN DESeq2 results
ext_df <- list()
for (sx in c("male","female")) {
  for (comp_name in c("YEEYHC","MEEMHC","OEEOHC")) {
    f <- file.path(ext_deseq_dir, paste0(sx,"_ALL_",comp_name,".csv"))
    if (!file.exists(f)) next
    d <- read.csv(f, stringsAsFactors=FALSE) %>%
      dplyr::mutate(sex=sx, comparison=comp_name,
                    sig=!is.na(padj)&padj<PADJ_THR&abs(log2FoldChange)>=LFC_THR)
    ext_df[[paste(sx,comp_name)]] <- d
  }
}
ext_full <- dplyr::bind_rows(ext_df)
cat(sprintf("External DESeq2 loaded: %d rows\n\n", nrow(ext_full)))

# all_specific keyed as "{sx}__{comp}__{subtype}"
all_specific <- list()

for (sx_label in c("Male","Female","Combined")) {
  sx_key <- tolower(sx_label)
  cat(sprintf("── %s ──\n", sx_label))

  # ── Discover ARGs per age × subtype ────────────────────────────────────────
  for (comp in COMPS) {
    cat(sprintf("  %s\n", comp$label))
    if (sx_key=="combined") {
      comp_df   <- full_df %>% dplyr::filter(sex=="combined", comparison==comp$name)
      excit_sig <- ext_full %>%
        dplyr::filter(celltype %in% EXCIT_CTS, sig, comparison==comp$name) %>%
        dplyr::pull(gene) %>% unique()
    } else {
      comp_df   <- full_df %>% dplyr::filter(sex==sx_key, comparison==comp$name)
      excit_sig <- ext_full %>%
        dplyr::filter(sex==sx_key, celltype %in% EXCIT_CTS, sig,
                       comparison==comp$name) %>%
        dplyr::pull(gene) %>% unique()
    }
    for (st in subtypes_present) {
      st_sig <- comp_df %>% dplyr::filter(celltype==st, sig) %>% dplyr::pull(gene)
      if (length(st_sig)==0) next
      other_in_sig <- comp_df %>%
        dplyr::filter(celltype!=st, sig) %>% dplyr::pull(gene) %>% unique()
      specific_genes <- setdiff(st_sig, union(excit_sig, other_in_sig))
      if (length(specific_genes)==0) {
        cat(sprintf("    %s: 0 specific ARGs\n", st)); next
      }
      summ <- comp_df %>%
        dplyr::filter(celltype==st, gene %in% specific_genes) %>%
        dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
        dplyr::mutate(sex=sx_key, subtype=st,
                       comparison=comp$name, age_label=comp$label)
      key <- paste(sx_key, comp$name, st, sep="__")
      all_specific[[key]] <- summ
      cat(sprintf("    %s: %d specific ARGs\n", st, nrow(summ)))
      write.csv(summ,
                file.path(out_dir,
                           paste0(fig_date,"_Sec3_SpecificARGs_",
                                  sx_label,"_",comp$name,"_",st,".csv")),
                row.names=FALSE)
    }
  }

  # ── top candidates: per age and union ──────────────────────────────────────
  top_cands_age <- setNames(lapply(COMPS, function(comp) {
    age_keys <- names(all_specific)[
      startsWith(names(all_specific), paste0(sx_key,"__",comp$name,"__"))]
    age_dat  <- dplyr::bind_rows(all_specific[age_keys])
    if (nrow(age_dat)==0) return(character(0))
    age_dat %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(mean_abs_lfc=mean(abs(log2FoldChange)), .groups="drop") %>%
      dplyr::arrange(dplyr::desc(mean_abs_lfc)) %>%
      dplyr::slice_head(n=TOP_N) %>% dplyr::pull(gene)
  }), sapply(COMPS, `[[`, "name"))

  top_cands_union <- unique(unlist(top_cands_age))
  cat(sprintf(
    "\n  Top candidates — Young: %d | Middle: %d | Old: %d | Union: %d\n\n",
    length(top_cands_age[["YEEYHC"]]), length(top_cands_age[["MEEMHC"]]),
    length(top_cands_age[["OEEOHC"]]), length(top_cands_union)))
  if (length(top_cands_union)==0) next

  # ── Plot 1: Scatter — per age ──────────────────────────────────────────────
  cat("  [Plot 1: Scatter — per age]\n")
  for (comp in COMPS) {
    stem <- file.path(out_dir,
                       paste0(fig_date,"_Sec3_Scatter_",comp$name,"_",sx_label))
    if (skip_fig(stem)) { cat("  Skip:", basename(stem), "\n"); next }
    if (sx_key=="combined") {
      comp_df <- full_df %>% dplyr::filter(sex=="combined", comparison==comp$name)
      excit_g <- ext_full %>%
        dplyr::filter(celltype %in% EXCIT_CTS, sig, comparison==comp$name) %>%
        dplyr::pull(gene) %>% unique()
    } else {
      comp_df <- full_df %>% dplyr::filter(sex==sx_key, comparison==comp$name)
      excit_g <- ext_full %>%
        dplyr::filter(sex==sx_key, celltype %in% EXCIT_CTS, sig,
                       comparison==comp$name) %>%
        dplyr::pull(gene) %>% unique()
    }
    p_list <- lapply(subtypes_present, function(st) {
      st_lfc <- comp_df %>%
        dplyr::filter(celltype==st, !is.na(log2FoldChange)) %>%
        dplyr::select(gene, st_lfc=log2FoldChange, st_sig=sig)
      ex_lfc <- comp_df %>%
        dplyr::filter(celltype!=st, !is.na(log2FoldChange)) %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(ex_lfc=mean(log2FoldChange,na.rm=TRUE),
                         ex_sig=any(sig), .groups="drop")
      sc <- dplyr::inner_join(st_lfc, ex_lfc, by="gene") %>%
        dplyr::mutate(category=dplyr::case_when(
          st_sig & !ex_sig & !gene %in% excit_g ~ "Specific",
          st_sig & (ex_sig | gene %in% excit_g)  ~ "Shared",
          TRUE ~ "NS"))
      ggplot(sc, aes(x=ex_lfc, y=st_lfc, color=category))+
        geom_hline(yintercept=0,linetype="dashed",color="grey50",linewidth=0.4)+
        geom_vline(xintercept=0,linetype="dashed",color="grey50",linewidth=0.4)+
        geom_point(data=sc%>%dplyr::filter(category=="NS"),size=0.5,alpha=0.25)+
        geom_point(data=sc%>%dplyr::filter(category!="NS"),size=1.2,alpha=0.85)+
        geom_text_repel(data=sc%>%dplyr::filter(category=="Specific")%>%
                          dplyr::slice_max(abs(st_lfc),n=12),
                        aes(label=gene),size=2.3,max.overlaps=20,
                        box.padding=0.3,seed=42,show.legend=FALSE)+
        scale_color_manual(values=c(Specific=SUBTYPE_COLS[st],
                                     Shared="#984EA3",NS="grey82"),name=NULL)+
        theme_classic(base_size=9)+
        theme(plot.title=element_text(face="bold",size=10,hjust=0.5,
                                       color=SUBTYPE_COLS[st]),
              legend.position="top",legend.text=element_text(size=7))+
        labs(title=st, x="Mean log2FC in other IN + ExNs",
             y=paste0("log2FC in ",st))
    })
    p <- wrap_plots(p_list, ncol=ceiling(sqrt(length(subtypes_present))))+
      plot_annotation(
        title=paste0("Subtype-specific ARGs — Scatter — ",comp$label," — ",sx_label),
        subtitle="Specific = sig in subtype, not sig in other INs or ExNs",
        theme=theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
                     plot.subtitle=element_text(size=8,color="grey40",hjust=0.5)))
    save_both(p, stem, w=14, h=12)
  }

  # Up / down split for the heatmaps: a gene's direction is the sign of its largest |log2FC| among the rows given;
  # returns genes ordered up (strongest first) then down (strongest first), and each gene's direction.
  updown_order <- function(d) {
    o <- d %>% dplyr::filter(!is.na(log2FoldChange)) %>% dplyr::group_by(gene) %>%
      dplyr::slice_max(abs(log2FoldChange), n=1, with_ties=FALSE) %>% dplyr::ungroup() %>%
      dplyr::mutate(direction=factor(ifelse(log2FoldChange > 0, "Up", "Down"), levels=c("Up","Down"))) %>%
      dplyr::arrange(direction, dplyr::desc(abs(log2FoldChange)))
    list(genes=o$gene, direction=setNames(o$direction, o$gene))
  }
  updown_theme <- theme(strip.text.y=element_text(face="bold", size=8, angle=0),
                        strip.background.y=element_rect(fill="grey92", colour=NA))

  # ── Plot 2a: LFC heatmap — per age ────────────────────────────────────────
  cat("  [Plot 2a: LFC heatmap — per age]\n")
  for (comp in COMPS) {
    cands <- top_cands_age[[comp$name]]
    if (length(cands)==0) next
    stem <- file.path(out_dir,
                       paste0(fig_date,"_Sec3_LFCHeatmap_",comp$name,"_",sx_label))
    if (skip_fig(stem)) { cat("  Skip:", basename(stem), "\n"); next }
    heat_df <- full_df %>%
      dplyr::filter(sex==(if(sx_key=="combined") "combined" else sx_key),
                     comparison==comp$name, gene %in% cands)
    if (nrow(heat_df)==0) next
    ud <- updown_order(heat_df)                                   # up genes on top, then down
    heat_df <- heat_df %>%
      dplyr::mutate(lfc_cap=pmax(pmin(log2FoldChange,3),-3), sig_dot=sig,
                    celltype=factor(celltype, levels=subtypes_present),
                    direction=ud$direction[gene],
                    gene=factor(gene, levels=rev(ud$genes)))
    p <- ggplot(heat_df, aes(x=celltype, y=gene, fill=lfc_cap))+
      geom_tile(color="white",linewidth=0.5)+
      geom_text(data=heat_df%>%dplyr::filter(sig_dot),
                aes(label="*"),size=3.5,color="black",vjust=0.75,hjust=0.5)+
      scale_fill_gradient2(low="#1E88E5",mid="grey97",high="#E53935",
                            midpoint=0,limits=c(-3,3),oob=scales::squish,
                            name="log2FC\n(EE vs HC)")+
      facet_grid(direction~., scales="free_y", space="free_y")+
      theme_classic(base_size=9)+
      theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
            axis.text.y=element_text(size=7.5,face="italic"),
            axis.title=element_blank(),
            plot.title=element_text(face="bold",size=11,hjust=0.5),
            legend.position="right")+ updown_theme+
      labs(title=paste0("Subtype-specific ARGs — ",comp$label," — ",sx_label),
           subtitle=paste0("* = padj<",PADJ_THR," | Top ",TOP_N," by |LFC|"))
    save_both(p, stem,
              w=max(6,length(subtypes_present)*0.8+2),
              h=max(5,length(cands)*0.22+2))
  }

  # ── Plot 2b: LFC heatmap — union across ages ──────────────────────────────
  cat("  [Plot 2b: LFC heatmap — union]\n")
  stem <- file.path(out_dir,
                     paste0(fig_date,"_Sec3_LFCHeatmap_Union_",sx_label))
  if (!skip_fig(stem)) {
    heat_u <- full_df %>%
      dplyr::filter(sex==(if(sx_key=="combined") "combined" else sx_key),
                     gene %in% top_cands_union)
    ud_u <- updown_order(heat_u)                                  # up genes on top, then down
    heat_u <- heat_u %>%
      dplyr::mutate(lfc_cap=pmax(pmin(log2FoldChange,3),-3), sig_dot=sig,
                    direction=ud_u$direction[gene],
                    age_lbl=factor(AGE_LABELS[comparison],
                                    levels=c("Young EE","Middle EE","Old EE")),
                    celltype=factor(celltype, levels=subtypes_present),
                    gene=factor(gene, levels=rev(ud_u$genes)))
    if (nrow(heat_u)>0) {
      p <- ggplot(heat_u, aes(x=age_lbl, y=gene, fill=lfc_cap))+
        geom_tile(color="white",linewidth=0.4)+
        geom_text(data=heat_u%>%dplyr::filter(sig_dot),
                  aes(label="*"),size=3,color="black",vjust=0.75,hjust=0.5)+
        scale_fill_gradient2(low="#1E88E5",mid="grey97",high="#E53935",
                              midpoint=0,limits=c(-3,3),oob=scales::squish,
                              name="log2FC\n(EE vs HC)")+
        facet_grid(direction~celltype, scales="free_y", space="free_y")+
        theme_classic(base_size=9)+
        theme(axis.text.x=element_text(angle=45,hjust=1,size=7),
              axis.text.y=element_text(size=7,face="italic"),
              axis.title=element_blank(),
              strip.text=element_text(face="bold",size=9),
              panel.spacing=unit(0.3,"lines"),
              plot.title=element_text(face="bold",size=11,hjust=0.5),
              legend.position="right")+ updown_theme+
        labs(title=paste0("Subtype-specific ARGs — Union — ",sx_label),
             subtitle=paste0("Union of top ",TOP_N,
                              " per age | * = padj<",PADJ_THR))
      save_both(p, stem,
                w=max(10,length(subtypes_present)*2.5+2),
                h=max(5,length(top_cands_union)*0.2+2))
    }
  } else cat("  Skip:", basename(stem), "\n")

  # ── Plot 3: Dotplot — per age + union ─────────────────────────────────────
  cat("  [Plot 3: Dotplot — per age + union]\n")
  cand_counts_u <- as(
    GetAssayData(in_obj, assay="RNA", layer="counts")[
      intersect(top_cands_union, rownames(in_obj)),, drop=FALSE],
    "dgCMatrix")
  cands_avail <- rownames(cand_counts_u)

  make_sec3_dotplot <- function(dot_df, cands, ttl, facet_age=FALSE) {
    if (is.null(dot_df)||nrow(dot_df)==0) return(invisible())
    dot_df <- dot_df %>%
      dplyr::mutate(gene=factor(gene, levels=rev(cands)),
                    celltype=factor(celltype, levels=subtypes_present))
    p <- ggplot(dot_df, aes(x=celltype, y=gene))+
      geom_point(data=dot_df%>%dplyr::filter(!sig_dot),
                 aes(size=pct_expr),color="grey88",alpha=0.7)+
      geom_point(data=dot_df%>%dplyr::filter(sig_dot),
                 aes(color=lfc_cap,size=pct_expr),alpha=0.92)+
      scale_color_gradient2(low="#1E88E5",mid="grey97",high="#E53935",
                             midpoint=0,limits=c(-3,3),oob=scales::squish,
                             name="log2FC\n(EE vs HC)")+
      scale_size_continuous(name="% EE\nexpressing",range=c(0.3,10),
                             limits=c(0,100),breaks=c(5,10,25,50,75))+
      theme_classic(base_size=9)+
      theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
            axis.text.y=element_text(size=7.5,face="italic"),
            axis.title=element_blank(),strip.text=element_text(face="bold",size=9),
            legend.position="right",
            panel.grid.major=element_line(color="grey92",linewidth=0.3))+
      guides(color=guide_colorbar(barwidth=0.8,barheight=5,title.position="top"),
             size=guide_legend(title.position="top",
                                override.aes=list(color="grey50")))+
      labs(title=ttl,
           subtitle="Color=DESeq2 log2FC | Size=% EE expressing | Grey=NS")
    if (facet_age) p <- p + facet_wrap(~age_lbl, nrow=1)
    p
  }

  for (comp in COMPS) {
    cands <- intersect(top_cands_age[[comp$name]], cands_avail)
    if (length(cands)==0) next
    stem <- file.path(out_dir,
                       paste0(fig_date,"_Sec3_Dotplot_",comp$name,"_",sx_label))
    if (skip_fig(stem)) { cat("  Skip:", basename(stem), "\n"); next }
    pct_rows <- lapply(subtypes_present, function(st) {
      m  <- meta_full %>%
        dplyr::filter(in_subtype==st, age_lower==comp$age, cond=="EE",
                       if(sx_key=="combined") TRUE else sex_norm==sx_key)
      bc <- intersect(m$barcode, colnames(cand_counts_u))
      if (length(bc)<5) return(NULL)
      pct <- rowMeans(cand_counts_u[cands, bc, drop=FALSE]>0)*100
      data.frame(celltype=st, gene=cands, pct_expr=pct, stringsAsFactors=FALSE)
    })
    pct_df <- dplyr::bind_rows(pct_rows)
    lfc_df <- full_df %>%
      dplyr::filter(sex==(if(sx_key=="combined") "combined" else sx_key),
                     comparison==comp$name, gene %in% cands) %>%
      dplyr::mutate(sig_dot=sig, lfc_cap=pmax(pmin(log2FoldChange,3),-3))
    dot_df <- dplyr::inner_join(lfc_df, pct_df, by=c("celltype","gene"))
    p <- make_sec3_dotplot(dot_df, cands,
                            paste0("Subtype-specific ARGs — ",comp$label," — ",sx_label))
    if (!is.null(p))
      save_both(p, stem,
                w=max(6,length(subtypes_present)*0.9+2),
                h=max(5,length(cands)*0.38+2))
  }

  stem <- file.path(out_dir, paste0(fig_date,"_Sec3_Dotplot_Union_",sx_label))
  if (!skip_fig(stem)) {
    cands_u <- intersect(top_cands_union, cands_avail)
    pct_rows_u <- list()
    for (comp in COMPS) {
      for (st in subtypes_present) {
        m  <- meta_full %>%
          dplyr::filter(in_subtype==st, age_lower==comp$age, cond=="EE",
                         if(sx_key=="combined") TRUE else sex_norm==sx_key)
        bc <- intersect(m$barcode, colnames(cand_counts_u))
        if (length(bc)<5) next
        pct <- rowMeans(cand_counts_u[cands_u, bc, drop=FALSE]>0)*100
        pct_rows_u[[paste(comp$name,st)]] <- data.frame(
          celltype=st, comparison=comp$name,
          gene=cands_u, pct_expr=pct, stringsAsFactors=FALSE)
      }
    }
    pct_df_u <- dplyr::bind_rows(pct_rows_u)
    lfc_df_u <- full_df %>%
      dplyr::filter(sex==(if(sx_key=="combined") "combined" else sx_key),
                     gene %in% cands_u) %>%
      dplyr::mutate(sig_dot=sig, lfc_cap=pmax(pmin(log2FoldChange,3),-3),
                    age_lbl=factor(AGE_LABELS[comparison],
                                    levels=c("Young EE","Middle EE","Old EE")))
    dot_df_u <- dplyr::inner_join(lfc_df_u, pct_df_u,
                                   by=c("celltype","comparison","gene"))
    p <- make_sec3_dotplot(dot_df_u, cands_u,
                            paste0("Subtype-specific ARGs — Union — ",sx_label),
                            facet_age=TRUE)
    if (!is.null(p))
      save_both(p, stem,
                w=max(8,length(subtypes_present)*0.9+3),
                h=max(5,length(cands_u)*0.38+3))
  } else cat("  Skip:", basename(stem), "\n")

  # ── Plot 4a: Ranked bars — per age ────────────────────────────────────────
  cat("  [Plot 4a: Ranked bars — per age]\n")
  for (comp in COMPS) {
    stem <- file.path(out_dir,
                       paste0(fig_date,"_Sec3_RankedBars_",comp$name,"_",sx_label))
    if (skip_fig(stem)) { cat("  Skip:", basename(stem), "\n"); next }
    age_keys <- names(all_specific)[
      startsWith(names(all_specific), paste0(sx_key,"__",comp$name,"__"))]
    age_keys <- age_keys[order(match(sapply(strsplit(age_keys, "__"), `[`, 3), subtypes_present))]   # PV, SST, VIP, LAMP5
    bar_list <- Filter(Negate(is.null), lapply(age_keys, function(k) {
      st  <- strsplit(k,"__")[[1]][3]
      dat <- all_specific[[k]] %>%
        dplyr::slice_max(abs(log2FoldChange), n=TOP_PER_CT) %>%
        dplyr::mutate(direction=ifelse(log2FoldChange>0,"Up","Down"),
                      gene=reorder(gene,log2FoldChange))
      if (nrow(dat)==0) return(NULL)
      ggplot(dat,aes(x=gene,y=log2FoldChange,fill=direction))+
        geom_col(width=0.7)+
        geom_hline(yintercept=0,linewidth=0.4,color="grey20")+
        coord_flip()+
        scale_fill_manual(values=c(Up="#E53935",Down="#1E88E5"),guide="none")+
        theme_classic(base_size=8)+
        theme(axis.text.y=element_text(size=7,face="italic"),
              axis.text.x=element_text(size=7),axis.title.y=element_blank(),
              plot.title=element_text(face="bold",size=9,hjust=0.5,
                                       color=SUBTYPE_COLS[st]))+
        labs(title=st, y="log2FC (EE vs HC)")
    }))
    if (length(bar_list)==0) next
    p <- wrap_plots(bar_list, nrow=1)+
      plot_annotation(
        title=paste0("Subtype-specific ARGs — ",comp$label," — ",
                      sx_label," | Top ",TOP_PER_CT),
        subtitle="Not sig in other IN subtypes or ExNs",
        theme=theme(plot.title=element_text(face="bold",size=11,hjust=0.5),
                     plot.subtitle=element_text(size=8,color="grey40",hjust=0.5)))
    save_both(p, stem, w=max(10,length(bar_list)*3), h=max(5,TOP_PER_CT*0.3+3))
  }

  # ── Plot 4b: Age overlap tile — top genes × ages per subtype ──────────────
  cat("  [Plot 4b: Age overlap tile]\n")
  stem <- file.path(out_dir, paste0(fig_date,"_Sec3_AgeOverlap_",sx_label))
  if (!skip_fig(stem)) {
    tile_list <- Filter(Negate(is.null), lapply(subtypes_present, function(st) {
      rows <- dplyr::bind_rows(lapply(COMPS, function(comp) {
        k <- paste(sx_key, comp$name, st, sep="__")
        if (!k %in% names(all_specific)) return(NULL)
        all_specific[[k]] %>%
          dplyr::select(gene, log2FoldChange) %>%
          dplyr::mutate(age_lbl=AGE_LABELS[comp$name])
      }))
      if (nrow(rows)==0) return(NULL)
      # top genes split evenly into up and down (direction = sign of the gene's largest |log2FC| across ages),
      # each ranked by age breadth then |log2FC|; up genes on top
      gene_ord <- rows %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(n_ages=dplyr::n(),
                         max_lfc=max(abs(log2FoldChange)),
                         direction=ifelse(log2FoldChange[which.max(abs(log2FoldChange))] > 0, "Up", "Down"),
                         .groups="drop") %>%
        dplyr::mutate(direction=factor(direction, levels=c("Up","Down"))) %>%
        dplyr::group_by(direction) %>%
        dplyr::arrange(dplyr::desc(n_ages), dplyr::desc(max_lfc), .by_group=TRUE) %>%
        dplyr::slice_head(n=ceiling(TOP_N/2)) %>% dplyr::ungroup()
      rows %>%
        dplyr::filter(gene %in% gene_ord$gene) %>%
        dplyr::mutate(direction=gene_ord$direction[match(gene, gene_ord$gene)],
                      gene=factor(gene, levels=rev(gene_ord$gene)),
                      age_lbl=factor(age_lbl,
                                      levels=c("Young EE","Middle EE","Old EE")),
                      lfc_cap=pmax(pmin(log2FoldChange,3),-3)) %>%
        ggplot(aes(x=age_lbl, y=gene, fill=lfc_cap))+
          geom_tile(color="white",linewidth=0.5)+
          scale_fill_gradient2(low="#1E88E5",mid="grey97",high="#E53935",
                                midpoint=0,limits=c(-3,3),oob=scales::squish,
                                name="log2FC\n(EE vs HC)")+
          facet_grid(direction~., scales="free_y", space="free_y")+
          theme_classic(base_size=8)+
          theme(axis.text.x=element_text(angle=45,hjust=1,size=7,face="bold"),
                axis.text.y=element_text(size=6.5,face="italic"),
                axis.title=element_blank(),
                plot.title=element_text(face="bold",size=9,hjust=0.5,
                                         color=SUBTYPE_COLS[st]),
                legend.position="right")+ updown_theme+
          labs(title=st)
    }))
    if (length(tile_list)>0) {
      p <- wrap_plots(tile_list, nrow=1, guides="collect")+   # one shared log2FC colour legend
        plot_annotation(
          title=paste0("Subtype-specific ARGs across ages — ",sx_label),
          subtitle=paste0("Top ",ceiling(TOP_N/2)," up and ",ceiling(TOP_N/2),
                           " down per subtype, ranked by age breadth then |log2FC|"),
          theme=theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
                       plot.subtitle=element_text(size=8,color="grey40",hjust=0.5)))
      save_both(p, stem, w=max(10,length(tile_list)*3), h=max(5,TOP_N*0.2+2))
    }
  } else cat("  Skip:", basename(stem), "\n")

  # ── Plot 5: Stacked violin — per age + union ──────────────────────────────
  cat("  [Plot 5: Stacked violin]\n")
  vln_base <- FetchData(in_obj,
                         vars=c("in_subtype","sex_norm","age",
                                intersect(top_cands_union, rownames(in_obj))),
                         layer="data") %>%
    dplyr::filter(in_subtype %in% subtypes_present) %>%
    dplyr::mutate(in_subtype=factor(in_subtype, levels=subtypes_present),
                  age_lower2=tolower(as.character(age)))
  if (sx_key != "combined")
    vln_base <- vln_base %>% dplyr::filter(sex_norm==sx_key)

  for (comp in COMPS) {
    cands_v <- intersect(top_cands_age[[comp$name]], rownames(in_obj))
    if (length(cands_v)==0) next
    stem <- file.path(out_dir,
                       paste0(fig_date,"_Sec3_StackedViolin_",
                               comp$name,"_",sx_label))
    if (skip_fig(stem)) { cat("  Skip:", basename(stem), "\n"); next }
    sub <- vln_base %>%
      dplyr::filter(age_lower2==comp$age) %>%
      tidyr::pivot_longer(cols=all_of(cands_v),
                           names_to="gene", values_to="expression") %>%
      dplyr::mutate(gene=factor(gene, levels=cands_v))
    if (nrow(sub)<200) next
    p <- ggplot(sub, aes(x=in_subtype, y=expression, fill=in_subtype))+
      geom_violin(scale="width",trim=TRUE,linewidth=0.12,color="grey30")+
      facet_grid(gene~., scales="free_y", switch="y")+
      scale_fill_manual(values=SUBTYPE_COLS)+
      scale_x_discrete(expand=expansion(add=0.4))+
      scale_y_continuous(expand=expansion(mult=c(0.05,0.15)))+
      vln_theme()+
      labs(title=paste0("Subtype-specific ARGs — ",comp$label," — ",sx_label),
           subtitle="All conditions | LogNorm RNA")
    save_both(p, stem,
              w=max(5,length(subtypes_present)*0.6+3),
              h=max(8,length(cands_v)*0.45+3))
  }

  cands_vu <- intersect(top_cands_union, rownames(in_obj))
  stem <- file.path(out_dir,
                     paste0(fig_date,"_Sec3_StackedViolin_Union_",sx_label))
  if (length(cands_vu)>0 && !skip_fig(stem)) {
    sub_u <- vln_base %>%
      tidyr::pivot_longer(cols=all_of(cands_vu),
                           names_to="gene", values_to="expression") %>%
      dplyr::mutate(gene=factor(gene, levels=cands_vu))
    p <- ggplot(sub_u, aes(x=in_subtype, y=expression, fill=in_subtype))+
      geom_violin(scale="width",trim=TRUE,linewidth=0.12,color="grey30")+
      facet_grid(gene~., scales="free_y", switch="y")+
      scale_fill_manual(values=SUBTYPE_COLS)+
      scale_x_discrete(expand=expansion(add=0.4))+
      scale_y_continuous(expand=expansion(mult=c(0.05,0.15)))+
      vln_theme()+
      labs(title=paste0("Subtype-specific ARGs — Union — ",sx_label),
           subtitle="All ages and conditions | LogNorm RNA")
    save_both(p, stem,
              w=max(5,length(subtypes_present)*0.6+3),
              h=max(8,length(cands_vu)*0.45+3))
  } else if (skip_fig(stem)) cat("  Skip:", basename(stem), "\n")
}

cat("\nSec 3 complete.\n\n")

### SEC 4: STRIP VOLCANO + ARG TRANSCRIPTOME SIMILARITY
################################################################################

cat("╔══════════════════════════════════════════════════╗\n")
cat("  Sec 4: Strip Volcano + ARG Similarity\n")
cat("╚══════════════════════════════════════════════════╝\n\n")

volcano_comps <- COMPS[c(1,3)]   # Young and Old

violin_offsets <- function(df, y_col, group_col, max_width=MAX_W, n_bins=80) {
  df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::mutate(.bin=as.integer(cut(.data[[y_col]],breaks=n_bins,
                                       include.lowest=TRUE)),
                  .bin=ifelse(is.na(.bin),1L,.bin)) %>%
    dplyr::group_by(.data[[group_col]],.bin) %>%
    dplyr::mutate(.bsize=dplyr::n(),.idx=seq_along(.bin)) %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::mutate(.max_bsize=max(.bsize),
                  .xoff=((.idx-1)-(.bsize-1)/2)/max(1,.bsize-1)*
                    max_width*(.bsize/.max_bsize)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(.xoff=pmax(pmin(.xoff,max_width),-max_width))
}

for (sx_label in c("Male","Female","Combined")) {
  sx_key <- tolower(sx_label)
  cat(sprintf("── %s ──\n", sx_label))

  vol_df <- full_df %>%
    dplyr::filter(sex==(if(sx_key=="combined") "combined" else sx_key),
                   !is.na(log2FoldChange)) %>%
    dplyr::mutate(
      color_cat = dplyr::case_when(
        sig & log2FoldChange>0 ~ "Up",
        sig & log2FoldChange<0 ~ "Down",
        TRUE ~ "NS"),
      lfc_plot  = pmax(pmin(log2FoldChange,LFC_CAP),-LFC_CAP),
      celltype  = factor(celltype, levels=subtypes_present)
    )

  # ── Strip volcano (Young & Old) ────────────────────────────────────────────
  for (comp in volcano_comps) {
    stem <- file.path(out_dir,
                       paste0(fig_date,"_Sec4_StripVolcano_",
                               comp$name,"_",sx_label))
    if (skip_fig(stem)) { cat("  Skip:", basename(stem), "\n"); next }

    sub <- vol_df %>% dplyr::filter(comparison==comp$name)
    if (nrow(sub)==0) next

    if (USE_BEESWARM) {
      p <- ggplot(mapping=aes(x=celltype, y=lfc_plot)) +
        geom_quasirandom(data=sub%>%dplyr::filter(color_cat=="NS"),
                          color=VOL_COLS["NS"],size=0.25,alpha=0.25,width=MAX_W)+
        geom_quasirandom(data=sub%>%dplyr::filter(color_cat!="NS"),
                          aes(color=color_cat),size=0.65,alpha=0.85,width=MAX_W)
    } else {
      sub2 <- violin_offsets(sub,"lfc_plot","celltype") %>%
        dplyr::mutate(.xnum=as.integer(celltype)+.xoff)
      p <- ggplot() +
        geom_point(data=sub2%>%dplyr::filter(color_cat=="NS"),
                   aes(x=.xnum,y=lfc_plot),
                   color=VOL_COLS["NS"],size=0.25,alpha=0.25)+
        geom_point(data=sub2%>%dplyr::filter(color_cat!="NS"),
                   aes(x=.xnum,y=lfc_plot,color=color_cat),
                   size=0.65,alpha=0.85)+
        scale_x_continuous(breaks=seq_along(subtypes_present),
                            labels=subtypes_present,
                            expand=expansion(add=0.6))
    }

    p <- p +
      geom_hline(yintercept=0,linewidth=0.5,color="grey10")+
      scale_color_manual(values=VOL_COLS,breaks=c("Up","Down"),
                          labels=c("up","down"),name=NULL)+
      scale_y_continuous(limits=c(-LFC_CAP-0.3,LFC_CAP+0.3),
                          breaks=c(-4,-2,0,2,4),
                          expand=expansion(mult=c(0,0)))+
      guides(color=guide_legend(override.aes=list(size=2.5,alpha=1),
                                 keywidth=0.4,keyheight=0.4,nrow=1))+
      theme_classic(base_size=10)+
      theme(axis.text.x=element_text(angle=45,hjust=1,size=9,color="grey10"),
            axis.text.y=element_text(size=8.5),
            axis.title.y=element_text(size=10),
            axis.title.x=element_blank(),
            panel.grid=element_blank(),
            plot.title=element_text(size=11,hjust=0.5,face="bold"),
            plot.subtitle=element_text(size=8,hjust=0.5,color="grey40"),
            legend.position="top",legend.justification="right",
            legend.text=element_text(size=8))+
      labs(title   =paste0("IN subtypes — ",comp$label," — ",sx_label),
           subtitle=paste0("DESeq2 EE vs HC | padj<",PADJ_THR,
                            " & |LFC|\u2265",LFC_THR),
           y=bquote(log[2](FC)))

    save_both(p, stem,
              w=max(5, length(subtypes_present)*0.8+2), h=5)
  }

  # ── ARG transcriptome similarity heatmap ─────────────────────────────────
  for (comp in volcano_comps) {
    stem <- file.path(out_dir,
                       paste0(fig_date,"_Sec4_ARGCorrelation_",
                               comp$name,"_",sx_label))
    if (skip_fig(stem)) { cat("  Skip:", basename(stem), "\n"); next }

    lfc_wide <- vol_df %>%
      dplyr::filter(comparison==comp$name, !is.na(log2FoldChange)) %>%
      dplyr::select(gene, celltype, log2FoldChange) %>%
      tidyr::pivot_wider(names_from=celltype, values_from=log2FoldChange,
                          values_fill=0) %>%
      tibble::column_to_rownames("gene")

    if (ncol(lfc_wide)<2) next
    cor_mat <- cor(lfc_wide, method="pearson", use="pairwise.complete.obs")
    cl_ord  <- intersect(subtypes_present, colnames(cor_mat))
    cor_mat <- cor_mat[cl_ord, cl_ord]

    cor_df <- as.data.frame(cor_mat) %>%
      tibble::rownames_to_column("x") %>%
      tidyr::pivot_longer(-x, names_to="y", values_to="correlation") %>%
      dplyr::mutate(x=factor(x,levels=cl_ord),
                    y=factor(y,levels=rev(cl_ord)))

    p <- ggplot(cor_df,aes(x=x,y=y,fill=correlation))+
      geom_tile(color="white",linewidth=0.5)+
      geom_text(aes(label=sprintf("%.2f",correlation)),
                size=3.5,color="grey10")+
      scale_fill_gradient2(low="#1E88E5",mid="grey97",high="#E53935",
                            midpoint=0.5,limits=c(-1,1),
                            name="Pearson\ncorr.")+
      theme_classic(base_size=10)+
      theme(axis.text.x=element_text(angle=45,hjust=1,size=9),
            axis.text.y=element_text(size=9),
            axis.title=element_blank(),
            plot.title=element_text(face="bold",size=11,hjust=0.5),
            plot.subtitle=element_text(size=8,color="grey40",hjust=0.5),
            panel.border=element_rect(color="grey80",fill=NA,linewidth=0.5))+
      labs(title   =paste0("IN subtype ARG similarity — ",
                            comp$label," — ",sx_label),
           subtitle="Pearson r of log2FC profiles | Low off-diagonal = unique ARG program")

    n_cl <- length(cl_ord)
    save_both(p, stem, w=max(4,n_cl*0.8+2), h=max(4,n_cl*0.8+2))
  }
}


################################################################################
### SEC 5: EE AGE TRAJECTORY — OEEYEE / MEEYEE
###   5a. DESeq2 pseudobulk: Old EE vs Young EE, Middle EE vs Young EE
###       (EE cells only, age as contrast)
###   5b. Strip volcano + ARG similarity (as Sec 4)
###   5c. Mirrored stacked barchart: Lost / Maintained ARGs per IN subtype
################################################################################

cat("╔══════════════════════════════════════════════════╗\n")
cat("  Sec 5: EE Age Trajectory (OEEYEE / MEEYEE)\n")
cat("╚══════════════════════════════════════════════════╝\n\n")

TRAJ_COMPS <- list(
  list(name="MEEYEE", age_num="middleage", age_ref="young",
       label="Middle EE vs Young EE"),
  list(name="OEEYEE", age_num="old",       age_ref="young",
       label="Old EE vs Young EE")
)

# ── 5a. DESeq2 pseudobulk for EE trajectory comparisons ─────────────────────
cat("── 5a: DESeq2 (EE trajectory) ──\n")

traj_deseq <- list()

for (sx in c("male","female")) {
  cat(sprintf("  %s\n", sx))
  for (comp in TRAJ_COMPS) {
    cat(sprintf("    %s\n", comp$label))
    for (st in subtypes_present) {
      csv_path <- file.path(pb_dir,
                             sprintf("%s_%s_%s.csv", sx, st, comp$name))
      if (skip_csv(csv_path)) {
        cat(sprintf("      %s — load CSV\n", st))
        res <- read.csv(csv_path, stringsAsFactors=FALSE)
      } else {
        # EE cells only from the two age groups
        meta_sub <- meta_full %>%
          dplyr::filter(sex_norm==sx,
                         age_lower %in% c(comp$age_num, comp$age_ref),
                         cond=="EE",
                         in_subtype==st)
        cat(sprintf("      %s: %d EE cells", st, nrow(meta_sub)))

        uniq_oids <- unique(meta_sub$orig.ident)
        if (length(uniq_oids) < (MIN_PB * 2)) {
          cat(" — skip (too few samples)\n"); next
        }
        n_genes <- nrow(counts_mat)
        pb_mat  <- vapply(uniq_oids, function(oid) {
          bc <- intersect(meta_sub$barcode[meta_sub$orig.ident == oid],
                          colnames(counts_mat))
          if (length(bc)==0) return(rep(0.0, n_genes))
          as.numeric(Matrix::rowSums(counts_mat[, bc, drop=FALSE]))
        }, numeric(n_genes))
        rownames(pb_mat) <- rownames(counts_mat)
        colnames(pb_mat) <- uniq_oids

        col_data <- data.frame(
          sample  = uniq_oids,
          age_grp = factor(
            sapply(uniq_oids, function(oid)
              meta_sub$age_lower[meta_sub$orig.ident == oid][1]),
            levels=c(comp$age_ref, comp$age_num)),  # ref = young
          row.names = uniq_oids,
          stringsAsFactors=FALSE
        )
        cond_tab <- table(col_data$age_grp)
        if (any(cond_tab < MIN_PB)) { cat(" — skip (too few per age)\n"); next }

        dds <- tryCatch(
          DESeqDataSetFromMatrix(round(pb_mat), col_data, ~ age_grp),
          error=function(e){ message("  DESeq2 setup: ", e$message); NULL })
        if (is.null(dds)) { cat(" — DESeq2 setup failed\n"); next }
        keep <- rowSums(counts(dds) >= 5) >= MIN_PB
        dds  <- dds[keep,]
        if (nrow(dds) < 100) { cat(" — too few genes\n"); next }
        dds <- tryCatch(DESeq(dds, quiet=TRUE),
                         error=function(e){ message("  DESeq: ",e$message); NULL })
        if (is.null(dds)) { cat(" — DESeq2 failed\n"); next }

        res <- results(dds,
                        contrast=c("age_grp", comp$age_num, comp$age_ref),
                        independentFiltering=TRUE)
        res <- as.data.frame(res) %>%
          tibble::rownames_to_column("gene") %>%
          dplyr::select(gene, log2FoldChange, padj, baseMean, lfcSE) %>%
          dplyr::filter(!is.na(log2FoldChange))

        cat(sprintf(" — %d genes | %d sig\n",
                    nrow(res),
                    sum(!is.na(res$padj) & res$padj < PADJ_THR &
                          abs(res$log2FoldChange) >= LFC_THR)))
        write.csv(res, csv_path, row.names=FALSE)
      }
      res$celltype   <- st
      res$sex        <- sx
      res$comparison <- comp$name
      traj_deseq[[paste(sx, st, comp$name)]] <- res
    }
  }
}

# Combined: sex as covariate
cat("  combined\n")
for (comp in TRAJ_COMPS) {
  cat(sprintf("    %s\n", comp$label))
  for (st in subtypes_present) {
    csv_path <- file.path(pb_dir,
                           sprintf("combined_%s_%s.csv", st, comp$name))
    if (skip_csv(csv_path)) {
      cat(sprintf("      %s — load CSV\n", st))
      res <- read.csv(csv_path, stringsAsFactors=FALSE)
    } else {
      meta_sub <- meta_full %>%
        dplyr::filter(age_lower %in% c(comp$age_num, comp$age_ref),
                       cond=="EE", in_subtype==st)
      cat(sprintf("      %s: %d EE cells", st, nrow(meta_sub)))

      uniq_oids <- unique(meta_sub$orig.ident)
      if (length(uniq_oids) < (MIN_PB * 2)) {
        cat(" — skip\n"); next
      }
      n_genes <- nrow(counts_mat)
      pb_mat  <- vapply(uniq_oids, function(oid) {
        bc <- intersect(meta_sub$barcode[meta_sub$orig.ident == oid],
                        colnames(counts_mat))
        if (length(bc)==0) return(rep(0.0, n_genes))
        as.numeric(Matrix::rowSums(counts_mat[, bc, drop=FALSE]))
      }, numeric(n_genes))
      rownames(pb_mat) <- rownames(counts_mat)
      colnames(pb_mat) <- uniq_oids

      col_data <- data.frame(
        sample  = uniq_oids,
        age_grp = factor(
          sapply(uniq_oids, function(oid)
            meta_sub$age_lower[meta_sub$orig.ident == oid][1]),
          levels=c(comp$age_ref, comp$age_num)),
        sex = factor(
          sapply(uniq_oids, function(oid)
            meta_sub$sex_norm[meta_sub$orig.ident == oid][1])),
        row.names = uniq_oids,
        stringsAsFactors=FALSE
      )
      cond_tab <- table(col_data$age_grp)
      if (any(cond_tab < MIN_PB)) { cat(" — skip\n"); next }

      n_sexes <- length(levels(droplevels(col_data$sex)))
      design_f <- if (n_sexes > 1) ~ sex + age_grp else ~ age_grp

      dds <- tryCatch(
        DESeqDataSetFromMatrix(round(pb_mat), col_data, design_f),
        error=function(e){ NULL })
      if (is.null(dds)) { cat(" — setup failed\n"); next }
      keep <- rowSums(counts(dds) >= 5) >= MIN_PB
      dds  <- dds[keep,]
      if (nrow(dds) < 100) { cat(" — too few genes\n"); next }
      dds <- tryCatch(DESeq(dds, quiet=TRUE), error=function(e){ NULL })
      if (is.null(dds)) { cat(" — failed\n"); next }

      res <- results(dds,
                      contrast=c("age_grp", comp$age_num, comp$age_ref),
                      independentFiltering=TRUE)
      res <- as.data.frame(res) %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::select(gene, log2FoldChange, padj, baseMean, lfcSE) %>%
        dplyr::filter(!is.na(log2FoldChange))
      cat(sprintf(" — %d genes | %d sig\n", nrow(res),
                  sum(!is.na(res$padj) & res$padj < PADJ_THR &
                        abs(res$log2FoldChange) >= LFC_THR)))
      write.csv(res, csv_path, row.names=FALSE)
    }
    res$celltype   <- st
    res$sex        <- "combined"
    res$comparison <- comp$name
    traj_deseq[[paste("combined", st, comp$name)]] <- res
  }
}

traj_df <- dplyr::bind_rows(traj_deseq) %>%
  dplyr::mutate(sig = !is.na(padj) & padj < PADJ_THR &
                        abs(log2FoldChange) >= LFC_THR)
cat(sprintf("\nTraj DESeq2: %d rows | %d sig\n\n",
            nrow(traj_df), sum(traj_df$sig, na.rm=TRUE)))

# ── 5b. Strip volcano + ARG similarity (OEEYEE / MEEYEE) ────────────────────
cat("── 5b: Strip volcano + ARG similarity ──\n")

for (sx_label in c("Male","Female","Combined")) {
  sx_key <- tolower(sx_label)
  cat(sprintf("  %s\n", sx_label))

  vol_df <- traj_df %>%
    dplyr::filter(sex==(if(sx_key=="combined") "combined" else sx_key),
                   !is.na(log2FoldChange)) %>%
    dplyr::mutate(
      color_cat = dplyr::case_when(
        sig & log2FoldChange>0 ~ "Up",
        sig & log2FoldChange<0 ~ "Down",
        TRUE ~ "NS"),
      lfc_plot  = pmax(pmin(log2FoldChange, LFC_CAP), -LFC_CAP),
      celltype  = factor(celltype, levels=subtypes_present))

  for (comp in TRAJ_COMPS) {
    sub <- vol_df %>% dplyr::filter(comparison==comp$name)
    if (nrow(sub)==0) next

    # Strip volcano
    stem <- file.path(out_dir,
                       paste0(fig_date,"_Sec5_StripVolcano_",
                               comp$name,"_",sx_label))
    if (!skip_fig(stem)) {
      if (USE_BEESWARM) {
        p <- ggplot(mapping=aes(x=celltype, y=lfc_plot)) +
          geom_quasirandom(data=sub%>%dplyr::filter(color_cat=="NS"),
                            color=VOL_COLS["NS"],size=0.25,alpha=0.25,width=MAX_W)+
          geom_quasirandom(data=sub%>%dplyr::filter(color_cat!="NS"),
                            aes(color=color_cat),size=0.65,alpha=0.85,width=MAX_W)
      } else {
        sub2 <- violin_offsets(sub,"lfc_plot","celltype") %>%
          dplyr::mutate(.xnum=as.integer(celltype)+.xoff)
        p <- ggplot()+
          geom_point(data=sub2%>%dplyr::filter(color_cat=="NS"),
                     aes(x=.xnum,y=lfc_plot),
                     color=VOL_COLS["NS"],size=0.25,alpha=0.25)+
          geom_point(data=sub2%>%dplyr::filter(color_cat!="NS"),
                     aes(x=.xnum,y=lfc_plot,color=color_cat),
                     size=0.65,alpha=0.85)+
          scale_x_continuous(breaks=seq_along(subtypes_present),
                              labels=subtypes_present,
                              expand=expansion(add=0.6))
      }
      p <- p +
        geom_hline(yintercept=0,linewidth=0.5,color="grey10")+
        scale_color_manual(values=VOL_COLS,breaks=c("Up","Down"),
                            labels=c("up","down"),name=NULL)+
        scale_y_continuous(limits=c(-LFC_CAP-0.3,LFC_CAP+0.3),
                            breaks=c(-4,-2,0,2,4),
                            expand=expansion(mult=c(0,0)))+
        guides(color=guide_legend(override.aes=list(size=2.5,alpha=1),
                                   keywidth=0.4,keyheight=0.4,nrow=1))+
        theme_classic(base_size=10)+
        theme(axis.text.x=element_text(angle=45,hjust=1,size=9,color="grey10"),
              axis.text.y=element_text(size=8.5),
              axis.title.y=element_text(size=10),
              axis.title.x=element_blank(),
              panel.grid=element_blank(),
              plot.title=element_text(size=11,hjust=0.5,face="bold"),
              plot.subtitle=element_text(size=8,hjust=0.5,color="grey40"),
              legend.position="top",legend.justification="right")+
        labs(title=paste0("IN subtypes — ",comp$label," — ",sx_label),
             subtitle=paste0("DESeq2 EE cells only | padj<",PADJ_THR,
                              " & |LFC|\u2265",LFC_THR),
             y=bquote(log[2](FC)))
      save_both(p, stem, w=max(5,length(subtypes_present)*0.8+2), h=5)
    } else cat("  Skip:", basename(stem), "\n")

    # ARG transcriptome similarity
    stem_cor <- file.path(out_dir,
                           paste0(fig_date,"_Sec5_ARGCorrelation_",
                                   comp$name,"_",sx_label))
    if (!skip_fig(stem_cor)) {
      lfc_wide <- sub %>%
        dplyr::filter(!is.na(log2FoldChange)) %>%
        dplyr::select(gene, celltype, log2FoldChange) %>%
        tidyr::pivot_wider(names_from=celltype, values_from=log2FoldChange,
                            values_fill=0) %>%
        tibble::column_to_rownames("gene")
      if (ncol(lfc_wide) >= 2) {
        cor_mat <- cor(lfc_wide, method="pearson", use="pairwise.complete.obs")
        cl_ord  <- intersect(subtypes_present, colnames(cor_mat))
        cor_mat <- cor_mat[cl_ord, cl_ord]
        cor_df  <- as.data.frame(cor_mat) %>%
          tibble::rownames_to_column("x") %>%
          tidyr::pivot_longer(-x, names_to="y", values_to="correlation") %>%
          dplyr::mutate(x=factor(x,levels=cl_ord), y=factor(y,levels=rev(cl_ord)))
        p <- ggplot(cor_df,aes(x=x,y=y,fill=correlation))+
          geom_tile(color="white",linewidth=0.5)+
          geom_text(aes(label=sprintf("%.2f",correlation)),size=3.5,color="grey10")+
          scale_fill_gradient2(low="#1E88E5",mid="grey97",high="#E53935",
                                midpoint=0.5,limits=c(-1,1),name="Pearson\ncorr.")+
          theme_classic(base_size=10)+
          theme(axis.text.x=element_text(angle=45,hjust=1,size=9),
                axis.text.y=element_text(size=9),axis.title=element_blank(),
                plot.title=element_text(face="bold",size=11,hjust=0.5),
                plot.subtitle=element_text(size=8,color="grey40",hjust=0.5),
                panel.border=element_rect(color="grey80",fill=NA,linewidth=0.5))+
          labs(title=paste0("IN subtype ARG similarity — ",comp$label," — ",sx_label),
               subtitle="Pearson r of log2FC profiles | EE age trajectory")
        n_cl <- length(cl_ord)
        save_both(p, stem_cor, w=max(4,n_cl*0.8+2), h=max(4,n_cl*0.8+2))
      }
    } else cat("  Skip:", basename(stem_cor), "\n")
  }
}

# ── 5c. Mirrored stacked barchart ────────────────────────────────────────────
cat("\n── 5c: Mirrored stacked barchart ──\n")
cat("    (YEEYHC anchor → fate at MEEMHC / OEEOHC, per IN subtype)\n\n")

SEG_ORDER_BAR <- c(
  "Shared Up","Shared Down","Flipped",
  "Lost Up","Lost Down"
)
SEG_COLS_BAR <- c(
  "Shared Up"   = "#238B45",
  "Shared Down" = "#74C476",
  "Flipped"     = "#FD8D3C",
  "Lost Up"     = "#CB181D",
  "Lost Down"   = "#FC3D58"
)
FATE_COMPS <- list(
  list(name="MEEMHC", label="Young → Middle"),
  list(name="OEEOHC", label="Young → Old")
)

classify_fate <- function(yee_df, comp_df, st) {
  y    <- yee_df  %>% dplyr::filter(celltype==st, !is.na(padj))
  d    <- comp_df %>% dplyr::filter(celltype==st, !is.na(padj))
  y_up <- y$gene[y$padj<PADJ_THR & y$log2FoldChange>0]
  y_dn <- y$gene[y$padj<PADJ_THR & y$log2FoldChange<0]
  y_sig<- union(y_up,y_dn)
  d_up <- d$gene[d$padj<PADJ_THR & d$log2FoldChange>0]
  d_dn <- d$gene[d$padj<PADJ_THR & d$log2FoldChange<0]
  d_sig<- union(d_up,d_dn)
  data.frame(
    segment=SEG_ORDER_BAR,
    n_genes=c(
      length(intersect(y_up,d_up)),
      length(intersect(y_dn,d_dn)),
      length(intersect(y_up,d_dn))+length(intersect(y_dn,d_up)),
      length(setdiff(y_up,d_sig)),
      length(setdiff(y_dn,d_sig))),
    stringsAsFactors=FALSE)
}

# Build classification table
bar_rows <- list()
for (sx in c("male","female")) {
  sx_full <- if (sx=="male") "male" else "female"
  yee_sub <- full_df %>%
    dplyr::filter(sex==sx_full, comparison=="YEEYHC")

  for (fc in FATE_COMPS) {
    comp_sub <- full_df %>%
      dplyr::filter(sex==sx_full, comparison==fc$name)
    for (st in subtypes_present) {
      cl <- classify_fate(yee_sub, comp_sub, st)
      cl$celltype   <- st
      cl$comparison <- fc$name
      cl$comp_label <- fc$label
      cl$sex        <- sx
      bar_rows[[paste(sx,fc$name,st)]] <- cl
    }
  }
}

bar_df <- dplyr::bind_rows(bar_rows) %>%
  dplyr::mutate(
    segment    = factor(segment, levels=SEG_ORDER_BAR),
    comp_label = factor(comp_label, levels=c("Young → Middle","Young → Old")),
    celltype   = factor(celltype, levels=subtypes_present)
  )

write.csv(bar_df,
          file.path(out_dir,
                     paste0(fig_date,"_Sec5_BarChart_ClassifiedARGs.csv")),
          row.names=FALSE)

# Build rectangles
build_rect_bar <- function(df) {
  df %>%
    dplyr::filter(n_genes>0) %>%
    dplyr::group_by(sex, celltype, comp_label) %>%
    dplyr::arrange(segment, .by_group=TRUE) %>%
    dplyr::mutate(
      abs_end   = cumsum(n_genes),
      abs_start = abs_end - n_genes,
      xmin      = ifelse(sex=="male", -abs_end,  abs_start),
      xmax      = ifelse(sex=="male", -abs_start, abs_end)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(y_ctr = length(subtypes_present) + 1L - as.integer(celltype))   # PV on top, LAMP5 at bottom
}

# Count plot
stem_bar <- file.path(out_dir,
                       paste0(fig_date,"_Sec5_MirroredBar_Counts"))
if (!skip_fig(stem_bar)) {
  rect_df <- build_rect_bar(bar_df)
  totals  <- rect_df %>%
    dplyr::group_by(sex, celltype, comp_label, y_ctr) %>%
    dplyr::summarise(total=sum(n_genes),
                     tip=ifelse(sex[1]=="male",-max(abs_end),max(abs_end)),
                     .groups="drop")
  x_max <- max(totals$total, na.rm=TRUE)
  x_lim <- x_max * 1.25

  p <- ggplot()+
    geom_rect(data=rect_df,
              aes(xmin=xmin,xmax=xmax,
                  ymin=y_ctr-0.38,ymax=y_ctr+0.38,fill=segment),
              color="white",linewidth=0.2)+
    geom_text(data=totals,
              aes(x=tip+ifelse(sex=="male",-x_max*0.02,x_max*0.02),
                  y=y_ctr,label=total,
                  hjust=ifelse(sex=="male",1,0)),
              size=2.4,fontface="bold")+
    geom_vline(xintercept=0,linewidth=0.5,color="grey20")+
    annotate("text",x=-x_lim*0.6,y=length(subtypes_present)+0.7,
             label="Male",fontface="bold",size=3.5,color="grey20")+
    annotate("text",x= x_lim*0.6,y=length(subtypes_present)+0.7,
             label="Female",fontface="bold",size=3.5,color="grey20")+
    scale_fill_manual(values=SEG_COLS_BAR,name=NULL,breaks=SEG_ORDER_BAR,
                       guide=guide_legend(ncol=2,
                         override.aes=list(size=4,shape=22,stroke=0)))+
    scale_x_continuous(limits=c(-x_lim,x_lim),
                        breaks=pretty(c(-x_max,x_max),n=6),
                        labels=abs,
                        expand=expansion(mult=c(0,0)))+
    scale_y_continuous(breaks=seq_along(subtypes_present),
                        labels=rev(subtypes_present),            # PV on top
                        expand=expansion(add=c(0.6,1.0)))+
    facet_wrap(~comp_label,ncol=1)+
    theme_classic(base_size=10)+
    theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
          plot.subtitle=element_text(size=8,color="grey45",hjust=0.5),
          axis.text.y=element_text(size=10,face="bold",
                                    color=SUBTYPE_COLS[rev(subtypes_present)]),
          axis.text.x=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.line.y=element_blank(),axis.ticks.y=element_blank(),
          strip.text=element_text(face="bold",size=10,color="white"),
          strip.background=element_rect(fill="grey35",color=NA),
          legend.position="bottom",legend.text=element_text(size=8),
          panel.grid.major.x=element_line(color="grey92",linewidth=0.3),
          panel.spacing=unit(0.6,"lines"))+
    labs(title="IN subtype ARG fate — Male vs Female",
         subtitle=paste0(
           "\u2190 Male  |  Female \u2192\n",
           "Green=maintained | Orange=flipped | Red=lost"),
         x="Number of DEGs (Young EE vs HC anchor)", y=NULL)

  save_both(p, stem_bar,
            w=10, h=max(6, length(subtypes_present)*0.9+5))
} else cat("  Skip:", basename(stem_bar), "\n")

# Percentage plot
stem_pct <- file.path(out_dir, paste0(fig_date,"_Sec5_MirroredBar_Pct"))
if (!skip_fig(stem_pct)) {
  pct_df <- bar_df %>%
    dplyr::group_by(sex, celltype, comp_label) %>%
    dplyr::mutate(total_g=sum(n_genes),
                  n_genes=ifelse(total_g>0, n_genes/total_g*100, 0)) %>%
    dplyr::ungroup()
  rect_pct <- build_rect_bar(pct_df)
  totals_pct <- pct_df %>%
    dplyr::group_by(sex, celltype, comp_label) %>%
    dplyr::summarise(total_g=sum(n_genes), .groups="drop") %>%
    dplyr::mutate(total  = sprintf("%.0f%%", total_g),
                  y_ctr  = length(subtypes_present) + 1L - as.integer(factor(celltype,levels=subtypes_present)),
                  tip    = ifelse(sex=="male",-100,100))
  p_pct <- ggplot()+
    geom_rect(data=rect_pct,
              aes(xmin=xmin,xmax=xmax,
                  ymin=y_ctr-0.38,ymax=y_ctr+0.38,fill=segment),
              color="white",linewidth=0.2)+
    geom_text(data=totals_pct,
              aes(x=tip+ifelse(sex=="male",-2,2),y=y_ctr,label=total,
                  hjust=ifelse(sex=="male",1,0)),
              size=2.4,fontface="bold")+
    geom_vline(xintercept=0,linewidth=0.5,color="grey20")+
    annotate("text",x=-60,y=length(subtypes_present)+0.7,
             label="Male",fontface="bold",size=3.5,color="grey20")+
    annotate("text",x= 60,y=length(subtypes_present)+0.7,
             label="Female",fontface="bold",size=3.5,color="grey20")+
    scale_fill_manual(values=SEG_COLS_BAR,name=NULL,breaks=SEG_ORDER_BAR,
                       guide=guide_legend(ncol=2,
                         override.aes=list(size=4,shape=22,stroke=0)))+
    scale_x_continuous(limits=c(-115,115),
                        breaks=c(-100,-75,-50,-25,0,25,50,75,100),
                        labels=function(x) paste0(abs(x),"%"),
                        expand=expansion(mult=c(0,0)))+
    scale_y_continuous(breaks=seq_along(subtypes_present),
                        labels=rev(subtypes_present),            # PV on top
                        expand=expansion(add=c(0.6,1.0)))+
    facet_wrap(~comp_label,ncol=1)+
    theme_classic(base_size=10)+
    theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
          plot.subtitle=element_text(size=8,color="grey45",hjust=0.5),
          axis.text.y=element_text(size=10,face="bold",
                                    color=SUBTYPE_COLS[rev(subtypes_present)]),
          axis.text.x=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.line.y=element_blank(),axis.ticks.y=element_blank(),
          strip.text=element_text(face="bold",size=10,color="white"),
          strip.background=element_rect(fill="grey35",color=NA),
          legend.position="bottom",legend.text=element_text(size=8),
          panel.grid.major.x=element_line(color="grey92",linewidth=0.3),
          panel.spacing=unit(0.6,"lines"))+
    labs(title="IN subtype ARG fate (%) — Male vs Female",
         subtitle=paste0(
           "\u2190 Male  |  Female \u2192\n",
           "Each bar = 100% of Young EE vs HC ARGs\n"),
         x="% of Young EE vs HC DEGs", y=NULL)

  save_both(p_pct, stem_pct, w=10, h=max(6,length(subtypes_present)*0.9+5))
} else cat("  Skip:", basename(stem_pct), "\n")

cat("\nSec 5 complete.\n\n")


cat("\n\u2554══════════════════════════════════════════════════\u2557\n")
cat("  ALL DONE\n")
cat("\u255a══════════════════════════════════════════════════\u255d\n")
cat("Outputs:", out_dir, "\n\n")
cat("Sec 1: DESeq2 CSVs (EE vs HC)     →", pb_dir, "\n")
cat("Sec 2: IEG DotViolin              → *_Sec2_*\n")
cat("Sec 3: ARG Discovery              → *_Sec3_*\n")
cat("Sec 4: Strip Volcano + Similarity → *_Sec4_*\n")
cat("Sec 5: EE Trajectory + Barchart   → *_Sec5_*\n")
