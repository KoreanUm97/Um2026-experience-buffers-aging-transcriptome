################################################################################
### ED8f_ARG_Classes.R
###
### Extended Data Fig. 8f: number of DREADD ARGs per cell type in five classes (female only), covariate approach.
###   Young        YC21YPBS padj < 0.05
###   Old          OC21OPBS padj < 0.05
###   Maintained   significant at both ages, same direction
###   Lost         significant in young, not significant in old (either direction)
###   Age-emergent significant in old, not significant in young
### No fold-change cutoff (as for every other ARG set in the paper).
###
### Reads the DESeq2 tables written by DREADD_DESeq2.R:
###   <DESEQ_ROOT>/DESeq2_mScarletPos_Approach<X>_MergedIN/<date>/<date>_{YC21YPBS,OC21OPBS}_ALL.csv  (latest date used)
### APPROACHES = "B" for ED8f; add "D" for the same plot with the trimming + downsampling approach.
### Run:  Rscript ED8f_ARG_Classes.R     (set DESEQ_ROOT below first)
################################################################################

suppressPackageStartupMessages({
  library(conflicted); library(dplyr); library(tidyr); library(ggplot2)
  conflicts_prefer(dplyr::select, .quiet = TRUE); conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicts_prefer(dplyr::mutate, .quiet = TRUE); conflicts_prefer(base::intersect, .quiet = TRUE)
  conflicts_prefer(base::setdiff, .quiet = TRUE); conflicts_prefer(base::union, .quiet = TRUE)
})

################################################################################
### SETTINGS
################################################################################

DESEQ_ROOT <- "PATH_TO_DREADD_OUTPUT"                     # same as OUT_ROOT in DREADD_DESeq2.R
OUT_DIR    <- file.path(DESEQ_ROOT, "Paper_Figures", "Figure4", "ARG_Classes")
APPROACHES <- c("B")                                      # B = covariate (ED8f); c("B", "D") for both
APPROACH_DIRS <- c(B = "DESeq2_mScarletPos_ApproachB_MergedIN", D = "DESeq2_mScarletPos_ApproachD_MergedIN")
APPROACH_LABS <- c(B = "Covariate", D = "Downsampling")
PADJ_THR <- 0.05

CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN","Astrocytes","Oligodendrocytes","OPC",
              "Microglia","Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell","Macrophage","Epithelial","T Cell")
CLASS_COLS <- c(Young = "#A6CEE3", Old = "#33A02C", Maintained = "#2166AC", Lost = "#D6604D", `Age-Emergent` = "#4DAF4A")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

latest_table <- function(dir, comp) {
  f <- list.files(dir, pattern = paste0("_", comp, "_ALL\\.csv$"), recursive = TRUE, full.names = TRUE)
  if (!length(f)) stop("No ", comp, "_ALL.csv under ", dir)
  f <- sort(f, decreasing = TRUE)[1]; cat("  ", comp, ":", f, "\n"); read.csv(f, stringsAsFactors = FALSE)
}
sig_genes <- function(d, ct, dir) {
  d <- d[d$celltype == ct & !is.na(d$padj) & !is.na(d$log2FoldChange) & d$padj < PADJ_THR, ]
  if (dir == "up") d$gene[d$log2FoldChange > 0] else if (dir == "down") d$gene[d$log2FoldChange < 0] else d$gene
}

################################################################################
### COUNTS + DOT PLOT
################################################################################

for (ap in APPROACHES) {
  cat("\n=== Approach", ap, "===\n")
  yc21 <- latest_table(file.path(DESEQ_ROOT, APPROACH_DIRS[[ap]]), "YC21YPBS")
  oc21 <- latest_table(file.path(DESEQ_ROOT, APPROACH_DIRS[[ap]]), "OC21OPBS")
  cts  <- c(intersect(CT_ORDER, unique(yc21$celltype)), sort(setdiff(unique(yc21$celltype), CT_ORDER)))

  class_df <- dplyr::bind_rows(lapply(cts, function(ct) {
    yu <- sig_genes(yc21, ct, "up"); yd <- sig_genes(yc21, ct, "down"); ya <- c(yu, yd)
    ou <- sig_genes(oc21, ct, "up"); od <- sig_genes(oc21, ct, "down"); oa <- c(ou, od)
    data.frame(celltype = ct,
               n_young_all = length(ya), n_young_up = length(yu), n_young_down = length(yd),
               n_old_all   = length(oa), n_old_up   = length(ou), n_old_down   = length(od),
               n_mnt_up  = length(intersect(yu, ou)), n_mnt_down  = length(intersect(yd, od)),
               n_lost_up = length(setdiff(yu, oa)),   n_lost_down = length(setdiff(yd, oa)),
               n_gain_up = length(setdiff(ou, ya)),   n_gain_down = length(setdiff(od, ya)),
               stringsAsFactors = FALSE)
  })) %>% dplyr::mutate(n_mnt_all = n_mnt_up + n_mnt_down, n_lost_all = n_lost_up + n_lost_down,
                        n_gain_all = n_gain_up + n_gain_down)
  utils::write.csv(class_df, file.path(OUT_DIR, paste0("ARG_Class_Counts_Approach", ap, ".csv")), row.names = FALSE)
  print(class_df %>% dplyr::select(celltype, n_young_all, n_old_all, n_mnt_all, n_lost_all, n_gain_all), row.names = FALSE)

  pd <- class_df %>%
    dplyr::select(celltype, Young = n_young_all, Old = n_old_all, Maintained = n_mnt_all, Lost = n_lost_all,
                  `Age-Emergent` = n_gain_all) %>%
    tidyr::pivot_longer(-celltype, names_to = "class", values_to = "n") %>%
    dplyr::mutate(class = factor(class, levels = names(CLASS_COLS)), celltype = factor(celltype, levels = rev(cts)))
  p <- ggplot(pd, aes(x = class, y = celltype, size = n, color = class)) +
    geom_point(alpha = 0.85) +
    geom_text(aes(label = n), size = 2.5, color = "grey20", vjust = -1.2) +
    scale_size_continuous(range = c(1, 12), name = "n genes") +
    scale_color_manual(values = CLASS_COLS, guide = "none") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 8.5, color = "grey40", hjust = 0.5), legend.position = "top",
          axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 9)) +
    labs(title = paste0(APPROACH_LABS[[ap]], "  DREADD: ARG counts per class per cell type"),
         subtitle = paste0("padj < ", PADJ_THR, " | size and label = n genes | Female only"), x = NULL, y = NULL)
  f <- file.path(OUT_DIR, paste0("ARG_Class_Dotplot_Approach", ap))
  ggsave(paste0(f, ".pdf"), p, width = 9, height = max(6, length(cts) * 0.5 + 3), device = "pdf")
  ggsave(paste0(f, ".png"), p, width = 9, height = max(6, length(cts) * 0.5 + 3), dpi = 600, bg = "white")
  cat("  saved:", basename(f), "\n")
}

cat("\n=== DONE ===\nOutputs:", OUT_DIR, "\n")
