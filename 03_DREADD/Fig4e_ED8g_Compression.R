################################################################################
### Fig4e_ED8g_Compression.R
###
### Dynamic range compression of the young DREADD ARG program (female only).
###   Fig. 4e   Approach D (trimming + downsampling)
###   ED8g      Approach B (mScarlet covariate)
### Young ARGs = YC21YPBS padj < 0.05, Up and Down separately, per cell type (a cell type is used if either direction has
### >= MIN_GENES genes). Each gene's log2FC (C21 vs PBS) in young (YC21YPBS) and old (OC21OPBS) mice.
### Plots (FitWindow): box per age (median, IQR, whiskers 1.5 x IQR, no outliers) + mean line; y axis = whisker range
### (and 0); grey lines only for genes whose whole trajectory lies inside that window (boxes and means use all genes).
### Stars (SHOW_STARS) = fdr_loc: * < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001; non-significant = no bracket.
### Stats (Compression_Stats_Approach<X>.csv): two-sided paired Wilcoxon signed-rank test, young vs old log2FC, on all
### genes in the set (p_loc_wilcox); BH-corrected across all cell types x both directions within each approach (fdr_loc).
### change = compressed (|mean old| < |mean young|) or increased; the test is two-sided, so check it with the star.
### Also IQR ratio young/old with 1000-gene bootstrap and an F test (not shown in the figures).
###
### Reads the DESeq2 tables written by DREADD_DESeq2.R:
###   <DESEQ_ROOT>/DESeq2_mScarletPos_Approach{B,D}_MergedIN/<date>/<date>_{YC21YPBS,OC21OPBS}_ALL.csv  (latest date used)
### Run:  Rscript Fig4e_ED8g_Compression.R     (set DESEQ_ROOT and OUT_ROOT below first)
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(conflicted)
  conflicts_prefer(dplyr::select, .quiet = TRUE); conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicts_prefer(base::intersect, .quiet = TRUE); conflicts_prefer(base::setdiff, .quiet = TRUE)
  conflicts_prefer(base::union, .quiet = TRUE)
})

################################################################################
### SETTINGS
################################################################################

DESEQ_ROOT <- "PATH_TO_DREADD_OUTPUT"                       # same as OUT_ROOT in DREADD_DESeq2.R
OUT_ROOT   <- file.path(DESEQ_ROOT, "Paper_Figures", "Figure4", "Compression")
APPROACHES <- c(B = "DESeq2_mScarletPos_ApproachB_MergedIN",   # covariate (ED8g)
                D = "DESeq2_mScarletPos_ApproachD_MergedIN")   # trimming + downsampling (Fig. 4e)

PADJ_THR  <- 0.05
MIN_GENES <- 10
N_BOOT    <- 1000
ZOOM_INCLUDE_ZERO <- TRUE        # keep log2FC = 0 (dashed line) in view
SHOW_STARS <- TRUE               # draw the BH-corrected star (fdr_loc) as a Young-Old bracket above the boxes; FALSE = no bracket
COL_UP <- "#D6604D"; COL_DN <- "#4393C3"
SKIP_EXISTING <- TRUE

CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN","Astrocytes","Oligodendrocytes","OPC",
              "Microglia","Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell","Macrophage","Epithelial","T Cell")

################################################################################
### HELPERS
################################################################################

latest_table <- function(dir, comp) {
  f <- list.files(dir, pattern = paste0("_", comp, "_ALL\\.csv$"), recursive = TRUE, full.names = TRUE)
  if (!length(f)) stop("No ", comp, "_ALL.csv under ", dir)
  f <- sort(f, decreasing = TRUE)[1]; cat("  ", comp, ":", f, "\n"); read.csv(f, stringsAsFactors = FALSE)
}

# Box statistics exactly as geom_boxplot computes them (quantile type 7, whiskers at 1.5 x IQR)
box_stats <- function(v) {
  v <- v[!is.na(v)]; q <- stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE); iqr <- q[3] - q[1]
  c(ymin = min(v[v >= q[1] - 1.5 * iqr]), ymax = max(v[v <= q[3] + 1.5 * iqr]))
}

# df: gene, x (1 = young, 2 = old), y (log2FC)
sig_label <- function(p) if (is.na(p)) "ns" else if (p < 1e-4) "****" else if (p < 1e-3) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"

make_window_plot <- function(df, col, title, subtitle, star = NA_character_) {
  df <- df[!is.na(df$y), ]
  if (!nrow(df)) return(NULL)
  bs <- sapply(split(df$y, df$x), box_stats)
  lo <- min(bs["ymin", ]); hi <- max(bs["ymax", ])
  if (ZOOM_INCLUDE_ZERO) { lo <- min(lo, 0); hi <- max(hi, 0) }
  hair  <- (hi - lo) * 0.002                     # hairline so whisker tips are never clipped
  draw  <- SHOW_STARS && !is.na(star) && star != "ns"
  y_br  <- hi + (hi - lo) * 0.05; top <- if (draw) hi + (hi - lo) * 0.14 else hi + hair
  means <- stats::aggregate(y ~ x, data = df, FUN = mean)
  inside  <- tapply(seq_len(nrow(df)), df$gene, function(i) length(unique(df$x[i])) == 2 && all(df$y[i] >= lo & df$y[i] <= hi))
  line_df <- df[df$gene %in% names(inside)[inside], ]
  p <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.6, color = "grey30", linetype = "dashed") +
    geom_line(data = line_df, aes(x = x, y = y, group = gene), color = "grey70", alpha = 0.35, linewidth = 0.25) +
    geom_boxplot(data = df, aes(x = x, y = y, group = x), fill = adjustcolor(col, 0.30), color = col,
                 linewidth = 0.7, outlier.shape = NA, width = 0.45) +
    geom_line(data = means, aes(x = x, y = y), color = col, linewidth = 1.0) +
    geom_point(data = means, aes(x = x, y = y), color = col, size = 2.4) +
    scale_x_continuous(breaks = 1:2, labels = c("Young C21", "Old C21"), expand = expansion(add = 0.45)) +
    scale_y_continuous(expand = expansion(0)) +
    coord_cartesian(ylim = c(lo - hair, top), clip = "on") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 7.5, color = "grey40", hjust = 0.5),
          axis.text = element_text(color = "black")) +
    labs(title = title, x = NULL, y = expression(log[2]*FC~"(C21 vs PBS)"),
         subtitle = paste0(subtitle, "\nn = ", length(unique(df$gene)), " genes (",
                           length(unique(line_df$gene)), " lines shown)"))
  if (draw) {
    tk <- (hi - lo) * 0.02
    p <- p + annotate("segment", x = 1, xend = 2, y = y_br, yend = y_br, linewidth = 0.5) +
      annotate("segment", x = c(1, 2), xend = c(1, 2), y = y_br, yend = y_br - tk, linewidth = 0.5) +
      annotate("text", x = 1.5, y = y_br + tk * 0.6, label = star, size = 3.4, vjust = 0)
  }
  p
}

compute_stats <- function(y_df, o_df, ct, dir) {
  y_lfc <- y_df$log2FoldChange; o_lfc <- o_df$log2FoldChange
  if (length(y_lfc) < 5 || length(o_lfc) < 5) return(NULL)
  common <- intersect(y_df$gene, o_df$gene)
  p_loc <- if (length(common) >= 5)
    tryCatch(wilcox.test(y_df$log2FoldChange[match(common, y_df$gene)], o_df$log2FoldChange[match(common, o_df$gene)],
                         paired = TRUE)$p.value, error = function(e) NA_real_) else NA_real_
  set.seed(42)
  obs  <- IQR(y_lfc) / max(IQR(o_lfc), 1e-9)
  boot <- replicate(N_BOOT, IQR(sample(y_lfc, replace = TRUE)) / max(IQR(sample(o_lfc, replace = TRUE)), 1e-9))
  data.frame(celltype = ct, direction = dir, n_young = length(y_lfc), n_old = length(o_lfc), n_paired = length(common),
             mean_young = mean(y_lfc), mean_old = mean(o_lfc), p_loc_wilcox = p_loc,
             iqr_young = IQR(y_lfc), iqr_old = IQR(o_lfc), iqr_ratio = obs,
             iqr_ci_lo = unname(quantile(boot, 0.025)), iqr_ci_hi = unname(quantile(boot, 0.975)),
             p_boot_iqr = mean(boot <= 1),
             p_ftest = tryCatch(var.test(y_lfc, o_lfc, alternative = "greater")$p.value, error = function(e) NA_real_),
             stringsAsFactors = FALSE)
}

################################################################################
### MAIN
################################################################################

for (ap in names(APPROACHES)) {
  cat("\n=== Approach", ap, "===\n")
  yc21 <- latest_table(file.path(DESEQ_ROOT, APPROACHES[[ap]]), "YC21YPBS")
  oc21 <- latest_table(file.path(DESEQ_ROOT, APPROACHES[[ap]]), "OC21OPBS")
  OUT  <- file.path(OUT_ROOT, paste0("Approach", ap)); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

  cts <- c(intersect(CT_ORDER, unique(yc21$celltype)), sort(setdiff(unique(yc21$celltype), CT_ORDER)))
  sets <- list(); stats <- list()
  for (ct in cts) {
    y_ct <- yc21[yc21$celltype == ct & !is.na(yc21$log2FoldChange), ]
    o_ct <- oc21[oc21$celltype == ct & !is.na(oc21$log2FoldChange), ]
    up <- y_ct$gene[!is.na(y_ct$padj) & y_ct$padj < PADJ_THR & y_ct$log2FoldChange > 0]
    dn <- y_ct$gene[!is.na(y_ct$padj) & y_ct$padj < PADJ_THR & y_ct$log2FoldChange < 0]
    if (length(up) < MIN_GENES && length(dn) < MIN_GENES) next
    for (dir in c("Up", "Down")) {
      g <- if (dir == "Up") up else dn
      if (!length(g)) next
      yd <- y_ct[y_ct$gene %in% g, ]; od <- o_ct[o_ct$gene %in% g, ]
      sets[[paste(ct, dir)]] <- list(ct = ct, dir = dir, y = yd, o = od)
      stats[[paste(ct, dir)]] <- compute_stats(yd, od, ct, dir)
    }
  }
  st <- dplyr::bind_rows(stats) %>%
    dplyr::mutate(fdr_loc = p.adjust(p_loc_wilcox, "BH"), fdr_boot = p.adjust(p_boot_iqr, "BH"),
                  fdr_ftest = p.adjust(p_ftest, "BH"),
                  # direction of the age change: compressed = old response closer to 0 than young
                  change = ifelse(abs(mean_old) < abs(mean_young), "compressed", "increased"))
  utils::write.csv(st, file.path(OUT, paste0("Compression_Stats_Approach", ap, ".csv")), row.names = FALSE)

  for (s in sets) {
    out_f <- file.path(OUT, paste0(gsub("[^A-Za-z0-9]+", "_", s$ct), "_", s$dir, "_GreyLines_FitWindow.pdf"))
    if (SKIP_EXISTING && file.exists(out_f)) { cat("  exists, skipped:", basename(out_f), "\n"); next }
    if (!nrow(s$y) || !nrow(s$o)) { cat("  no genes at one age, skipped:", s$ct, s$dir, "\n"); next }
    df <- rbind(data.frame(gene = s$y$gene, x = rep(1L, nrow(s$y)), y = s$y$log2FoldChange),
                data.frame(gene = s$o$gene, x = rep(2L, nrow(s$o)), y = s$o$log2FoldChange))
    fl <- st$fdr_loc[st$celltype == s$ct & st$direction == s$dir]
    p <- make_window_plot(df, if (s$dir == "Up") COL_UP else COL_DN, s$ct, paste0(s$dir, " in young C21 (YC21YPBS)"),
                          star = if (length(fl)) sig_label(fl) else NA_character_)
    if (is.null(p)) next
    ggsave(out_f, p, width = 3.2, height = 4, device = "pdf")
    cat("  saved:", basename(out_f), "\n")
  }
  cat("Stats:\n"); print(st %>% dplyr::filter(celltype %in% c("DG ExNs","CA1 ExNs","CA3 ExNs","IN")) %>%
                           dplyr::select(celltype, direction, n_paired, mean_young, mean_old, change, fdr_loc), row.names = FALSE)
}

cat("\n=== DONE ===\nOutputs:", OUT_ROOT, "\n")
