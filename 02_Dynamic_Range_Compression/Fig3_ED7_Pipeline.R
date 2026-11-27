################################################################################
### Fig3_ED7_Pipeline.R
###
### Figure 3a and Extended Data Fig. 7a-c: dynamic-range compression of the young NE program.
### Reads only the per-sex DESeq2 tables written by 01_DESeq2/Fig2_ED4_Pipeline.R (DESeq2_Results/); no Seurat object,
### DESeq2 is not rerun.
### Run:  Rscript Fig3_ED7_Pipeline.R     (set BASE below first)
###
###   Step      Panel        What it does
###   --------  -----------  ------------------------------------------------------------------------------------------
###   3b        Fig. 3a      mean log2FC trajectory of young NE genes (YEEYHC padj < 0.05, Up / Down separately) at
###                          young, middle and old age, one line per cell type
###                          (3b_Trajectory_Boxplots/NoLines/Fig3b_AllCelltypes_MeanTrajectory_<sex>.pdf)
###             ED7a,b       per cell type and sex: log2FC box per age + mean line for the same genes
###                          (NoLines/Fig3b_<celltype>_<sex>_UP_only.pdf / _DOWN_only.pdf)
###                          Paired two-sided Wilcoxon between ages (>= 5 paired genes), BH-corrected across all
###                          cell types x the 3 age comparisons within each sex x direction
###                          (NoLines/Fig3b_Wilcoxon_Stats_BH.csv; same tests and correction as Fig3b_Boxplots_BH.R)
###   compress  ED7c         compression score = mean|log2FC young| / mean|log2FC aged| over the young NE genes, per
###                          Up / Down arm, averaged over arms; Young->Middle and Young->Old; 95% CI from N_BOOT gene
###                          bootstraps (Compression_LFCRatio/Compression_<age>_SexColours.pdf)
###
### Each sex is analysed on its own. Every step writes to its own subfolder of Paper_Figures/Figure3/, skips files that
### already exist (SKIP_EXISTING), and a failing step does not stop the others. Log: Figure3/Logs/.
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(conflicted)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::mutate, .quiet=TRUE); conflicts_prefer(dplyr::arrange, .quiet=TRUE)
  conflicts_prefer(base::intersect, .quiet=TRUE); conflicts_prefer(base::setdiff, .quiet=TRUE)
  conflicts_prefer(base::union, .quiet=TRUE)
})

################################################################################
### SETTINGS
################################################################################

RUN <- c("3b", "compress")
BASE       <- "PATH_TO_OUTPUT"                                  # same BASE as 01_DESeq2/Fig2_ED4_Pipeline.R
DESEQ_DIR  <- file.path(BASE, "DESeq2_Results")                  # <sex>_ALL_<comparison>.csv
FIG_ROOT   <- file.path(BASE, "Paper_Figures", "Figure3")        # one subfolder per step
dir.create(FIG_ROOT, recursive=TRUE, showWarnings=FALSE)

SEXES     <- c("male", "female")
PADJ_THR  <- 0.05         # BH, young NE genes
MIN_GENES <- 10           # a cell type is plotted if either Up or Down arm has >= this many young NE genes

# -- 3b ------------------------------------------------------------------------------------------------------------
# y axis of the arm plots: "whisker" = fitted to the box whiskers (used for ED7a,b);
#                          "original" = 95th percentile x 1.4, lower edge at the most extreme gene
NOLINE_Y_ZOOM <- "whisker"
# Significance brackets on the arm plots (BH-corrected stars). FALSE = no brackets and no space kept for them; the
# tests are still written to Fig3b_Wilcoxon_Stats_BH.csv. (ED7a,b were drawn with FALSE and the stars added from the CSV.)
SHOW_BRACKETS <- FALSE
ZOOM_INCLUDE_ZERO <- TRUE # zoomed axes keep log2FC = 0 (dashed line) in view
BRACKET_FIRST <- 0.03     # gap between the data and the first bracket row (fraction of the plotted range)
BRACKET_ROW   <- 0.07     # height of each bracket row

# -- compress ------------------------------------------------------------------------------------------------------
N_BOOT <- 1000            # bootstrap resamples (genes) for the 95% CI

CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN",
              "Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal",
              "Cajal-Retzius Cell","Macrophage","Epithelial","T Cell")
CT_COLS <- c("DG ExNs"="#E53935","CA1 ExNs"="#CC0000","CA2 ExNs"="#8E24AA","CA3 ExNs"="#5E35B1","Sub ExNs"="#80D8FF",
             "IN"="#00897B","Astrocytes"="#43A047","Oligodendrocytes"="#1E88E5","OPC"="#00ACC1","Microglia"="#7CB342",
             "Endothelial"="#FB8C00","Pericyte"="#6D4C41","APC"="#546E7A","Stromal"="#27408B","Macrophage"="#C62828",
             "Cajal-Retzius Cell"="#FDD835","T Cell"="#757575","Ambiguous"="#BDBDBD")
AGE_COLS <- c(Middle="#00897B", Old="#FF1744")
SEX_COLS <- c(male="#1565C0", female="#FF1744")
NEEDED_COMPS <- c("YEEYHC","MEEMHC","OEEOHC")

SECTION_INFO <- c(
  "3b"       = "Fig. 3a mean trajectories + ED7a,b per-cell-type boxplots (young NE genes at middle / old age)",
  "compress" = "ED7c compression score (mean|LFC| young / aged)")

# -- logging, skip-if-exists ---------------------------------------------------------------------------------------
LOG_DIR <- file.path(FIG_ROOT, "Logs"); dir.create(LOG_DIR, recursive=TRUE, showWarnings=FALSE)
LOG <- file.path(LOG_DIR, paste0("Figure3_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
say <- function(...) { m <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", ...); cat(m, "\n"); cat(m, "\n", file=LOG, append=TRUE) }
banner <- function(id) { line <- strrep("=", 90)
  cat("\n", line, "\n", sep=""); say(sprintf("STEP %s  |  %s", id, SECTION_INFO[[id]])); cat(line, "\n\n", sep="") }
SKIP_EXISTING <- TRUE      # FALSE = redraw and overwrite everything
ggsave <- function(filename, plot=ggplot2::last_plot(), ...) {
  if (SKIP_EXISTING && file.exists(filename)) { cat("  exists, skipped:", basename(filename), "\n"); return(invisible(filename)) }
  ggplot2::ggsave(filename, plot, ...); cat("  saved:", basename(filename), "\n"); invisible(filename)
}
write.csv <- function(x, file="", ...) {
  if (SKIP_EXISTING && is.character(file) && nzchar(file) && file.exists(file)) {
    cat("  exists, skipped:", basename(file), "\n"); return(invisible(NULL)) }
  utils::write.csv(x, file, ...)
}
save_both <- function(p, stem, w, h, device_pdf="pdf") {      # PDF + 600 dpi PNG
  ggsave(paste0(stem, ".pdf"), p, width=w, height=h, device=device_pdf, limitsize=FALSE)
  if (!requireNamespace("ragg", quietly=TRUE) && isTRUE(capabilities("cairo")))   # ggsave uses ragg when installed
    ggsave(paste0(stem, ".png"), p, width=w, height=h, device="png", dpi=600, bg="white", limitsize=FALSE, type="cairo")
  else ggsave(paste0(stem, ".png"), p, width=w, height=h, device="png", dpi=600, bg="white", limitsize=FALSE)
}
ARROW <- if (isTRUE(l10n_info()[["UTF-8"]])) "\u2192" else "->"
step_dir <- function(...) { d <- file.path(FIG_ROOT, ...); dir.create(d, recursive=TRUE, showWarnings=FALSE); d }

# -- shared data helpers -------------------------------------------------------------------------------------------
tables_ok <- function() all(file.exists(file.path(DESEQ_DIR, paste0(rep(SEXES, each=length(NEEDED_COMPS)), "_ALL_", NEEDED_COMPS, ".csv"))))
.deseq_cache <- new.env()
load_deseq <- function(sx, comp) {          # one DESeq2 table (sex x comparison), cached for the whole run
  key <- paste(sx, comp)
  if (!exists(key, envir=.deseq_cache, inherits=FALSE)) {
    f <- file.path(DESEQ_DIR, paste0(sx, "_ALL_", comp, ".csv"))
    assign(key, if (file.exists(f)) read.csv(f, stringsAsFactors=FALSE) else NULL, envir=.deseq_cache)
  }
  get(key, envir=.deseq_cache)
}
ct_levels_of <- function(cts) c(intersect(CT_ORDER, cts), sort(setdiff(cts, CT_ORDER)))
sig_genes <- function(d, dir=0) {            # dir: +1 Up, -1 Down, 0 either
  ok <- !is.na(d$padj) & d$padj < PADJ_THR
  if (dir > 0) ok <- ok & d$log2FoldChange > 0
  if (dir < 0) ok <- ok & d$log2FoldChange < 0
  unique(d$gene[ok])
}
sig_label <- function(p) {
  if (is.na(p)) return("ns"); if (p < 1e-4) return("****"); if (p < 1e-3) return("***")
  if (p < 0.01) return("**"); if (p < 0.05) return("*"); "ns"
}

AGE_COMPS <- list(Middle=c(ee="MEEMHC"), Old=c(ee="OEEOHC"))


################################################################################
### PRE-FLIGHT
################################################################################

say("Figure 3 pipeline started. Steps: ", paste(intersect(names(SECTION_INFO), RUN), collapse=", "), "   Log: ", LOG)
if (!tables_ok()) stop("Missing DESeq2 tables in ", DESEQ_DIR, ": need ", paste(NEEDED_COMPS, collapse=", "),
                       " for both sexes (run the deseq2 step of 01_DESeq2/Fig2_ED4_Pipeline.R).")

SECTIONS <- list()

################################################################################
### STEP 3b: Fig. 3a MEAN TRAJECTORIES + ED7a,b BOXPLOTS
###   Anchor = YEEYHC-significant genes (Up and Down separately) tracked through MEEMHC and OEEOHC.
################################################################################

SECTIONS[["3b"]] <- function() {
  OUT_DIR  <- step_dir("3b_Trajectory_Boxplots", "NoLines")
  fig_date <- "Fig3b"
  age_order <- c("YEEYHC", "MEEMHC", "OEEOHC")
  COL_UP <- "#D6604D"; COL_DN <- "#4393C3"

  # Box statistics exactly as geom_boxplot computes them (quantile type 7, whiskers at 1.5 x IQR)
  box_stats <- function(v) { v <- v[!is.na(v)]; q <- stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE); iqr <- q[3] - q[1]
    c(ymin = min(v[v >= q[1] - 1.5 * iqr]), lower = q[1], middle = q[2], upper = q[3], ymax = max(v[v <= q[3] + 1.5 * iqr])) }
  bracket_rows <- function(lbls) {                  # row (1, 2) per bracket 1-2, 2-3, 1-3; NA = not drawn ("ns")
    draw <- SHOW_BRACKETS & !is.na(lbls) & lbls != "ns"; r <- rep(NA_integer_, 3)
    r[1:2][draw[1:2]] <- 1L
    if (draw[3]) r[3] <- if (any(draw[1:2])) 2L else 1L
    r }
  bracket_x <- list(c(1, 1.97), c(2.03, 3), c(1, 3))   # 1-2 and 2-3 share a row, so they stop short of each other at x = 2
  whisker_range <- function(df, grp) {              # lowest lower whisker and highest upper whisker across the ages
    b <- sapply(split(df$log2FC, df[[grp]]), box_stats); c(min(b["ymin", ]), max(b["ymax", ])) }

  make_arm_plot <- function(df_arm, direction, ct, sx, COL_LINE, COL_FILL, pv_star) {
    age_labels <- c("YEEYHC\n(Young)", "MEEMHC\n(Middle)", "OEEOHC\n(Old)")
    df_arm <- df_arm %>% dplyr::filter(!is.na(log2FC)) %>% dplyr::mutate(age_f = factor(age, levels = age_order, labels = age_labels))
    if (nrow(df_arm) == 0) return(NULL)
    n_genes <- length(unique(df_arm$gene))
    means <- df_arm %>% dplyr::group_by(age_f) %>% dplyr::summarise(mean_lfc = mean(log2FC, na.rm = TRUE), .groups = "drop")
    flip <- direction == "Down"
    if (identical(NOLINE_Y_ZOOM, "whisker")) {      # fit the boxes: whiskers + room for the brackets
      wr <- whisker_range(df_arm, "age_f")
      ymax <- if (!flip) wr[2] else -wr[1]; if (ZOOM_INCLUDE_ZERO) ymax <- max(ymax, 0); tick <- abs(ymax) * 0.025 + 1e-6
      far <- if (!flip) wr[1] else wr[2]                                    # whisker on the other side
      if (ZOOM_INCLUDE_ZERO) far <- if (!flip) min(far, 0) else max(far, 0)
      span_n <- ymax - (if (!flip) far else -far)
    } else {
      ymax <- max(quantile(abs(df_arm$log2FC), 0.95, na.rm = TRUE) * 1.4, 1); tick <- ymax * 0.025
      far <- if (!flip) min(df_arm$log2FC, na.rm = TRUE) * 1.05 else max(df_arm$log2FC, na.rm = TRUE) * 1.05
      span_n <- ymax - (if (!flip) far else -far)
    }
    lbl3 <- sapply(pv_star, sig_label); rows <- bracket_rows(lbl3)          # BH-corrected P
    n_rows <- max(c(0, rows), na.rm = TRUE)
    h_mag <- ymax + span_n * (BRACKET_FIRST + (rows - 1) * BRACKET_ROW)          # bracket heights (only the drawn ones)
    h <- (if (!flip) 1 else -1) * h_mag
    y_pad <- if (n_rows > 0) span_n * (BRACKET_FIRST + n_rows * BRACKET_ROW) else span_n * 0.002   # no brackets = no space kept
    edge_pad <- span_n * 0.002                                              # hairline so whisker tips are never clipped
    y_lim <- if (!flip) c(far - edge_pad, ymax + y_pad) else c(-ymax - y_pad, far + edge_pad)
    p <- ggplot() +
      geom_hline(yintercept = 0, linewidth = 0.6, color = "grey30", linetype = "dashed") +
      geom_boxplot(data = df_arm, aes(x = age_f, y = log2FC), fill = adjustcolor(COL_FILL, 0.20), color = COL_LINE,
                   linewidth = 0.8, outlier.shape = NA, width = 0.50) +
      geom_line(data = means, aes(x = age_f, y = mean_lfc, group = 1), color = COL_LINE, linewidth = 2.4, alpha = 0.95) +
      geom_point(data = means, aes(x = age_f, y = mean_lfc), color = COL_LINE, size = 4.5, shape = 16) +
      { if (identical(NOLINE_Y_ZOOM, "whisker")) scale_y_continuous(expand = expansion(0)) } +   # zoomed: no y margin
      coord_cartesian(ylim = y_lim, clip = "on") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            axis.text.x = element_text(size = 11), panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3)) +
      labs(title = paste0(ct, " | ", direction, " genes | ", tools::toTitleCase(sx)),
           subtitle = paste0("n=", n_genes, " YEEYHC-significant genes", "  |  Box = IQR (outliers hidden)  |  Line = mean",
                             if (identical(NOLINE_Y_ZOOM, "whisker")) "  |  y-axis: box whiskers" else "  |  y-axis: 95th pct"),
           x = NULL, y = "log2FC (EE vs HC)")
    add_bracket <- function(plt, x1, x2, y, lbl, col, flip = FALSE) {
      if (is.na(lbl) || lbl == "ns") return(plt)
      s <- if (!flip) -1 else 1
      plt + annotate("segment", x = x1, xend = x2, y = y, yend = y, color = col, linewidth = 0.6) +
        annotate("segment", x = x1, xend = x1, y = y + s * abs(tick), yend = y, color = col, linewidth = 0.6) +
        annotate("segment", x = x2, xend = x2, y = y + s * abs(tick), yend = y, color = col, linewidth = 0.6) +
        annotate("text", x = (x1 + x2) / 2, y = y - s * abs(tick) * 0.9, label = lbl, size = 3.8, color = col, fontface = "bold") }
    for (k in 1:3) if (!is.na(rows[k]))
      p <- add_bracket(p, bracket_x[[k]][1], bracket_x[[k]][2], h[k], lbl3[k], COL_LINE, flip)
    p
  }

  paired_p <- function(d, a, b) {                   # two-sided paired Wilcoxon, >= 5 paired genes
    g1 <- d[d$age == a, ]; g2 <- d[d$age == b, ]
    both <- intersect(g1$gene, g2$gene); if (length(both) < 5) return(NA_real_)
    tryCatch(wilcox.test(g1$log2FC[match(both, g1$gene)], g2$log2FC[match(both, g2$gene)], paired = TRUE)$p.value,
             error = function(e) NA_real_)
  }

  ## ---- 1. gene sets per sex x cell type x arm, and BH-corrected Wilcoxon P -----------------------------------------
  arm_data <- list(); stats_rows <- list()
  for (sx in SEXES) {
    tabs <- setNames(lapply(age_order, function(a) load_deseq(sx, a)), age_order)
    yee  <- tabs[["YEEYHC"]]
    for (ct in ct_levels_of(unique(yee$celltype))) {
      yee_ct <- yee[yee$celltype == ct, ]
      sig_up <- sig_genes(yee_ct, 1); sig_dn <- sig_genes(yee_ct, -1)
      if (length(sig_up) < MIN_GENES && length(sig_dn) < MIN_GENES) next
      for (arm in c("Up", "Down")) {
        genes <- if (arm == "Up") sig_up else sig_dn
        d <- dplyr::bind_rows(lapply(age_order, function(a) {
          s <- tabs[[a]]; s <- s[s$celltype == ct & s$gene %in% genes & !is.na(s$log2FoldChange), ]
          data.frame(gene = s$gene, log2FC = s$log2FoldChange, age = rep(a, nrow(s)), stringsAsFactors = FALSE) }))
        if (!nrow(d)) next
        arm_data[[paste(sx, ct, arm, sep = "|")]] <- d
        stats_rows[[paste(sx, ct, arm, sep = "|")]] <- data.frame(
          sex = sx, celltype = ct, arm = arm, n_genes = length(unique(d$gene)),
          test = c("Young_vs_Middle", "Middle_vs_Old", "Young_vs_Old"),
          p = c(paired_p(d, "YEEYHC", "MEEMHC"), paired_p(d, "MEEMHC", "OEEOHC"), paired_p(d, "YEEYHC", "OEEOHC")),
          stringsAsFactors = FALSE)
      }
    }
  }
  st <- dplyr::bind_rows(stats_rows) %>%
    dplyr::group_by(sex, arm) %>% dplyr::mutate(p_BH = p.adjust(p, "BH")) %>% dplyr::ungroup()
  write.csv(st, file.path(OUT_DIR, paste0(fig_date, "_Wilcoxon_Stats_BH.csv")), row.names = FALSE)
  P_BH <- lapply(split(st, paste(st$sex, st$celltype, st$arm, sep = "|")),
                 function(x) x$p_BH[match(c("Young_vs_Middle", "Middle_vs_Old", "Young_vs_Old"), x$test)])

  ## ---- 2. ED7a,b: one plot per sex x cell type x arm -----------------------------------------------------------
  for (key in names(arm_data)) {
    k <- strsplit(key, "|", fixed = TRUE)[[1]]; sx <- k[1]; ct <- k[2]; arm <- k[3]
    col <- if (arm == "Up") COL_UP else COL_DN
    p <- make_arm_plot(arm_data[[key]], arm, ct, sx, col, col, P_BH[[key]])
    if (!is.null(p)) ggsave(file.path(OUT_DIR, paste0(fig_date, "_", gsub("[^A-Za-z0-9]+", "_", ct), "_", sx, "_",
                                                      toupper(arm), "_only.pdf")), p, width = 7, height = 6, device = "pdf")
  }

  ## ---- 3. Fig. 3a: all-cell-types mean trajectory ----------------------------------------------------------------

  AGE_LABELS <- c("Young EE\n(YEEYHC)", "Middle EE\n(MEEMHC)", "Old EE\n(OEEOHC)")
  for (sx in SEXES) {
    yee2 <- load_deseq(sx, "YEEYHC"); mee2 <- load_deseq(sx, "MEEMHC"); oee2 <- load_deseq(sx, "OEEOHC")
    mean_rows <- list()
    for (ct in sort(unique(yee2$celltype))) {
      yee_ct <- yee2[yee2$celltype == ct, ]
      for (arm in c("Up", "Down")) {
        anchor <- sig_genes(yee_ct, if (arm == "Up") 1 else -1); if (length(anchor) < MIN_GENES) next
        f <- function(d) mean(d$log2FoldChange[d$celltype == ct & d$gene %in% anchor], na.rm = TRUE)
        mean_rows[[paste(ct, arm)]] <- data.frame(celltype = ct, arm = arm, age_f = factor(AGE_LABELS, levels = AGE_LABELS),
                                                  mean_lfc = c(f(yee2), f(mee2), f(oee2)), stringsAsFactors = FALSE)
      }
    }
    if (!length(mean_rows)) next
    mean_df <- dplyr::bind_rows(mean_rows)
    all_cts <- unique(mean_df$celltype)
    ct_cols <- c(CT_COLS[intersect(names(CT_COLS), all_cts)], setNames(rep("grey60", length(setdiff(all_cts, names(CT_COLS)))), setdiff(all_cts, names(CT_COLS))))
    arm_labels <- c(Up = "Upregulated in Young EE", Down = "Downregulated in Young EE")
    mean_df$arm_lbl <- factor(arm_labels[mean_df$arm], levels = arm_labels)
    p_mean <- ggplot(mean_df, aes(x = age_f, y = mean_lfc, color = celltype, group = celltype)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.7) +
      geom_line(linewidth = 0.9, alpha = 0.85) + geom_point(size = 2.5, shape = 16) +
      scale_color_manual(values = ct_cols, name = "Cell type") + facet_wrap(~arm_lbl, ncol = 1, scales = "free_y") +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5), plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
            axis.text.x = element_text(size = 10.5), strip.text = element_text(face = "bold", size = 11), legend.text = element_text(size = 8),
            panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3)) +
      guides(color = guide_legend(ncol = 1, override.aes = list(linewidth = 2, size = 3))) +
      labs(title = paste0("Mean log2FC trajectories — All cell types (", tools::toTitleCase(sx), ")"),
           subtitle = "Anchor = YEEYHC-significant genes tracked through Middle and Old EE\nOne line per cell type  |  Converging toward 0 = dynamic range compression",
           x = NULL, y = "Mean log2FC (EE vs HC)")
    ggsave(file.path(OUT_DIR, paste0(fig_date, "_AllCelltypes_MeanTrajectory_", sx, ".pdf")), p_mean, width = 9, height = 10, device = "pdf")
  }
  invisible(NULL)
}

################################################################################
### STEP compress: ED7c COMPRESSION SCORE (LFC ratio)
###   Per sex x cell type x arm (Up / Down young-NE genes, >= MIN_GENES):
###     arm ratio = mean|log2FC in YEEYHC| / mean|log2FC at the older age| over the same genes (capped at 10)
###   Cell-type score = mean of the Up and Down arm ratios (one arm if only one qualifies). > 1 = compressed.
###   Young->Middle (MEEMHC) and Young->Old (OEEOHC). 95% CI: N_BOOT resamples of genes within each arm.
###   Plots:
###     Compression_<Sex>      cell types ranked within that sex, Middle and Old points + CI (main version)
###     Compression_BothSexes  the two sex panels side by side
###     Compression_Old_SexColours / _Middle_  one age, males blue + females red on shared rows (earlier ED style)
################################################################################

SECTIONS[["compress"]] <- function() {
  OUT <- step_dir("Compression_LFCRatio")
  ratio <- function(y, a) min(mean(abs(y)) / max(mean(abs(a)), 1e-6), 10)
  rows <- list(); arm_rows <- list()
  set.seed(42)
  for (sx in SEXES) {
    yee <- load_deseq(sx, "YEEYHC")
    for (age in c("Middle", "Old")) {
      aged <- load_deseq(sx, AGE_COMPS[[age]][["ee"]])
      for (ct in unique(yee$celltype)) {
        y <- yee[yee$celltype == ct & !is.na(yee$log2FoldChange), ]; a <- aged[aged$celltype == ct & !is.na(aged$log2FoldChange), ]
        arms <- list()
        for (arm in c("Up", "Down")) {
          g <- intersect(sig_genes(y, if (arm == "Up") 1 else -1), a$gene)
          if (length(g) < MIN_GENES) next
          arms[[arm]] <- cbind(y$log2FoldChange[match(g, y$gene)], a$log2FoldChange[match(g, a$gene)])
          arm_rows[[paste(sx, age, ct, arm)]] <- data.frame(sex = sx, age = age, celltype = ct, arm = arm, n_genes = length(g),
            mean_abs_young = mean(abs(arms[[arm]][, 1])), mean_abs_aged = mean(abs(arms[[arm]][, 2])),
            lfc_ratio = ratio(arms[[arm]][, 1], arms[[arm]][, 2]))
        }
        if (!length(arms)) next
        score <- mean(sapply(arms, function(m) ratio(m[, 1], m[, 2])))
        boot <- replicate(N_BOOT, mean(sapply(arms, function(m) { i <- sample.int(nrow(m), replace = TRUE); ratio(m[i, 1], m[i, 2]) })))
        rows[[paste(sx, age, ct)]] <- data.frame(sex = sx, age = age, celltype = ct, score = score,
          ci_lo = unname(quantile(boot, 0.025)), ci_hi = unname(quantile(boot, 0.975)), p_compress = mean(boot <= 1),
          n_up = if (is.null(arms$Up)) 0L else nrow(arms$Up), n_down = if (is.null(arms$Down)) 0L else nrow(arms$Down),
          arms_used = paste(names(arms), collapse = "+"))
      }
    }
  }
  sc <- dplyr::bind_rows(rows)
  if (!nrow(sc)) { say("  no cell type had enough young-NE genes"); return(invisible(NULL)) }
  sc <- sc %>% dplyr::group_by(sex, age) %>% dplyr::mutate(fdr_compress = p.adjust(p_compress, "BH")) %>% dplyr::ungroup()
  write.csv(sc, file.path(OUT, "Compression_LFCRatio_Scores.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(arm_rows), file.path(OUT, "Compression_LFCRatio_ByArm.csv"), row.names = FALSE)
  XLAB <- "Compression score: LFC ratio only (mean|LFC| Y/A)"

  one_sex <- function(sx) {
    d <- sc %>% dplyr::filter(sex == sx) %>% dplyr::mutate(age = factor(age, levels = c("Middle", "Old")))
    ord <- d %>% dplyr::filter(age == "Old") %>% dplyr::arrange(score) %>% dplyr::pull(celltype)
    ord <- c(setdiff(unique(d$celltype), ord), ord)          # cell types without an old score go to the bottom
    d$celltype <- factor(d$celltype, levels = ord)
    seg <- d %>% dplyr::select(celltype, age, score) %>% tidyr::pivot_wider(names_from = age, values_from = score)
    ggplot(d, aes(y = celltype)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
      { if (all(c("Middle", "Old") %in% names(seg))) geom_segment(data = seg, aes(x = Middle, xend = Old, yend = celltype), color = "grey75", linewidth = 0.6) } +
      geom_linerange(aes(xmin = ci_lo, xmax = ci_hi, color = age), linewidth = 1.6, alpha = 0.35) +
      geom_point(aes(x = score, color = age), size = 2.6) +
      scale_color_manual(values = AGE_COLS, name = NULL, labels = c(Middle = paste("Young", ARROW, "Middle"), Old = paste("Young", ARROW, "Old"))) +
      theme_classic(base_size = 10) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text = element_text(color = "black"),
            legend.position = "bottom", panel.grid.major.y = element_line(color = "grey94", linewidth = 0.3)) +
      labs(title = tools::toTitleCase(sx), x = XLAB, y = NULL,
           subtitle = NULL)
  }
  dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"
  n_ct <- length(unique(sc$celltype)); H <- max(4, n_ct * 0.27 + 1.8)
  plots <- lapply(SEXES, one_sex); names(plots) <- SEXES
  for (sx in SEXES) save_both(plots[[sx]], file.path(OUT, paste0("Compression_", tools::toTitleCase(sx))), 4.6, H, dev)
  if (requireNamespace("patchwork", quietly = TRUE))
    save_both(patchwork::wrap_plots(plots, nrow = 1, guides = "collect") & theme(legend.position = "bottom"),
              file.path(OUT, "Compression_BothSexes"), 8.8, H, dev)

  for (age in c("Old", "Middle")) {        # earlier style: one age, sexes coloured on shared rows (side by side, no test)
    d <- sc %>% dplyr::filter(age == !!age)
    ord <- d %>% dplyr::filter(sex == "female") %>% dplyr::arrange(score) %>% dplyr::pull(celltype)
    d$celltype <- factor(d$celltype, levels = c(setdiff(unique(d$celltype), ord), ord))
    p <- ggplot(d, aes(x = score, y = celltype, color = sex)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
      geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.5, alpha = 0.55, position = position_dodge(width = 0.5)) +
      geom_point(size = 2.6, position = position_dodge(width = 0.5)) +
      scale_color_manual(values = SEX_COLS, breaks = SEXES, labels = c("Male", "Female"), name = NULL) +
      theme_classic(base_size = 10) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text = element_text(color = "black"), legend.position = "bottom",
            panel.grid.major.y = element_line(color = "grey94", linewidth = 0.3)) +
      labs(title = paste("Young", ARROW, age), x = XLAB, y = NULL)
    save_both(p, file.path(OUT, paste0("Compression_", age, "_SexColours")), 4.8, H, dev)
  }
  invisible(NULL)
}

################################################################################
### RUN
################################################################################

status <- data.frame(section = character(0), result = character(0), minutes = numeric(0), stringsAsFactors = FALSE)
for (id in intersect(names(SECTION_INFO), RUN)) {
  banner(id)
  t0 <- Sys.time()
  res <- tryCatch({ SECTIONS[[id]](); "done" },
                  error = function(e) { say("  ERROR in step ", id, ": ", conditionMessage(e)); paste("FAILED:", conditionMessage(e)) })
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  invisible(gc())
  status <- rbind(status, data.frame(section = id, result = res, minutes = mins))
  say(sprintf("STEP %s %s in %.1f min", id, if (res == "done") "finished" else "FAILED", mins))
}

cat("\n", strrep("=", 90), "\n", sep = "")
say("Figure 3 pipeline finished. Summary:")
for (i in seq_len(nrow(status)))
  say(sprintf("  %-9s %-5s %5.1f min  %s", status$section[i], ifelse(status$result[i] == "done", "ok", "FAIL"),
              status$minutes[i], SECTION_INFO[[status$section[i]]]))
if (any(status$result != "done")) say("Failed steps can be rerun alone, e.g. RUN <- c(\"compress\").")
say("Outputs: ", FIG_ROOT, " (one subfolder per panel). Log: ", LOG)
