################################################################################
### Fig3b_Boxplots_BH.R
###
### Standalone Fig 3b trajectory boxplots (grey-line version), with significance stars from
### BH-corrected paired Wilcoxon tests.
###   Anchor  = YEEYHC-significant genes (Up / Down separately), tracked to MEEMHC and OEEOHC
###   Test    = two-sided paired Wilcoxon signed-rank on gene log2FC between ages (>= 5 paired genes)
###   Stars   = BH-corrected P, corrected across all cell types x the 3 age comparisons within each sex x arm
###             (* < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001; non-significant brackets are not drawn)
### Fig. 3b uses the "fit" version (FitWindow/): y axis fitted to the box whiskers; grey lines drawn only for genes whose
### whole trajectory lies inside that window (boxes and means use all genes). Panels in Fig. 3b: male DG ExNs, female
### CA3 ExNs, female astrocytes (PANEL_3B). Other axis versions can be added to GREY_MODES.
### Reads the per-sex DESeq2 tables written by 01_DESeq2/Fig2_ED4_Pipeline.R; DESeq2 is not rerun. Overwrites its outputs.
### Run:  Rscript Fig3b_Boxplots_BH.R     (set BASE below first)
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(conflicted)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::mutate, .quiet=TRUE); conflicts_prefer(base::intersect, .quiet=TRUE)
  conflicts_prefer(base::setdiff, .quiet=TRUE)
})

## ---- settings ------------------------------------------------------------------------------------------------------
BASE      <- "PATH_TO_OUTPUT"                                  # same BASE as 01_DESeq2/Fig2_ED4_Pipeline.R
DESEQ_DIR <- file.path(BASE, "DESeq2_Results")                   # <sex>_ALL_<comparison>.csv
OUT_DIR   <- file.path(BASE, "Paper_Figures", "Figure3", "3b_Trajectory_Boxplots", "GreyLines_BH")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SEXES     <- c("male", "female")
PADJ_THR  <- 0.05      # anchor genes: BH padj < 0.05 in YEEYHC
MIN_GENES <- 10        # a cell type is plotted if either arm has >= 10 anchor genes
# Axis versions, each in its own folder: "fit" -> FitWindow/ (Fig. 3b); also available:
#   "full" -> Full/ (axis spans every gene line), "zoom" -> Zoom/ (axis at the whiskers, extreme lines run off),
#   "break" -> BrokenAxis/ (linear to the whiskers, then broken and compressed so every line stays visible)
GREY_MODES <- c("fit")
GREY_Y_QUANTILE   <- 1     # "full": 1 = no ceiling
GREY_BREAK_SHRINK <- 0.15  # "break": compression beyond the break
SHOW_BRACKETS     <- TRUE
ZOOM_INCLUDE_ZERO <- TRUE
BRACKET_FIRST <- 0.03
BRACKET_ROW   <- 0.07
SEX_SYMBOL    <- TRUE
SYMBOL_FONT   <- "DejaVu Sans"
PANEL_3B <- list(c("DG ExNs", "male"), c("CA3 ExNs", "female"), c("Astrocytes", "female"))
COL_UP <- "#D6604D"; COL_DN <- "#4393C3"
age_order <- c("YEEYHC", "MEEMHC", "OEEOHC")
CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN",
              "Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal",
              "Cajal-Retzius Cell","Macrophage","Epithelial","T Cell")

## ---- helpers -------------------------------------------------------------------------------------------------------
load_deseq <- function(sx, comp) {
  f <- file.path(DESEQ_DIR, paste0(sx, "_ALL_", comp, ".csv"))
  if (!file.exists(f)) stop("Missing DESeq2 table: ", f)
  read.csv(f, stringsAsFactors = FALSE)
}
ct_levels_of <- function(cts) c(intersect(CT_ORDER, cts), sort(setdiff(cts, CT_ORDER)))
sig_genes <- function(d, dir = 0) {
  ok <- !is.na(d$padj) & d$padj < PADJ_THR
  if (dir > 0) ok <- ok & d$log2FoldChange > 0
  if (dir < 0) ok <- ok & d$log2FoldChange < 0
  unique(d$gene[ok])
}
sig_label <- function(p) {
  if (is.na(p)) return("ns"); if (p < 1e-4) return("****"); if (p < 1e-3) return("***")
  if (p < 0.01) return("**"); if (p < 0.05) return("*"); "ns"
}
save_both <- function(p, stem, w, h, device_pdf = "pdf") {      # PDF + 600 dpi PNG
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = w, height = h, device = device_pdf, limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, ".png"), p, width = w, height = h, device = "png", dpi = 600, bg = "white", limitsize = FALSE)
  cat("  saved:", basename(stem), "\n")
}
pull_arm <- function(tabs, ct, genes) {                           # log2FC of the anchor genes at each age
  dplyr::bind_rows(lapply(age_order, function(a) {
    s <- tabs[[a]]; s <- s[s$celltype == ct & s$gene %in% genes & !is.na(s$log2FoldChange), ]
    data.frame(gene = s$gene, log2FC = s$log2FoldChange, age = rep(a, nrow(s)), stringsAsFactors = FALSE)
  }))
}
paired_p <- function(d, a, b) {                                   # two-sided paired Wilcoxon, >= 5 paired genes
  g1 <- d[d$age == a, ]; g2 <- d[d$age == b, ]
  both <- intersect(g1$gene, g2$gene); if (length(both) < 5) return(NA_real_)
  tryCatch(wilcox.test(g1$log2FC[match(both, g1$gene)], g2$log2FC[match(both, g2$gene)], paired = TRUE)$p.value,
           error = function(e) NA_real_)
}

## ---- box / bracket helpers -----------------------------------------------------------
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


## ---- grey-line arm plot (BH-corrected stars) ------------------------------------------
  use_cairo <- SEX_SYMBOL && isTRUE(capabilities("cairo")) && isTRUE(l10n_info()[["UTF-8"]])   # symbols need a UTF-8 session
  pdf_dev   <- if (use_cairo) grDevices::cairo_pdf else "pdf"
  sex_tag   <- function(sx) if (use_cairo) c(male = "\u2642", female = "\u2640")[[sx]] else tools::toTitleCase(sx)
  make_arm_plot_lines <- function(df_arm, direction, ct, sx, COL, mode = "full") {
    labs3 <- c("Young", "Middle", "Old")
    df_arm <- df_arm %>% dplyr::filter(!is.na(log2FC)) %>%
      dplyr::mutate(age_f = factor(age, levels = age_order, labels = labs3), x = as.integer(age_f))
    if (nrow(df_arm) == 0) return(NULL)
    means <- df_arm %>% dplyr::group_by(x) %>% dplyr::summarise(m = mean(log2FC), .groups = "drop")
    pw <- function(a, b) { g1 <- df_arm[df_arm$x == a, ]; g2 <- df_arm[df_arm$x == b, ]
      both <- intersect(g1$gene, g2$gene); if (length(both) < 5) return(NA_real_)
      tryCatch(wilcox.test(g1$log2FC[match(both, g1$gene)], g2$log2FC[match(both, g2$gene)], paired = TRUE)$p.value,
               error = function(e) NA_real_) }
    pv <- c(pw(1, 2), pw(2, 3), pw(1, 3))
    q  <- GREY_Y_QUANTILE
    lo <- min(0, quantile(df_arm$log2FC, 1 - q)); hi <- max(0, quantile(df_arm$log2FC, q))
    up <- direction == "Up"; s <- if (up) 1 else -1
    wr <- whisker_range(df_arm, "x")
    # axis transform: identity, except in "break" mode (linear out to the whiskers, compressed beyond)
    fwd <- function(y) y; inv <- function(y) y; b_lo <- NA_real_; b_hi <- NA_real_
    if (mode %in% c("zoom", "fit")) { lo <- if (ZOOM_INCLUDE_ZERO) min(0, wr[1]) else wr[1]; hi <- if (ZOOM_INCLUDE_ZERO) max(0, wr[2]) else wr[2] }
    # "fit": grey lines only for genes whose whole trajectory sits inside the window (boxes / means keep all genes)
    n_all <- length(unique(df_arm$gene))
    line_df <- if (mode == "fit") {
      keep <- df_arm %>% dplyr::group_by(gene) %>%
        dplyr::summarise(inside = dplyr::n() == 3 && all(log2FC >= lo & log2FC <= hi), .groups = "drop") %>%
        dplyr::filter(inside) %>% dplyr::pull(gene)
      df_arm[df_arm$gene %in% keep, ]
    } else df_arm
    if (mode == "break") {
      k_br <- GREY_BREAK_SHRINK   # (not 'k': the bracket loop below reuses k)
      b_hi <- max(0, wr[2]); b_lo <- min(0, wr[1])                     # break right at the whisker ends
      if (hi <= b_hi) b_hi <- NA_real_; if (lo >= b_lo) b_lo <- NA_real_       # break only on a side that needs it
      fwd <- function(y) { o <- y
        if (is.finite(b_hi)) o <- ifelse(y > b_hi, b_hi + (y - b_hi) * k_br, o)
        if (is.finite(b_lo)) o <- ifelse(y < b_lo, b_lo + (y - b_lo) * k_br, o); o }
      inv <- function(y) { o <- y
        if (is.finite(b_hi)) o <- ifelse(y > b_hi, b_hi + (y - b_hi) / k_br, o)
        if (is.finite(b_lo)) o <- ifelse(y < b_lo, b_lo + (y - b_lo) / k_br, o); o }
    }
    # brackets and limits are laid out on the drawn (transformed) scale, then mapped back to log2FC
    t_lo <- fwd(lo); t_hi <- fwd(hi); span <- t_hi - t_lo
    edge <- if (up) t_hi else t_lo
    pv_star <- P_BH[[paste(sx, ct, direction, sep = "|")]]; if (is.null(pv_star)) pv_star <- pv   # BH-corrected stars
    lbl3 <- sapply(pv_star, sig_label); rows <- bracket_rows(lbl3); n_rows <- max(c(0, rows), na.rm = TRUE)
    yb_t <- edge + s * span * (BRACKET_FIRST + (rows - 1) * BRACKET_ROW)   # bracket heights (drawn brackets only)
    tick_t <- span * 0.02
    yb <- inv(yb_t); yb_tick <- inv(yb_t - s * tick_t); yb_lab <- inv(yb_t + s * tick_t * 1.2)
    room <- if (n_rows > 0) span * (BRACKET_FIRST + n_rows * BRACKET_ROW) else span * 0.002   # no brackets = no space kept
    hair <- span * 0.002                                                   # hairline so whisker tips are never clipped
    ylim <- inv(if (up) c(t_lo - hair, t_hi + room) else c(t_lo - room, t_hi + hair))   # coord_cartesian takes log2FC units
    box_layer <- if (mode == "break") {                  # stats fixed in log2FC units so the transform cannot alter whiskers
      bs <- do.call(rbind, lapply(split(df_arm$log2FC, df_arm$x), box_stats)); bs <- data.frame(x = as.integer(rownames(bs)), bs)
      geom_boxplot(data = bs, aes(x = x, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, group = x),
                   stat = "identity", fill = adjustcolor(COL, 0.30), color = COL, linewidth = 0.7, width = 0.45)
    } else geom_boxplot(data = df_arm, aes(x = x, y = log2FC, group = x), fill = adjustcolor(COL, 0.30), color = COL,
                        linewidth = 0.7, outlier.shape = NA, width = 0.45)
    y_scale <- if (mode == "break" && (is.finite(b_lo) || is.finite(b_hi))) {
      inner <- scales::pretty_breaks(n = 5)(c(if (is.finite(b_lo)) b_lo else lo, if (is.finite(b_hi)) b_hi else hi))
      inner <- inner[inner >= (if (is.finite(b_lo)) b_lo else -Inf) & inner <= (if (is.finite(b_hi)) b_hi else Inf)]
      outer <- c(if (is.finite(b_lo)) ceiling(lo), if (is.finite(b_hi)) floor(hi))
      outer <- outer[(outer > b_hi & is.finite(b_hi)) | (outer < b_lo & is.finite(b_lo))]
      scale_y_continuous(trans = scales::trans_new("axis_break", fwd, inv), breaks = sort(unique(c(inner, outer))),
                         labels = function(x) as.character(x), expand = expansion(0))
    } else if (mode != "full") scale_y_continuous(expand = expansion(0)) else NULL      # zoom / break: no y margin
    break_marks <- if (mode == "break") {
      bb <- c(b_lo, b_hi); bb <- bb[is.finite(bb)]
      if (length(bb)) list(geom_hline(yintercept = bb, linetype = "dotted", colour = "grey55", linewidth = 0.3),
                           annotate("text", x = -Inf, y = bb, label = "//", size = 3.2, colour = "grey20"))
    } else NULL
    p <- ggplot() +
      geom_hline(yintercept = 0, linewidth = 0.6, color = "grey30", linetype = "dashed") +
      geom_line(data = line_df, aes(x = x, y = log2FC, group = gene), color = "grey70", alpha = 0.35, linewidth = 0.25) +
      box_layer + break_marks +
      geom_line(data = means, aes(x = x, y = m), color = COL, linewidth = 1.0) +
      geom_point(data = means, aes(x = x, y = m), color = COL, size = 2.4) +
      scale_x_continuous(breaks = 1:3, labels = labs3, expand = expansion(add = 0.45)) + y_scale +
      coord_cartesian(ylim = ylim, clip = if (mode %in% c("zoom", "fit")) "on" else "off") +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5, family = if (use_cairo) SYMBOL_FONT else ""),
            axis.text = element_text(color = "black"),
            plot.subtitle = element_text(size = 7.5, color = "grey40", hjust = 0.5)) +
      labs(title = paste(sex_tag(sx), ct),
           subtitle = if (mode == "fit") paste0(direction, " in young NE, n = ", n_all, " genes (", length(unique(line_df$gene)), " lines shown)")
                      else paste0(direction, " in young NE, n = ", n_all, " genes"),
           x = NULL, y = expression(log[2]*FC~"(NE vs HC)"))
    for (k in 1:3) { lbl <- lbl3[k]; if (is.na(rows[k])) next
      a <- bracket_x[[k]][1]; b <- bracket_x[[k]][2]
      p <- p + annotate("segment", x = a, xend = b, y = yb[k], yend = yb[k], linewidth = 0.5) +
        annotate("segment", x = c(a, b), xend = c(a, b), y = yb[k], yend = yb_tick[k], linewidth = 0.5) +
        annotate("text", x = (a + b) / 2, y = yb_lab[k], label = lbl, size = 3.4, vjust = if (up) 0.3 else 0.9) }
    attr(p, "pvals") <- pv
    p
  }

## ---- 1. collect gene sets and BH-corrected Wilcoxon P ---------------------------------------------------------------
arm_data <- list(); stats_rows <- list()
for (sx in SEXES) {
  tabs <- setNames(lapply(age_order, function(a) load_deseq(sx, a)), age_order)
  yee  <- tabs[["YEEYHC"]]
  for (ct in ct_levels_of(unique(yee$celltype))) {
    yee_ct <- yee[yee$celltype == ct, ]
    sig_up <- sig_genes(yee_ct, 1); sig_dn <- sig_genes(yee_ct, -1)
    if (length(sig_up) < MIN_GENES && length(sig_dn) < MIN_GENES) next
    for (arm in c("Up", "Down")) {
      d <- pull_arm(tabs, ct, if (arm == "Up") sig_up else sig_dn)
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
utils::write.csv(st, file.path(OUT_DIR, "Fig3b_Wilcoxon_Stats_BH.csv"), row.names = FALSE)
P_BH <- lapply(split(st, paste(st$sex, st$celltype, st$arm, sep = "|")),
               function(x) x$p_BH[match(c("Young_vs_Middle", "Middle_vs_Old", "Young_vs_Old"), x$test)])
cat("Wilcoxon tests:", nrow(st), "| significant after BH:", sum(st$p_BH < 0.05, na.rm = TRUE), "\n\n")

## ---- 2. draw every cell type x sex x arm, in each axis mode ---------------------------------------------------------
mode_dir <- function(md) { d <- file.path(OUT_DIR, c(full = "Full", zoom = "Zoom", "break" = "BrokenAxis", fit = "FitWindow")[[md]])
  dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
store <- list()
for (key in names(arm_data)) {
  k <- strsplit(key, "|", fixed = TRUE)[[1]]; sx <- k[1]; ct <- k[2]; arm <- k[3]
  col <- if (arm == "Up") COL_UP else COL_DN
  slug <- gsub("[^A-Za-z0-9]+", "_", ct)
  for (md in GREY_MODES) {
    g <- make_arm_plot_lines(arm_data[[key]], arm, ct, sx, col, mode = md)
    if (is.null(g)) next
    store[[paste(ct, sx, arm, md, sep = "|")]] <- g
    save_both(g, file.path(mode_dir(md), paste0("Fig3b_", slug, "_", sx, "_", toupper(arm), "_only")), 3.2, 3.6, pdf_dev)
  }
}

## ---- 3. assembled main panel (Up row on top, Down row below) -------------------------------------------------------
if (requireNamespace("patchwork", quietly = TRUE)) {
  for (md in GREY_MODES) {
    pick <- function(arm) lapply(PANEL_3B, function(e) store[[paste(e[1], e[2], arm, md, sep = "|")]])
    up <- pick("Up"); dn <- pick("Down")
    if (any(sapply(c(up, dn), is.null))) { cat("Main panel (", md, "): a requested cell type / sex has no plot; skipped\n"); next }
    strip <- function(ps, ttl) patchwork::wrap_plots(ps, nrow = 1) +
      patchwork::plot_annotation(title = ttl, theme = theme(plot.title = element_text(hjust = 0.5, size = 13)))
    panel <- patchwork::wrap_elements(strip(up, "Upregulated in Young NE")) /
             patchwork::wrap_elements(strip(dn, "Downregulated in Young NE"))
    save_both(panel, file.path(mode_dir(md), "Fig3b_MainPanel"), 3.0 * length(PANEL_3B) + 0.6, 7.6, pdf_dev)
  }
}
cat("\nDone. Output:", OUT_DIR, "\n")
