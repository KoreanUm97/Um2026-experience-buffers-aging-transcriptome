################################################################################
### AgeEmergent_ARG_StackedBar.R
###
### Mirrored stacked bar: age-emergent ARGs (novel EE genes). Reads DESeq2_Results.
###   Genes significant in MEEMHC or OEEOHC but NOT in YEEYHC.
###
### Layout: Male bars LEFT, Female bars RIGHT
###         Segments: Novel Up (dark blue) | Novel Down (light blue)
###
### Outputs
###   1. _ED8a_Mirrored.pdf            — combined (MEEMHC + OEEOHC, shared scale)
###   2. _ED8a_Mirrored_MEEMHC.pdf     — middle-age only, scaled to its own max
###   3. _ED8a_Mirrored_OEEOHC.pdf     — old-age only,    scaled to its own max
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

suppressPackageStartupMessages({
  library(conflicted)
  conflicts_prefer(dplyr::select,   .quiet = TRUE)
  conflicts_prefer(dplyr::filter,   .quiet = TRUE)
  conflicts_prefer(base::intersect, .quiet = TRUE)
  conflicts_prefer(base::setdiff,   .quiet = TRUE)
  conflicts_prefer(base::union,     .quiet = TRUE)
})

################################################################################
### PATHS
################################################################################

base_path <- "PATH_TO_NE_OUTPUT/"
deseq_dir <- paste0(base_path, "DESeq2_Results/")
fig_date  <- format(Sys.Date(), "%Y-%m-%d")
out_dir   <- paste0(base_path, "Paper_Figures/Figure5/ED8a_AgeEmergent_ARG_StackedBar/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output:", out_dir, "\n\n")

PADJ_THR <- 0.05

CT_ORDER <- c(
  "DG ExNs", "CA3 ExNs", "CA1 ExNs", "CA2 ExNs", "Sub ExNs",
  "IN",
  "Astrocytes", "Oligodendrocytes", "OPC", "Microglia",
  "Endothelial", "Pericyte", "APC", "Stromal",
  "Cajal-Retzius Cell", "Macrophage", "Epithelial", "T Cell"
)

SEG_ORDER <- c("Novel Up", "Novel Down")
SEG_COLS  <- c(
  "Novel Up"   = "#2171B5",
  "Novel Down" = "#6BAED6"
)

should_skip <- function(f) {
  if (file.exists(f)) { cat("  Skipping:", basename(f), "\n"); return(TRUE) }
  FALSE
}

################################################################################
### CLASSIFY
################################################################################

classify_ct <- function(yee_df, comp_df, ct) {
  y     <- yee_df  %>% dplyr::filter(celltype == ct)
  d     <- comp_df %>% dplyr::filter(celltype == ct)
  y_sig <- y$gene[!is.na(y$padj) & y$padj < PADJ_THR]
  d_up  <- d$gene[!is.na(d$padj) & d$padj < PADJ_THR & d$log2FoldChange > 0]
  d_dn  <- d$gene[!is.na(d$padj) & d$padj < PADJ_THR & d$log2FoldChange < 0]
  data.frame(
    segment = SEG_ORDER,
    n_genes = c(length(setdiff(d_up, y_sig)),
                length(setdiff(d_dn, y_sig))),
    stringsAsFactors = FALSE
  )
}

################################################################################
### BUILD DATA
################################################################################

all_rows <- list()

for (sx in c("male", "female")) {
  cat(sprintf("\n### %s ###\n", toupper(sx)))
  load_df <- function(comp) {
    f <- file.path(deseq_dir, paste0(sx, "_ALL_", comp, ".csv"))
    if (!file.exists(f)) { cat("  Not found:", comp, "\n"); return(NULL) }
    read.csv(f, stringsAsFactors = FALSE)
  }
  yee <- load_df("YEEYHC")
  mee <- load_df("MEEMHC")
  oee <- load_df("OEEOHC")
  if (is.null(yee) || is.null(mee) || is.null(oee)) next

  for (ct in sort(unique(yee$celltype))) {
    for (comp_nm in c("MEEMHC", "OEEOHC")) {
      comp_df <- if (comp_nm == "MEEMHC") mee else oee
      cl <- classify_ct(yee, comp_df, ct)
      cl$celltype   <- ct
      cl$comparison <- comp_nm
      cl$sex        <- sx
      all_rows[[paste(sx, ct, comp_nm)]] <- cl
      cat(sprintf("  %-22s %-8s  up=%d  dn=%d\n",
                  ct, comp_nm, cl$n_genes[1], cl$n_genes[2]))
    }
  }
}

plot_df <- dplyr::bind_rows(all_rows) %>%
  dplyr::mutate(
    segment    = factor(segment, levels = SEG_ORDER),
    comparison = factor(comparison,
                        levels = c("MEEMHC", "OEEOHC"),
                        labels = c("Middle EE (MEEMHC)", "Old EE (OEEOHC)"))
  )

write.csv(plot_df,
          file.path(out_dir, paste0(fig_date, "_ED8a_data.csv")),
          row.names = FALSE)

################################################################################
### PLOT FUNCTION
###
### df        — already-filtered subset of plot_df (one or both comparisons)
### ct_levels — ordered cell-type factor levels (reversed for horizontal bars)
### title     — plot title string
### subtitle  — plot subtitle string
### facet     — logical; TRUE = facet_wrap(~comparison), FALSE = single panel
### x_lim     — passed-in x limit (so caller controls shared vs. independent)
################################################################################

make_mirrored_plot <- function(df, ct_levels, title, subtitle,
                               facet = TRUE, x_lim = NULL) {

  rect_df <- df %>%
    dplyr::filter(n_genes > 0) %>%
    dplyr::group_by(sex, celltype, comparison) %>%
    dplyr::arrange(segment, .by_group = TRUE) %>%
    dplyr::mutate(
      abs_end   = cumsum(n_genes),
      abs_start = abs_end - n_genes,
      xmin = ifelse(sex == "male", -abs_end,   abs_start),
      xmax = ifelse(sex == "male", -abs_start,  abs_end)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      celltype   = factor(celltype, levels = ct_levels),
      comparison = factor(as.character(comparison),
                          levels = c("Middle EE (MEEMHC)", "Old EE (OEEOHC)")),
      segment    = factor(segment, levels = SEG_ORDER),
      y_ctr      = as.integer(celltype)
    )

  totals <- rect_df %>%
    dplyr::group_by(sex, celltype, comparison, y_ctr) %>%
    dplyr::summarise(
      total = sum(n_genes),
      tip   = ifelse(sex[1] == "male", -max(abs_end), max(abs_end)),
      .groups = "drop"
    )

  x_max <- max(totals$total, na.rm = TRUE)
  if (is.null(x_lim)) x_lim <- x_max * 1.25
  BAR_H <- 0.38

  legend_pts <- data.frame(
    x       = rep(0, length(SEG_ORDER)),
    y       = rep(0, length(SEG_ORDER)),
    segment = factor(SEG_ORDER, levels = SEG_ORDER)
  )

  p <- ggplot() +
    geom_rect(
      data = rect_df,
      aes(xmin = xmin, xmax = xmax,
          ymin = y_ctr - BAR_H, ymax = y_ctr + BAR_H,
          fill = segment),
      color = "white", linewidth = 0.2
    ) +
    geom_text(
      data    = totals,
      mapping = aes(
        x     = tip + ifelse(sex == "male", -x_max * 0.02, x_max * 0.02),
        y     = y_ctr,
        label = total,
        hjust = ifelse(sex == "male", 1, 0)
      ),
      size = 2.6, color = "black", fontface = "bold"
    ) +
    geom_point(
      data  = legend_pts,
      aes(x = x, y = y, fill = segment),
      shape = 22, size = 0, color = NA
    ) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey20") +
    annotate("text", x = -x_lim * 0.55,
             y = length(ct_levels) + 0.75,
             label = "Male", fontface = "bold",
             size = 3.5, color = "grey20") +
    annotate("text", x =  x_lim * 0.55,
             y = length(ct_levels) + 0.75,
             label = "Female", fontface = "bold",
             size = 3.5, color = "grey20") +
    scale_fill_manual(
      values = SEG_COLS,
      name   = NULL,
      breaks = SEG_ORDER,
      labels = c("Novel Up"   = "Novel Up (induced in older only)",
                 "Novel Down" = "Novel Down (repressed in older only)"),
      guide  = guide_legend(override.aes = list(size = 4, shape = 22, stroke = 0))
    ) +
    scale_x_continuous(
      limits = c(-x_lim, x_lim),
      breaks = pretty(c(-x_max, x_max), n = 6),
      labels = function(x) abs(x),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = seq_along(ct_levels),
      labels = ct_levels,
      expand = expansion(add = c(0.6, 1.0))
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title         = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle      = element_text(size = 8, color = "grey45", hjust = 0.5),
      axis.text.y        = element_text(size = 9),
      axis.text.x        = element_text(size = 8),
      axis.title.x       = element_text(size = 9),
      axis.line.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      strip.text         = element_text(face = "bold", size = 10, color = "white"),
      strip.background   = element_rect(fill = "grey35", color = NA),
      legend.position    = "bottom",
      legend.text        = element_text(size = 9),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
      panel.spacing      = unit(0.6, "lines")
    ) +
    labs(title = title, subtitle = subtitle, x = "Number of novel ARGs", y = NULL)

  if (facet) {
    p <- p + facet_wrap(~comparison, ncol = 1)
  }

  p
}

################################################################################
### SHARED CELL-TYPE ORDERING (used by all three plots)
################################################################################

ct_p      <- intersect(CT_ORDER, unique(plot_df$celltype))
ct_e      <- sort(setdiff(unique(plot_df$celltype), CT_ORDER))
ct_levels <- rev(c(ct_p, ct_e))
n_cts     <- length(ct_levels)

################################################################################
### 1. COMBINED (original output — MEEMHC + OEEOHC, shared x scale)
################################################################################

out_combined <- file.path(out_dir, paste0(fig_date, "_ED8a_Mirrored.pdf"))
if (!should_skip(out_combined)) {

  # Shared x limit = max across both comparisons
  x_max_combined <- plot_df %>%
    dplyr::filter(n_genes > 0) %>%
    dplyr::group_by(sex, celltype, comparison) %>%
    dplyr::summarise(bar_total = sum(n_genes), .groups = "drop") %>%
    dplyr::pull(bar_total) %>%
    max(na.rm = TRUE)
  x_lim_combined <- x_max_combined * 1.25

  p_combined <- make_mirrored_plot(
    df        = plot_df,
    ct_levels = ct_levels,
    title     = "Age-emergent ARGs — Male vs Female",
    subtitle  = paste0(
      "\u2190 Male | Female \u2192     ",
      "Novel = significant in comparison but NOT in YEEYHC (padj<0.05)\n",
      "Numbers = total novel ARGs per bar"),
    facet  = TRUE,
    x_lim  = x_lim_combined
  )

  ggsave(out_combined, p_combined,
         width  = 10,
         height = max(8, n_cts * 0.45 + 5),
         device = "pdf")
  cat("Saved:", basename(out_combined), "\n")
}

################################################################################
### 2. MIDDLE-AGE ONLY (MEEMHC, own x scale)
################################################################################

out_mee <- file.path(out_dir, paste0(fig_date, "_ED8a_Mirrored_MEEMHC.pdf"))
if (!should_skip(out_mee)) {

  df_mee <- plot_df %>%
    dplyr::filter(as.character(comparison) == "Middle EE (MEEMHC)")

  p_mee <- make_mirrored_plot(
    df        = df_mee,
    ct_levels = ct_levels,
    title     = "Age-emergent ARGs — Middle Age (MEEMHC)",
    subtitle  = paste0(
      "\u2190 Male | Female \u2192     ",
      "Novel = significant in MEEMHC but NOT in YEEYHC (padj<0.05)\n",
      "Numbers = total novel ARGs per bar  |  x-axis scaled to MEEMHC"),
    facet  = FALSE,   # single comparison — no facet needed
    x_lim  = NULL     # auto-scale to MEEMHC data
  )

  ggsave(out_mee, p_mee,
         width  = 10,
         height = max(6, n_cts * 0.35 + 3),
         device = "pdf")
  cat("Saved:", basename(out_mee), "\n")
}

################################################################################
### 3. OLD-AGE ONLY (OEEOHC, own x scale)
################################################################################

out_oee <- file.path(out_dir, paste0(fig_date, "_ED8a_Mirrored_OEEOHC.pdf"))
if (!should_skip(out_oee)) {

  df_oee <- plot_df %>%
    dplyr::filter(as.character(comparison) == "Old EE (OEEOHC)")

  p_oee <- make_mirrored_plot(
    df        = df_oee,
    ct_levels = ct_levels,
    title     = "Age-emergent ARGs — Old Age (OEEOHC)",
    subtitle  = paste0(
      "\u2190 Male | Female \u2192     ",
      "Novel = significant in OEEOHC but NOT in YEEYHC (padj<0.05)\n",
      "Numbers = total novel ARGs per bar  |  x-axis scaled to OEEOHC"),
    facet  = FALSE,   # single comparison — no facet needed
    x_lim  = NULL     # auto-scale to OEEOHC data
  )

  ggsave(out_oee, p_oee,
         width  = 10,
         height = max(6, n_cts * 0.35 + 3),
         device = "pdf")
  cat("Saved:", basename(out_oee), "\n")
}

################################################################################
### SUMMARY
################################################################################

summ <- plot_df %>%
  dplyr::group_by(sex, comparison, segment) %>%
  dplyr::summarise(total = sum(n_genes), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = segment, values_from = total, values_fill = 0)
cat("\n=== Novel ARG totals ===\n")
print(summ)
write.csv(summ,
          file.path(out_dir, paste0(fig_date, "_ED8a_Summary.csv")),
          row.names = FALSE)

cat("\n=== DONE ===\nOutputs:", out_dir, "\n")
