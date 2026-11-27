################################################################################
### ED2_QC_Panel.R
###
### Extended Data Fig. 2: QC of the eNE dataset.
### Needs the integrated eNE Seurat object loaded as df_harmony (deposited in GEO);
### panel e reads the DESeq2 tables written by Fig2_ED4_Pipeline.R (step deseq2).
###
### PANELS:
###   a  — UMAP split by condition (12: sex×age×EE/HC) + per sample + batch
###   b  — Condition proportion by cell type (12 conditions averaged)
###   c  — Sample proportion by cell type (all individual samples)
###   d  — Cell type count by sample (stacked bar)
###   e  — DE gene counts across all pairwise comparisons per CT
###   f  — RNA QC violins per sample (nFeature, nCount, %mito)
###
### METADATA KEYS:
###   age:       "age"                    → "young","Middleage","old"
###   condition: "experimental.condition" → "EE","HC"
###   sex:       "sex"                    → "female","male"
###   sample:    "orig.ident"
###   batch:     "batch"
###   cell type: "celltype"
###   % mito:    "percent.mt"
###
### OUTPUT: <PATH_TO_OUTPUT>/Figures/ED2_QC/
################################################################################

suppressPackageStartupMessages({
  library(conflicted)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  conflicts_prefer(dplyr::filter,  .quiet = TRUE)
  conflicts_prefer(dplyr::select,  .quiet = TRUE)
  conflicts_prefer(dplyr::rename,  .quiet = TRUE)
  conflicts_prefer(dplyr::mutate,  .quiet = TRUE)
})

################################################################################
### PATHS
################################################################################

fig_date <- format(Sys.Date(), "%Y-%m-%d")

BASE <- "PATH_TO_OUTPUT"     # same folder as BASE in Fig2_ED4_Pipeline.R
if (grepl("^PATH_TO_", BASE)) stop("Set BASE at the top of this script")

OUT_DIR <- paste0(file.path(BASE, "Figures", "ED2_QC"), "/")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("Output:", OUT_DIR, "\n\n")

DESeq2_DIR <- paste0(file.path(BASE, "DESeq2_Results"), "/")

################################################################################
### PARAMETERS
################################################################################

AGE_KEY    <- "age"
COND_KEY   <- "experimental.condition"
SEX_KEY    <- "sex"
CT_KEY     <- "celltype"
SAMPLE_KEY <- "orig.ident"
BATCH_KEY  <- "batch"
MITO_KEY   <- "percent.mt"

# Age values → display labels (O not A)
AGE_LABS <- c("young" = "Y", "Middleage" = "M", "old" = "O")
SEX_LABS <- c("female" = "F", "male" = "M")

# All pairwise DESeq2 comparisons
DE_COMPS <- c(
  "YEEYHC", "MEEMHC", "OEEOHC",
  "OHCYHC", "MHCYHC", "OHCMHC",
  "OEEYHC", "MEEYHC",
  "OEEYEE", "MEEYEE", "OEEMEE"
)
DE_COMPS <- unique(DE_COMPS)

CT_ORDER <- c(
  "DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN",
  "Astrocytes","Oligodendrocytes","OPC","Microglia",
  "Endothelial","Pericyte","APC","Stromal",
  "Cajal-Retzius Cell","Macrophage","T Cell"
)

CT_COLS <- c(
  "DG ExNs"           = "#E53935",
  "CA1 ExNs"          = "#CC0000",
  "CA2 ExNs"          = "#8E24AA",
  "CA3 ExNs"          = "#5E35B1",
  "Sub ExNs"          = "#80D8FF",
  "IN"                = "#00897B",
  "Astrocytes"        = "#43A047",
  "Oligodendrocytes"  = "#1E88E5",
  "OPC"               = "#00ACC1",
  "Microglia"         = "#7CB342",
  "Endothelial"       = "#FB8C00",
  "Pericyte"          = "#6D4C41",
  "APC"               = "#546E7A",
  "Stromal"           = "#27408B",
  "Macrophage"        = "#C62828",
  "Cajal-Retzius Cell"= "#FDD835",
  "T Cell"            = "#757575"
)

# 12 condition colors: age hue (blue=Young, green=Middle, red=Old)
#   HC = lighter/pastel shade  |  EE/NE = saturated/darker shade
#   Female = slightly warmer within hue; Male = slightly cooler
COND_COLS <- c(
  # ── Young (blue) ──────────────────────────────────────────
  "YF_HC" = "#90CAF9",   # light blue,       female HC
  "YF_EE" = "#1565C0",   # deep blue,         female NE
  "YM_HC" = "#B3E5FC",   # pale sky blue,     male HC
  "YM_EE" = "#0277BD",   # medium-dark blue,  male NE
  # ── Middle (green) ────────────────────────────────────────
  "MF_HC" = "#A5D6A7",   # light green,       female HC
  "MF_EE" = "#2E7D32",   # forest green,      female NE
  "MM_HC" = "#C8E6C9",   # pale mint,         male HC
  "MM_EE" = "#388E3C",   # medium green,      male NE
  # ── Old (red) ─────────────────────────────────────────────
  "OF_HC" = "#EF9A9A",   # light red,         female HC
  "OF_EE" = "#C62828",   # deep red,          female NE
  "OM_HC" = "#FFCDD2",   # pale rose,         male HC
  "OM_EE" = "#B71C1C"    # dark crimson,      male NE
)

save_both <- function(p, stem, w, h, dpi = 300) {
  w <- min(w, 50); h <- min(h, 50)
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h,
         device = "pdf", limitsize = FALSE)
  ggsave(paste0(stem, ".png"), p, width = w, height = h,
         device = "png", dpi = dpi, bg = "white", limitsize = FALSE)
  cat("  Saved:", basename(stem), "\n")
}

################################################################################
### CHECK SEURAT OBJECT + BUILD METADATA
################################################################################

if (!exists("df_harmony")) {
  stop("df_harmony not found. Load the EE Seurat object first.")
}
miss_qc <- setdiff(c("nFeature_RNA", "nCount_RNA", MITO_KEY), colnames(df_harmony@meta.data))
if (length(miss_qc)) stop("Panel F needs metadata column(s) missing from this object: ", paste(miss_qc, collapse=", "))

meta <- df_harmony@meta.data %>%
  dplyr::select(-dplyr::any_of("cell_id")) %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::mutate(
    age_lab        = dplyr::recode(.data[[AGE_KEY]],  !!!AGE_LABS),
    sex_lab        = dplyr::recode(.data[[SEX_KEY]],  !!!SEX_LABS),
    cond_lab       = .data[[COND_KEY]],
    condition_full = paste0(age_lab, sex_lab, "_", cond_lab),
    celltype       = factor(.data[[CT_KEY]], levels = CT_ORDER)
  )

# Order samples: female → male, Young → Middle → Old, HC → EE (NE)
sample_order_df <- meta %>%
  dplyr::select(sample = .data[[SAMPLE_KEY]],
                sex    = .data[[SEX_KEY]],
                age    = .data[[AGE_KEY]],
                cond   = .data[[COND_KEY]]) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    sex_ord  = match(tolower(sex),  c("female", "male")),
    age_ord  = match(tolower(age),  c("young", "middleage", "old")),
    cond_ord = match(cond,          c("HC", "EE"))
  ) %>%
  dplyr::arrange(sex_ord, age_ord, cond_ord)

all_samples <- sample_order_df$sample
n_samples   <- length(all_samples)
cat(sprintf("Samples  : %d\n", n_samples))
cat(sprintf("Cells    : %d\n", nrow(meta)))
cat(sprintf("CTs      : %s\n",
            paste(intersect(CT_ORDER, unique(as.character(meta$celltype))),
                  collapse = ", ")))
cat(sprintf("Conditions: %s\n",
            paste(sort(unique(meta$condition_full)), collapse = ", ")))

# UMAP coords
umap_mat <- tryCatch(Embeddings(df_harmony, "umap"),
                     error = function(e) Embeddings(df_harmony, "UMAP"))
umap_df  <- as.data.frame(umap_mat) %>%
  setNames(c("UMAP1","UMAP2")) %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::left_join(meta, by = "cell_id")

pt_size <- max(0.01, 50000 / nrow(umap_df))

# Condition display order
cond_order <- c(
  "YF_HC","YF_EE","MF_HC","MF_EE","OF_HC","OF_EE",
  "YM_HC","YM_EE","MM_HC","MM_EE","OM_HC","OM_EE"
)
cond_order <- intersect(cond_order, unique(umap_df$condition_full))

################################################################################
### PANEL A — UMAP SPLITS
################################################################################

cat("\n=== Panel A: UMAP splits ===\n")

umap_theme <- theme_void(base_size = 7) +
  theme(
    plot.title      = element_text(face = "bold", size = 7.5, hjust = 0.5),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA)
  )

# A1: by condition
umap_cond_plots <- lapply(cond_order, function(cond) {
  col <- if (cond %in% names(COND_COLS)) COND_COLS[[cond]] else "grey50"
  ggplot() +
    geom_point(data = umap_df %>% dplyr::filter(condition_full != cond),
               aes(UMAP1, UMAP2), color = "grey92",
               size = pt_size * 0.4, stroke = 0) +
    geom_point(data = umap_df %>% dplyr::filter(condition_full == cond),
               aes(UMAP1, UMAP2), color = col,
               size = pt_size, stroke = 0, alpha = 0.8) +
    umap_theme + labs(title = cond)
})
p_a1 <- patchwork::wrap_plots(umap_cond_plots, ncol = 6) +
  patchwork::plot_annotation(
    title = "Condition",
    theme = theme(plot.title = element_text(face="bold", size=9, hjust=0.5))
  )

# A2: by sample
sample_cols <- setNames(scales::hue_pal()(n_samples), all_samples)
umap_sample_plots <- lapply(all_samples, function(s) {
  ggplot() +
    geom_point(data = umap_df %>% dplyr::filter(.data[[SAMPLE_KEY]] != s),
               aes(UMAP1, UMAP2), color = "grey92",
               size = pt_size * 0.4, stroke = 0) +
    geom_point(data = umap_df %>% dplyr::filter(.data[[SAMPLE_KEY]] == s),
               aes(UMAP1, UMAP2), color = sample_cols[s],
               size = pt_size, stroke = 0, alpha = 0.8) +
    umap_theme + labs(title = s)
})
n_scols <- min(8, ceiling(sqrt(n_samples)))
p_a2 <- patchwork::wrap_plots(umap_sample_plots, ncol = n_scols) +
  patchwork::plot_annotation(
    title = "Sample",
    theme = theme(plot.title = element_text(face="bold", size=9, hjust=0.5))
  )

# A3: by batch (optional)
if (BATCH_KEY %in% colnames(umap_df)) {
  batches <- sort(unique(umap_df[[BATCH_KEY]]))
  n_b     <- max(3, length(batches))
  batch_pal <- setNames(
    RColorBrewer::brewer.pal(n_b, "Dark2")[seq_along(batches)],
    as.character(batches)
  )
  umap_batch_plots <- lapply(batches, function(b) {
    ggplot() +
      geom_point(data = umap_df %>%
                   dplyr::filter(.data[[BATCH_KEY]] != b),
                 aes(UMAP1, UMAP2), color = "grey92",
                 size = pt_size * 0.4, stroke = 0) +
      geom_point(data = umap_df %>%
                   dplyr::filter(.data[[BATCH_KEY]] == b),
                 aes(UMAP1, UMAP2),
                 color = batch_pal[as.character(b)],
                 size = pt_size, stroke = 0, alpha = 0.8) +
      umap_theme + labs(title = paste0("Batch ", b))
  })
  p_a3 <- patchwork::wrap_plots(umap_batch_plots,
                                  ncol = length(batches)) +
    patchwork::plot_annotation(
      title = "Batch",
      theme = theme(plot.title = element_text(face="bold", size=9,
                                               hjust=0.5))
    )
  p_a <- p_a1 / p_a2 / p_a3
} else {
  cat("  Batch key not found — skipping batch UMAP\n")
  p_a <- p_a1 / p_a2
}

save_both(p_a, file.path(OUT_DIR, paste0(fig_date, "_A_UMAP_Splits")),
          w = max(14, n_scols * 2.4), h = 20)

################################################################################
### PANEL B — CONDITION PROPORTION BY CELL TYPE (averaged, 12 conditions)
################################################################################

cat("\n=== Panel B: Condition proportion by cell type ===\n")

prop_b <- meta %>%
  dplyr::count(condition_full, celltype, .drop = FALSE) %>%
  dplyr::group_by(condition_full) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    condition_full = factor(condition_full, levels = cond_order),
    celltype       = factor(celltype, levels = rev(CT_ORDER))
  ) %>%
  dplyr::filter(!is.na(condition_full))

p_b <- ggplot(prop_b, aes(x = condition_full, y = prop, fill = celltype)) +
  geom_col(width = 0.85, color = NA) +
  scale_fill_manual(
    values = CT_COLS[intersect(names(CT_COLS),
                                levels(prop_b$celltype))],
    name   = "Cell type",
    guide  = guide_legend(reverse = TRUE, ncol = 1,
                           override.aes = list(size = 3))
  ) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle    = element_text(size = 8, color = "grey40", hjust = 0.5),
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.4, "cm"),
    legend.position  = "right"
  ) +
  labs(
    title    = "Condition proportion by major cell type",
    subtitle = "Averaged across samples within each sex × age × condition",
    x = NULL, y = "Proportion"
  )

save_both(p_b,
          file.path(OUT_DIR, paste0(fig_date, "_B_ConditionProportion")),
          w = max(10, length(cond_order) * 0.8 + 4), h = 6)

################################################################################
### PANEL C — SAMPLE PROPORTION BY CELL TYPE
################################################################################

cat("\n=== Panel C: Sample proportion by cell type ===\n")

prop_c <- meta %>%
  dplyr::count(.data[[SAMPLE_KEY]], celltype, .drop = FALSE) %>%
  dplyr::rename(sample = .data[[SAMPLE_KEY]]) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    sample   = factor(sample, levels = all_samples),
    celltype = factor(celltype, levels = rev(CT_ORDER))
  )

p_c <- ggplot(prop_c, aes(x = sample, y = prop, fill = celltype)) +
  geom_col(width = 0.85, color = NA) +
  scale_fill_manual(
    values = CT_COLS[intersect(names(CT_COLS), levels(prop_c$celltype))],
    name   = "Cell type",
    guide  = guide_legend(reverse = TRUE, ncol = 1,
                           override.aes = list(size = 3))
  ) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7.5),
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "right"
  ) +
  labs(title = "Sample proportion by major cell type",
       x = NULL, y = "Proportion")

save_both(p_c,
          file.path(OUT_DIR, paste0(fig_date, "_C_SampleProportion")),
          w = max(12, n_samples * 0.55 + 4), h = 6)

################################################################################
### PANEL D — CELL TYPE COUNT BY SAMPLE
################################################################################

cat("\n=== Panel D: Cell type count by sample ===\n")

count_d <- meta %>%
  dplyr::count(.data[[SAMPLE_KEY]], celltype) %>%
  dplyr::rename(sample = .data[[SAMPLE_KEY]]) %>%
  dplyr::mutate(
    sample   = factor(sample, levels = all_samples),
    celltype = factor(celltype, levels = rev(CT_ORDER))
  )

p_d <- ggplot(count_d, aes(x = sample, y = n, fill = celltype)) +
  geom_col(width = 0.85, color = NA) +
  scale_fill_manual(
    values = CT_COLS[intersect(names(CT_COLS), levels(count_d$celltype))],
    name   = "Cell type",
    guide  = guide_legend(reverse = TRUE, ncol = 1,
                           override.aes = list(size = 3))
  ) +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0)) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7.5),
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "right"
  ) +
  labs(title = "Cell type count by sample", x = NULL, y = "Count")

save_both(p_d,
          file.path(OUT_DIR, paste0(fig_date, "_D_CellTypeCount")),
          w = max(12, n_samples * 0.55 + 4), h = 6)

################################################################################
### PANEL E — DE GENE COUNTS PER CT × COMPARISON
################################################################################

cat("\n=== Panel E: DE gene counts ===\n")

de_rows <- list()
for (sx in c("female", "male")) {
  for (comp in DE_COMPS) {
    # Try two filename patterns
    candidates <- c(
      list.files(DESeq2_DIR,
                 pattern = paste0(sx, "_ALL_", comp, "\\.csv$"),
                 full.names = TRUE),
      list.files(DESeq2_DIR,
                 pattern = paste0("_", comp, "_ALL\\.csv$"),
                 full.names = TRUE)
    )
    if (length(candidates) == 0) next
    f  <- candidates[length(candidates)]
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (is.null(df) || !"padj" %in% colnames(df)) next

    cts_in <- intersect(CT_ORDER, unique(df$celltype))
    for (ct in cts_in) {
      n_de <- sum(df$padj[df$celltype == ct] < 0.05, na.rm = TRUE)
      de_rows[[paste(sx, comp, ct, sep = "__")]] <- data.frame(
        sex = sx, comp = comp, celltype = ct, n_de = n_de,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(de_rows) > 0) {
  de_df <- dplyr::bind_rows(de_rows) %>%
    dplyr::mutate(
      celltype = factor(celltype, levels = CT_ORDER),
      comp     = factor(comp, levels = rev(DE_COMPS))
    )

  # Group CTs for subplot panels
  ct_groups <- list(
    "Excitatory Neurons" = c("DG ExNs","CA3 ExNs","CA1 ExNs",
                              "CA2 ExNs","Sub ExNs"),
    "Inhibitory Neurons" = "IN",
    "Glia"               = c("Astrocytes","Oligodendrocytes",
                              "OPC","Microglia"),
    "Vascular / Other"   = c("Endothelial","Pericyte","APC","Stromal",
                              "Macrophage","Cajal-Retzius Cell","T Cell")
  )

  e_plots <- lapply(names(ct_groups), function(grp) {
    cts_g <- intersect(ct_groups[[grp]], CT_ORDER)
    df_g  <- de_df %>% dplyr::filter(celltype %in% cts_g)
    if (nrow(df_g) == 0 || all(df_g$n_de == 0)) return(NULL)
    strip_col <- CT_COLS[cts_g[1]]
    if (is.na(strip_col)) strip_col <- "grey35"

    ggplot(df_g, aes(x = n_de, y = comp, fill = sex)) +
      geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
      scale_fill_manual(
        values = c(female = "#C62828", male = "#0D47A1"), name = NULL
      ) +
      facet_wrap(~celltype, ncol = 2, scales = "free_x") +
      theme_classic(base_size = 8) +
      theme(
        strip.text       = element_text(face="bold", size=8, color="white"),
        strip.background = element_rect(fill = strip_col, color = NA),
        axis.text.y      = element_text(size = 7),
        plot.title       = element_text(face="bold", size=10, hjust=0.5),
        legend.position  = "top"
      ) +
      labs(title = grp, x = "N DE genes (padj < 0.05)", y = NULL)
  })

  e_plots <- e_plots[!sapply(e_plots, is.null)]
  if (length(e_plots) > 0) {
    p_e <- patchwork::wrap_plots(e_plots, ncol = 2) +
      patchwork::plot_annotation(
        title = "DE gene counts across pairwise comparisons",
        theme = theme(plot.title = element_text(face="bold", size=12,
                                                 hjust=0.5))
      )
    save_both(p_e,
              file.path(OUT_DIR, paste0(fig_date, "_E_DEGeneCounts")),
              w = 16, h = max(10, length(DE_COMPS) * 0.55 + 4))
  }
  write.csv(de_df,
            file.path(OUT_DIR, "DE_Gene_Counts_Summary.csv"),
            row.names = FALSE)
} else {
  cat("  No DE results found — skipping panel E\n")
}

################################################################################
### PANEL F — RNA QC VIOLINS
################################################################################

cat("\n=== Panel F: RNA QC violins ===\n")

qc_df <- meta %>%
  dplyr::rename(sample = .data[[SAMPLE_KEY]]) %>%
  dplyr::select(sample, condition_full,
                nFeature = nFeature_RNA,
                nCount   = nCount_RNA,
                pct_mito = .data[[MITO_KEY]]) %>%
  dplyr::mutate(sample = factor(sample, levels = all_samples))

make_qc_violin <- function(df, y_var, y_lab, log_y = FALSE) {
  p <- ggplot(df, aes(x = sample, y = .data[[y_var]],
                       fill = condition_full)) +
    geom_violin(scale = "width", linewidth = 0.3, alpha = 0.85) +
    geom_boxplot(width = 0.14, outlier.shape = NA,
                 fill = "white", linewidth = 0.4, alpha = 0.85) +
    scale_fill_manual(
      values = COND_COLS[intersect(names(COND_COLS),
                                    unique(df$condition_full))],
      name = NULL, guide = "none"
    ) +
    theme_classic(base_size = 9) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
      plot.title         = element_text(face="bold", size=10, hjust=0.5),
      panel.grid.major.y = element_line(color="grey93", linewidth=0.3)
    ) +
    labs(title = y_lab, x = NULL, y = y_lab)
  if (log_y) p <- p + scale_y_log10(labels = scales::comma)
  p
}

p_f1 <- make_qc_violin(qc_df, "nFeature", "N genes (nFeature_RNA)",
                        log_y = TRUE)
p_f2 <- make_qc_violin(qc_df, "nCount",   "N UMIs (nCount_RNA)",
                        log_y = TRUE)
p_f3 <- make_qc_violin(qc_df, "pct_mito", "% Mitochondrial reads")

p_f <- p_f1 / p_f2 / p_f3 +
  patchwork::plot_annotation(
    title = "RNA QC per sample",
    theme = theme(plot.title = element_text(face="bold", size=12, hjust=0.5))
  )

save_both(p_f,
          file.path(OUT_DIR, paste0(fig_date, "_F_RNA_QC_Violins")),
          w = max(14, n_samples * 0.55 + 4), h = 14)

################################################################################
### COMBINED PANEL
################################################################################

cat("\n=== Combined panel ===\n")

n_cts <- length(intersect(CT_ORDER, unique(as.character(meta[[CT_KEY]]))))

p_combined <- (p_b | p_d) /
  p_f +
  patchwork::plot_layout(heights = c(1, 1.5)) +
  patchwork::plot_annotation(
    title    = "NE Dataset — Extended QC Figure",
    subtitle = sprintf(
      "%d samples | %d conditions | %d cell types",
      n_samples, length(cond_order), n_cts
    ),
    tag_levels = "a",
    theme = theme(
      plot.title      = element_text(face="bold", size=14, hjust=0.5),
      plot.subtitle   = element_text(size=9, color="grey40", hjust=0.5),
      plot.background = element_rect(fill="white", color=NA)
    )
  )

save_both(p_combined,
          file.path(OUT_DIR, paste0(fig_date, "_Combined_QC_Panel")),
          w = 22, h = 26)

cat("\n=== DONE ===\nOutput:", OUT_DIR, "\n")
cat("Files:\n")
cat("  A_UMAP_Splits          — condition / sample / batch UMAPs\n")
cat("  B_ConditionProportion  — 12-condition averaged stacked bar\n")
cat("  C_SampleProportion     — per-sample stacked bar\n")
cat("  D_CellTypeCount        — per-sample count stacked bar\n")
cat("  E_DEGeneCounts         — pairwise DE counts per CT per sex\n")
cat("  F_RNA_QC_Violins       — nFeature / nCount / %mito\n")
cat("  Combined_QC_Panel      — b+d+f combined\n")
cat("  DE_Gene_Counts_Summary.csv\n")
