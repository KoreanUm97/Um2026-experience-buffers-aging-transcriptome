################################################################################
### ED8de_mScarlet_QC.R
###
### Extended Data Fig. 8d,e: is viral transduction (mScarlet) comparable across the four groups (young / old x PBS / C21)
### at each processing step? Neurons only (the cells the mScarlet steps apply to); female DREADD dataset.
###
###   Step                    Cells (neurons)
###   All mScarlet            every neuron
###   mScarlet >= 50          neurons with mScarlet raw UMI >= 50
###   Covariate               same cells as ">= 50" (Approach B adjusts for mScarlet in the model, not by removing cells)
###   Trim and downsample     exactly the Approach D neurons of DREADD_DESeq2.R: per cell type, 10-90th percentile
###                           window (median across mice), then every cell type downsampled to its smallest mouse x
###                           group count (same seed, same order). Or read from Approach D's _CellsUsed.csv (CELLS_D_CSV).
###
###   ED8d  per-mouse metrics (CV, IQR, median mScarlet UMI, number of neurons) at each step; one line per mouse
###   ED8e  mScarlet (log1p UMI) per group at each step: cells (max 200 per mouse), box, per-mouse medians (diamonds)
###   Test (both panels): Kruskal-Wallis across the four groups on per-mouse values (one value per mouse), per step
###   (and per metric in d). ns = no significant difference between groups.
###
### Run:  Rscript ED8de_mScarlet_QC.R     (set OBJ_RDS and OUT_DIR below first)
################################################################################

suppressPackageStartupMessages({
  library(Seurat); library(BPCells)
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(conflicted)
  conflicts_prefer(dplyr::select, .quiet = TRUE); conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicts_prefer(dplyr::mutate, .quiet = TRUE); conflicts_prefer(base::intersect, .quiet = TRUE)
  conflicts_prefer(base::setdiff, .quiet = TRUE)
})

################################################################################
### SETTINGS
################################################################################

OBJ_RDS     <- "PATH_TO_DREADD_OBJECT.rds"
OUT_DIR     <- file.path("PATH_TO_OUTPUT", "Paper_Figures", "Figure4", "mScarlet_QC")
CELLS_D_CSV <- NULL       # optional: path to <date>_CellsUsed.csv from DREADD_DESeq2.R (Approach D); NULL = recompute
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MSCARLET_THRESH <- 50
TRIM_LO <- 0.10; TRIM_HI <- 0.90
MIN_CELLS_FLOOR <- 5
SEED <- 42                # same as DREADD_DESeq2.R
JITTER_MAX <- 200         # cells drawn per mouse x step in ED8e (display only)

NEURONAL_CTS <- c("CA1 ExNs","CA2 ExNs","CA3 ExNs","DG ExNs","Sub ExNs","IN")
GROUP_COLS <- c(young_PBS = "#A6CEE3", young_C21 = "#1F78B4", old_PBS = "#FDBF6F", old_C21 = "#E31A1C")
GROUP_LABS <- c(young_PBS = "Young PBS", young_C21 = "Young C21", old_PBS = "Old PBS", old_C21 = "Old C21")
STEPS <- c("All mScarlet", "mScarlet >= 50", "Covariate", "Trim and downsample")

save_both <- function(p, stem, w, h) {
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h, device = "pdf", limitsize = FALSE)
  ggsave(paste0(stem, ".png"), p, width = w, height = h, device = "png", dpi = 600, bg = "white", limitsize = FALSE)
  cat("  Saved:", basename(stem), "\n")
}
sig_label <- function(p) dplyr::case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")

################################################################################
### LOAD + METADATA
################################################################################

cat("Loading", OBJ_RDS, "\n")
obj <- readRDS(OBJ_RDS)
in_present <- intersect(c("IN1", "IN2"), unique(obj@meta.data$celltype))
if (length(in_present)) obj@meta.data$celltype[obj@meta.data$celltype %in% in_present] <- "IN"
if (!"mScarlet_raw" %in% colnames(obj@meta.data))
  obj@meta.data$mScarlet_raw <- FetchData(obj, vars = "mScarlet", layer = "counts")$mScarlet

meta <- obj@meta.data[, !colnames(obj@meta.data) %in% "cell_id", drop = FALSE] %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::mutate(age_norm = dplyr::case_when(tolower(age) == "young" ~ "young", tolower(age) == "old" ~ "old",
                                            TRUE ~ tolower(age)),
                age_condition = paste0(age_norm, "_", experimental.condition)) %>%
  dplyr::select(cell_id, celltype, orig.ident, age_condition, mScarlet_raw)
rm(obj); invisible(gc())

################################################################################
### CELLS AT EACH STEP (neurons)
################################################################################

# Approach D cells, computed exactly as in DREADD_DESeq2.R (glia included so the random draw is the same), then neurons kept
approach_d_cells <- function(meta) {
  base <- meta %>% dplyr::filter((celltype %in% NEURONAL_CTS & mScarlet_raw >= MSCARLET_THRESH) | (!celltype %in% NEURONAL_CTS))
  trimmed <- list()
  for (ct in unique(base$celltype)) {
    cm <- base %>% dplyr::filter(celltype == ct)
    if (!ct %in% NEURONAL_CTS) { trimmed[[ct]] <- cm; next }
    b <- dplyr::bind_rows(lapply(unique(cm$orig.ident), function(m) {
      v <- cm$mScarlet_raw[cm$orig.ident == m]
      if (length(v) < 5) NULL else data.frame(lo = quantile(v, TRIM_LO), hi = quantile(v, TRIM_HI))
    }))
    if (!nrow(b) || median(b$lo) >= median(b$hi)) { trimmed[[ct]] <- cm; next }
    trimmed[[ct]] <- cm %>% dplyr::filter(mScarlet_raw >= median(b$lo) & mScarlet_raw <= median(b$hi))
  }
  mt <- dplyr::bind_rows(trimmed)
  tgt <- mt %>% dplyr::count(celltype, orig.ident, age_condition) %>% dplyr::group_by(celltype) %>%
    dplyr::summarise(target_n = min(n), .groups = "drop") %>% dplyr::filter(target_n >= MIN_CELLS_FLOOR)
  set.seed(SEED)
  dplyr::bind_rows(lapply(tgt$celltype, function(ct) {
    t <- tgt$target_n[tgt$celltype == ct]
    mt %>% dplyr::filter(celltype == ct) %>% dplyr::group_by(orig.ident, age_condition) %>%
      dplyr::group_modify(~dplyr::slice_sample(.x, n = min(t, nrow(.x)))) %>% dplyr::ungroup()
  }))$cell_id
}

neur    <- meta %>% dplyr::filter(celltype %in% NEURONAL_CTS)
d_cells <- if (!is.null(CELLS_D_CSV)) read.csv(CELLS_D_CSV, stringsAsFactors = FALSE)$cell_id else approach_d_cells(meta)
steps <- dplyr::bind_rows(
  neur %>% dplyr::mutate(step = STEPS[1]),
  neur %>% dplyr::filter(mScarlet_raw >= MSCARLET_THRESH) %>% dplyr::mutate(step = STEPS[2]),
  neur %>% dplyr::filter(mScarlet_raw >= MSCARLET_THRESH) %>% dplyr::mutate(step = STEPS[3]),
  neur %>% dplyr::filter(cell_id %in% d_cells) %>% dplyr::mutate(step = STEPS[4])) %>%
  dplyr::mutate(step = factor(step, levels = STEPS), grp_f = factor(age_condition, levels = names(GROUP_COLS)),
                mScarlet_log1p = log1p(mScarlet_raw))
cat("Neurons per step:\n"); print(table(steps$step))

################################################################################
### PER-MOUSE VALUES + KRUSKAL-WALLIS ACROSS GROUPS
################################################################################

mouse <- steps %>% dplyr::group_by(step, orig.ident, grp_f) %>%
  dplyr::summarise(median_msc = median(mScarlet_raw), iqr_msc = IQR(mScarlet_raw),
                   cv_msc = sd(mScarlet_raw) / max(mean(mScarlet_raw), 1), n_cells = dplyr::n(),
                   median_log1p = median(mScarlet_log1p), .groups = "drop")
METRICS <- c(cv_msc = "CV (SD/Mean)", iqr_msc = "IQR", median_msc = "Median mScarlet (UMI)", n_cells = "Number of mScarlet+ Neurons")
kw <- dplyr::bind_rows(lapply(STEPS, function(s) dplyr::bind_rows(lapply(names(METRICS), function(m) {
  d <- mouse[mouse$step == s, ]
  p <- tryCatch(kruskal.test(d[[m]] ~ d$grp_f)$p.value, error = function(e) NA_real_)
  data.frame(step = s, metric = m, n_mice = nrow(d), kw_p = p, sig = sig_label(p))
})))) %>% dplyr::mutate(step = factor(step, levels = STEPS))
utils::write.csv(mouse, file.path(OUT_DIR, "mScarlet_QC_PerMouse.csv"), row.names = FALSE)
utils::write.csv(kw, file.path(OUT_DIR, "mScarlet_QC_KruskalWallis_AcrossGroups.csv"), row.names = FALSE)
cat("\nKruskal-Wallis across the four groups (per-mouse values):\n"); print(kw, row.names = FALSE)

################################################################################
### ED8d: COMPARABILITY METRICS
################################################################################

ml <- mouse %>% tidyr::pivot_longer(names(METRICS), names_to = "metric", values_to = "value") %>%
  dplyr::mutate(metric_lab = factor(METRICS[metric], levels = METRICS))
ann <- kw %>% dplyr::left_join(ml %>% dplyr::group_by(metric) %>% dplyr::summarise(y = max(value) * 1.08, .groups = "drop"),
                               by = "metric") %>% dplyr::mutate(metric_lab = factor(METRICS[metric], levels = METRICS))
p_d <- ggplot(ml, aes(step, value, color = grp_f, group = interaction(orig.ident, grp_f))) +
  geom_line(linewidth = 0.65, alpha = 0.7) + geom_point(size = 2.5, alpha = 0.9) +
  geom_text(data = ann, aes(step, y, label = sig), inherit.aes = FALSE, size = 3.5, color = "grey20", vjust = 0) +
  scale_color_manual(values = GROUP_COLS, labels = GROUP_LABS, name = NULL) + scale_x_discrete(drop = FALSE) +
  facet_wrap(~metric_lab, scales = "free_y", ncol = 2) + theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8), strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "grey92"), plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8.5, color = "grey40", hjust = 0.5), legend.position = "top") +
  labs(title = "mScarlet Comparability Metrics",
       subtitle = "Neurons | each line = one mouse | label = Kruskal-Wallis across the four groups (per-mouse values)",
       x = NULL, y = "Value")
save_both(p_d, file.path(OUT_DIR, "ED8d_Comparability_Metrics"), 11, 9)

################################################################################
### ED8e: DISTRIBUTIONS + PER-MOUSE MEDIANS
################################################################################

set.seed(SEED)
jit <- steps %>% dplyr::group_by(orig.ident, step) %>%
  dplyr::group_modify(~ if (nrow(.x) <= JITTER_MAX) .x else dplyr::slice_sample(.x, n = JITTER_MAX)) %>% dplyr::ungroup()
kw_e <- kw %>% dplyr::filter(metric == "median_msc") %>% tidyr::crossing(grp_f = factor(names(GROUP_COLS), levels = names(GROUP_COLS))) %>%
  dplyr::mutate(y = max(steps$mScarlet_log1p) * 1.04)
p_e <- ggplot() +
  geom_jitter(data = jit, aes(step, mScarlet_log1p, color = grp_f), width = 0.22, size = 0.3, alpha = 0.15) +
  geom_boxplot(data = steps, aes(step, mScarlet_log1p), fill = NA, color = "grey20", width = 0.5, linewidth = 0.45, outlier.shape = NA) +
  geom_point(data = mouse, aes(step, median_log1p), shape = 18, size = 3.5, color = "grey10", alpha = 0.9) +
  geom_hline(yintercept = log1p(MSCARLET_THRESH), linewidth = 0.6, color = "grey20", linetype = "dashed") +
  geom_text(data = kw_e, aes(step, y, label = sig), size = 3.5, color = "grey20", vjust = 0) +
  scale_color_manual(values = GROUP_COLS, guide = "none") +
  facet_wrap(~grp_f, ncol = 2, labeller = labeller(grp_f = GROUP_LABS)) + theme_classic(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 11, color = "white"), strip.background = element_rect(fill = "grey35", color = NA),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8.5), plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8.5, color = "grey40", hjust = 0.5)) +
  labs(title = "mScarlet Distribution - Cells + Per-Mouse Medians",
       subtitle = paste0("Neurons | points = cells (max ", JITTER_MAX, "/mouse) | box = median + IQR | diamond = per-mouse median\n",
                         "Label = Kruskal-Wallis across the four groups on per-mouse medians | dashed = ", MSCARLET_THRESH, " UMI"),
       x = NULL, y = "mScarlet (log1p UMI)")
save_both(p_e, file.path(OUT_DIR, "ED8e_Distribution_PerMouseMedians"), 11, 9)

cat("\n=== DONE ===\nOutputs:", OUT_DIR, "\n")
