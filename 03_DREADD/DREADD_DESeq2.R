################################################################################
### DREADD_DESeq2.R
###
### Pseudobulk DESeq2 for the DREADD dataset (female only; IN1 + IN2 merged into "IN"), two approaches for matching
### viral transduction (mScarlet) across mice, plus the strip plots of Fig. 4d.
###
###   Approach B (covariate)
###     neurons with mScarlet raw UMI >= 50; all glial / vascular cells
###     neurons: ~ mScarlet_median_scaled + grp  (per-mouse, per-cell-type median mScarlet UMI, z-scored)
###     others:  ~ grp
###   Approach D (trimming + downsampling)
###     neurons with mScarlet raw UMI >= 50; all glial / vascular cells
###     neurons: per cell type, 10th and 90th percentile of mScarlet in each mouse, median of those across mice
###              = common window; only cells inside the window are kept
###     every cell type (neurons and glia): downsampled to the smallest mouse x group cell count of that cell type;
###              cell types whose smallest count is < MIN_CELLS_FLOOR are not tested
###     all cell types: ~ grp
###   Pseudobulk = summed raw counts per cell type x mouse; genes with >= 10 counts in >= 2 samples; cell types tested
###   when >= 4 samples and >= 2 mice per group.
###
###   Comparisons   YC21YPBS  young C21 vs young PBS          OC21OPBS  old C21 vs old PBS
###                 OPBSYPBS  old PBS vs young PBS (aging)    OC21YC21  old C21 vs young C21
###                 OC21YPBS  old C21 vs young PBS
###
### Outputs (read by the other 03_DREADD scripts):
###   <OUT_ROOT>/DESeq2_mScarletPos_Approach{B,D}_MergedIN/<date>/<date>_<COMP>_ALL.csv   (+ one file per cell type)
###   <OUT_ROOT>/Paper_Figures/Figure4/StripPlots/Approach{B,D}/<COMP>_StripVolcano.pdf/.png
###     Fig. 4d = Approach D, YC21YPBS (top) and OC21OPBS (bottom)
###
### Run:  Rscript DREADD_DESeq2.R     (set OBJ_RDS and OUT_ROOT below first)
################################################################################

suppressPackageStartupMessages({
  library(conflicted)
  library(Seurat)
  library(BPCells)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  conflicts_prefer(dplyr::select,       .quiet = TRUE)
  conflicts_prefer(dplyr::filter,       .quiet = TRUE)
  conflicts_prefer(dplyr::mutate,       .quiet = TRUE)
  conflicts_prefer(DESeq2::plotMA,      .quiet = TRUE)
  conflicts_prefer(base::intersect,     .quiet = TRUE)
  conflicts_prefer(base::setdiff,       .quiet = TRUE)
  conflicts_prefer(base::union,         .quiet = TRUE)
  conflicts_prefer(base::as.data.frame, .quiet = TRUE)
})
USE_BEESWARM <- requireNamespace("ggbeeswarm", quietly = TRUE)

################################################################################
### SETTINGS
################################################################################

OBJ_RDS  <- "PATH_TO_DREADD_OBJECT.rds"      # DREADD Seurat object (GEO)
OUT_ROOT <- "PATH_TO_DREADD_OUTPUT"
RUN      <- c("B", "D")
RUN_DATE <- format(Sys.Date(), "%Y-%m-%d")    # output subfolder; downstream scripts read the latest one

MSCARLET_THRESH <- 50
PADJ_THR        <- 0.05
TRIM_LO         <- 0.10
TRIM_HI         <- 0.90
MIN_CELLS_FLOOR <- 5
SEED            <- 42

NEURONAL_CELLTYPES <- c("CA1 ExNs","CA2 ExNs","CA3 ExNs","DG ExNs","Sub ExNs","IN")
CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN",
              "Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal",
              "Cajal-Retzius Cell","Macrophage","Epithelial","T Cell")
COMPARISONS <- list(
  list(name = "YC21YPBS", label = "Young C21 vs Young PBS",              age = "young", num = "C21", den = "PBS"),
  list(name = "OC21OPBS", label = "Old C21 vs Old PBS",                  age = "old",   num = "C21", den = "PBS"),
  list(name = "OPBSYPBS", label = "Old PBS vs Young PBS (aging)",        condition_filter = "PBS", num = "old", den = "young"),
  list(name = "OC21YC21", label = "Old C21 vs Young C21",                condition_filter = "C21", num = "old", den = "young"),
  list(name = "OC21YPBS", label = "Old C21 vs Young PBS",                cross = TRUE))
APPROACH_DIRS <- c(B = "DESeq2_mScarletPos_ApproachB_MergedIN", D = "DESeq2_mScarletPos_ApproachD_MergedIN")
APPROACH_LABS <- c(B = "Covariate | mScarlet >= 50 | ~ mScarlet_median_scaled + grp (neurons)",
                   D = "Trimming + downsampling | mScarlet >= 50 | range trim 10-90th pct + downsample | ~ grp")
col_vals <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey82")

should_skip <- function(f) { if (file.exists(f)) { cat("  Skipping:", basename(f), "\n"); return(TRUE) }; FALSE }
save_both <- function(p, stem, w, h, dpi = 300) {
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h, device = "pdf", limitsize = FALSE)
  ggsave(paste0(stem, ".png"), p, width = w, height = h, device = "png", dpi = dpi, bg = "white", limitsize = FALSE)
  cat("  Saved:", basename(stem), "\n")
}

################################################################################
### LOAD + METADATA
################################################################################

cat("Loading", OBJ_RDS, "\n")
obj <- readRDS(OBJ_RDS)

in_present <- intersect(c("IN1", "IN2"), unique(obj@meta.data$celltype))
if (length(in_present) > 0) {
  obj@meta.data$celltype[obj@meta.data$celltype %in% in_present] <- "IN"
  cat(sprintf("Merged %s -> IN\n", paste(in_present, collapse = "+")))
}
if (!"mScarlet_raw" %in% colnames(obj@meta.data))
  obj@meta.data$mScarlet_raw <- FetchData(obj, vars = "mScarlet", layer = "counts")$mScarlet
obj@meta.data <- obj@meta.data %>%
  dplyr::mutate(age_norm       = dplyr::case_when(tolower(age) == "young" ~ "young",
                                                  tolower(age) == "old"   ~ "old",
                                                  TRUE ~ tolower(age)),
                condition_norm = experimental.condition,
                age_condition  = paste0(age_norm, "_", experimental.condition))

meta_full <- obj@meta.data[, !colnames(obj@meta.data) %in% "cell_id", drop = FALSE] %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, celltype, orig.ident, age_norm, condition_norm, age_condition, mScarlet_raw)
cat("Cells by age x condition:\n"); print(table(meta_full$age_condition))

# Neurons with mScarlet >= threshold; all other cell types
meta_base <- meta_full %>%
  dplyr::filter((celltype %in% NEURONAL_CELLTYPES & mScarlet_raw >= MSCARLET_THRESH) |
                (!celltype %in% NEURONAL_CELLTYPES))
cat(sprintf("After mScarlet >= %d: neuronal %d | other %d\n", MSCARLET_THRESH,
            sum(meta_base$celltype %in% NEURONAL_CELLTYPES), sum(!meta_base$celltype %in% NEURONAL_CELLTYPES)))

################################################################################
### SHARED: PSEUDOBULK + DESeq2
################################################################################

# AggregateExpression column names: {celltype}_{orig.ident with "_" turned into "-"}_{condition}_{age}.
# The cell type is matched from the left with the original names; orig.ident hyphens are turned back into underscores
# so it joins with the per-mouse mScarlet covariate.
parse_pb_colname <- function(s) {
  ct_name <- CT_ORDER[sapply(CT_ORDER, function(ct) startsWith(s, paste0(ct, "_")))]
  if (length(ct_name) == 0) return(list(celltype = NA, orig.ident = NA, condition = NA, age = NA))
  ct_name <- ct_name[1]
  parts <- strsplit(sub(paste0("^", ct_name, "_"), "", s), "_")[[1]]; n <- length(parts)
  list(celltype = ct_name, orig.ident = gsub("-", "_", paste(parts[1:(n - 2)], collapse = "_")),
       condition = parts[n - 1], age = parts[n])
}

pseudobulk <- function(cells) {
  so <- obj[, cells$cell_id]
  idx <- match(colnames(so), cells$cell_id)
  so@meta.data$age_norm       <- cells$age_norm[idx]
  so@meta.data$condition_norm <- cells$condition_norm[idx]
  cat(sprintf("Aggregating %d cells...\n", ncol(so)))
  pb <- AggregateExpression(object = so, assays = "RNA", return.seurat = FALSE,
                            group.by = c("celltype", "orig.ident", "condition_norm", "age_norm"),
                            normalization.method = NULL, slot = "counts", verbose = FALSE)$RNA
  parsed  <- lapply(colnames(pb), parse_pb_colname)
  pb_meta <- data.frame(sample = colnames(pb),
                        celltype = sapply(parsed, "[[", "celltype"), orig.ident = sapply(parsed, "[[", "orig.ident"),
                        condition = sapply(parsed, "[[", "condition"), age = sapply(parsed, "[[", "age"),
                        stringsAsFactors = FALSE)
  keep <- !is.na(pb_meta$celltype)
  if (any(!keep)) cat("Dropping", sum(!keep), "pseudobulk columns of cell types outside CT_ORDER\n")
  pb_meta$age_condition <- paste0(pb_meta$age, "_", pb_meta$condition)
  list(pb = pb[, keep, drop = FALSE], meta = pb_meta[keep, ])
}

# mcov: per-mouse, per-cell-type median mScarlet (Approach B neurons) or NULL
run_deseq2 <- function(PB, out_dir, ap, mcov = NULL) {
  all_res <- list()
  for (comp in COMPARISONS) {
    cat(sprintf("\n  --- %s: %s ---\n", comp$name, comp$label))
    m <- PB$meta
    if (!is.null(comp$age)) {
      sel <- m$age == comp$age; m <- m[sel, ]; m$grp <- factor(m$condition, levels = c(comp$den, comp$num))
    } else if (isTRUE(comp$cross)) {
      sel <- m$age_condition %in% c("old_C21", "young_PBS"); m <- m[sel, ]
      m$grp <- factor(m$age_condition, levels = c("young_PBS", "old_C21"))
    } else {
      sel <- m$condition == comp$condition_filter; m <- m[sel, ]; m$grp <- factor(m$age, levels = c(comp$den, comp$num))
    }
    pb_comp <- PB$pb[, sel, drop = FALSE]
    comp_res <- list()
    for (ct in intersect(CT_ORDER, unique(PB$meta$celltype))) {
      out_f <- file.path(out_dir, paste0(RUN_DATE, "_", comp$name, "_", ct, ".csv"))
      if (should_skip(out_f)) { comp_res[[ct]] <- read.csv(out_f, stringsAsFactors = FALSE); next }
      k <- m$celltype == ct
      cnt <- pb_comp[, k, drop = FALSE]; sm <- m[k, ]
      if (ncol(cnt) < 4) next
      if (length(table(sm$grp)) < 2 || any(table(sm$grp) < 2)) { cat(sprintf("    %-22s skipped (group sizes)\n", ct)); next }
      cnt <- cnt[rowSums(cnt >= 10) >= 2, , drop = FALSE]
      if (nrow(cnt) < 10) next
      is_neuronal <- ct %in% NEURONAL_CELLTYPES
      use_cov <- FALSE
      if (!is.null(mcov) && is_neuronal) {
        sm$mScarlet_median <- mcov$mScarlet_median[match(paste(sm$orig.ident, ct), paste(mcov$orig.ident, mcov$celltype))]
        use_cov <- !any(is.na(sm$mScarlet_median)) && var(sm$mScarlet_median) > 0
        if (!use_cov) cat(sprintf("    %-22s covariate skipped (NA or zero variance)\n", ct))
        if (use_cov) sm$mScarlet_median_scaled <- as.numeric(scale(sm$mScarlet_median))
      }
      res <- tryCatch({
        col_df <- if (use_cov) data.frame(grp = sm$grp, mScarlet_median_scaled = sm$mScarlet_median_scaled, row.names = colnames(cnt))
                  else data.frame(grp = sm$grp, row.names = colnames(cnt))
        dds <- DESeqDataSetFromMatrix(cnt, col_df, design = if (use_cov) ~mScarlet_median_scaled + grp else ~grp)
        dds <- DESeq(dds, quiet = TRUE)
        as.data.frame(results(dds, alpha = PADJ_THR)) %>% tibble::rownames_to_column("gene") %>%
          dplyr::mutate(celltype = ct, comparison = comp$name, covariate_used = use_cov, is_neuronal = is_neuronal,
                        approach = ap)
      }, error = function(e) { cat("    DESeq2 error for", ct, ":", e$message, "\n"); NULL })
      if (!is.null(res)) {
        cat(sprintf("    %-22s genes: %d | sig: %d | covariate: %s\n", ct, nrow(res),
                    sum(!is.na(res$padj) & res$padj < PADJ_THR), ifelse(use_cov, "yes", "no")))
        write.csv(res, out_f, row.names = FALSE); comp_res[[ct]] <- res
      }
    }
    if (length(comp_res)) {
      all_res[[comp$name]] <- dplyr::bind_rows(comp_res)
      write.csv(all_res[[comp$name]], file.path(out_dir, paste0(RUN_DATE, "_", comp$name, "_ALL.csv")), row.names = FALSE)
    }
  }
  all_res
}

################################################################################
### SHARED: STRIP PLOTS (Fig. 4d)
################################################################################

violin_offsets <- function(df, max_w = 0.35) {     # used only when ggbeeswarm is not installed
  df %>% dplyr::group_by(celltype) %>%
    dplyr::mutate(.bin = as.integer(cut(lfc_plot, breaks = 80, include.lowest = TRUE)), .bin = ifelse(is.na(.bin), 1L, .bin)) %>%
    dplyr::group_by(celltype, .bin) %>%
    dplyr::mutate(.bsize = dplyr::n(), .idx = seq_along(.bin),
                  .xoff = ((.idx - 1) - (.bsize - 1) / 2) / max(1, .bsize - 1) * max_w * (.bsize / 20)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(.xoff = pmax(pmin(.xoff, max_w), -max_w), .xnum = as.integer(celltype) + .xoff)
}

strip_plots <- function(all_res, ap) {
  vdir <- file.path(OUT_ROOT, "Paper_Figures", "Figure4", "StripPlots", paste0("Approach", ap))
  dir.create(vdir, recursive = TRUE, showWarnings = FALSE)
  for (comp in COMPARISONS) {
    if (is.null(all_res[[comp$name]])) next
    df <- all_res[[comp$name]] %>% dplyr::filter(!is.na(log2FoldChange)) %>%
      dplyr::mutate(sig = !is.na(padj) & padj < PADJ_THR,
                    color_cat = dplyr::case_when(sig & log2FoldChange > 0 ~ "Up", sig & log2FoldChange < 0 ~ "Down", TRUE ~ "NS"),
                    lfc_plot = log2FoldChange)
    ct_levels <- c(intersect(CT_ORDER, unique(df$celltype)), sort(setdiff(unique(df$celltype), CT_ORDER)))
    df$celltype <- factor(df$celltype, levels = ct_levels)
    if (USE_BEESWARM) {
      p <- ggplot(mapping = aes(x = celltype, y = lfc_plot)) +
        ggbeeswarm::geom_quasirandom(data = df %>% dplyr::filter(color_cat == "NS"), color = col_vals["NS"], size = 0.25, alpha = 0.25, width = 0.35) +
        ggbeeswarm::geom_quasirandom(data = df %>% dplyr::filter(color_cat != "NS"), aes(color = color_cat), size = 0.65, alpha = 0.80, width = 0.35)
    } else {
      df <- violin_offsets(df)
      p <- ggplot() +
        geom_point(data = df %>% dplyr::filter(color_cat == "NS"), aes(x = .xnum, y = lfc_plot), color = col_vals["NS"], size = 0.25, alpha = 0.25) +
        geom_point(data = df %>% dplyr::filter(color_cat != "NS"), aes(x = .xnum, y = lfc_plot, color = color_cat), size = 0.65, alpha = 0.80) +
        scale_x_continuous(breaks = seq_along(ct_levels), labels = ct_levels, expand = expansion(add = 0.6))
    }
    p <- p +
      geom_hline(yintercept = 0, linewidth = 0.5, color = "grey10") +
      scale_color_manual(values = col_vals, breaks = c("Up", "Down"), labels = c("up-regulated", "down-regulated"), name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
      guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1), keywidth = 0.4, keyheight = 0.4, nrow = 1)) +
      theme_classic(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5, color = "grey10"), axis.text.y = element_text(size = 8.5),
            axis.title.y = element_text(size = 10), axis.title.x = element_blank(),
            axis.line = element_line(color = "grey20", linewidth = 0.35), axis.ticks = element_line(color = "grey20", linewidth = 0.3),
            panel.grid = element_blank(), plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 4)),
            plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40", margin = margin(b = 6)),
            legend.position = "top", legend.justification = "right", legend.text = element_text(size = 8),
            legend.spacing.x = unit(0.15, "cm"), plot.margin = margin(6, 8, 4, 4)) +
      labs(title = paste0("DREADD  –  ", comp$label, "  –  Female"),
           subtitle = paste0("Approach ", ap, ": ", APPROACH_LABS[[ap]], " | padj < ", PADJ_THR), y = bquote(log[2](FC)))
    stem <- file.path(vdir, paste0(comp$name, "_StripVolcano"))
    if (!should_skip(paste0(stem, ".pdf"))) save_both(p, stem, max(7, length(ct_levels) * 0.52 + 2), 5)
  }
}

################################################################################
### APPROACH B: COVARIATE
################################################################################

if ("B" %in% RUN) {
  cat("\n========== Approach B (covariate) ==========\n")
  out_dir <- file.path(OUT_ROOT, APPROACH_DIRS[["B"]], RUN_DATE); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  mcov <- meta_base %>% dplyr::filter(celltype %in% NEURONAL_CELLTYPES) %>%
    dplyr::group_by(orig.ident, celltype) %>% dplyr::summarise(mScarlet_median = median(mScarlet_raw), .groups = "drop")
  write.csv(mcov, file.path(out_dir, paste0(RUN_DATE, "_mScarlet_Covariate.csv")), row.names = FALSE)
  write.csv(meta_base %>% dplyr::count(celltype, orig.ident, age_condition, name = "n_cells"),
            file.path(out_dir, paste0(RUN_DATE, "_CellCounts.csv")), row.names = FALSE)
  PB <- pseudobulk(meta_base)
  cat(sprintf("Pseudobulk samples: %d | orig.ident matched to covariate: %d\n", nrow(PB$meta),
              sum(PB$meta$orig.ident %in% mcov$orig.ident)))
  res_B <- run_deseq2(PB, out_dir, "B", mcov = mcov)
  strip_plots(res_B, "B")
}

################################################################################
### APPROACH D: TRIMMING + DOWNSAMPLING
################################################################################

if ("D" %in% RUN) {
  cat("\n========== Approach D (trimming + downsampling) ==========\n")
  out_dir <- file.path(OUT_ROOT, APPROACH_DIRS[["D"]], RUN_DATE); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Range trimming (neurons only)
  trimmed_list <- list()
  for (ct in unique(meta_base$celltype)) {
    ct_meta <- meta_base %>% dplyr::filter(celltype == ct)
    if (!ct %in% NEURONAL_CELLTYPES) { trimmed_list[[ct]] <- ct_meta; next }
    bounds <- lapply(unique(ct_meta$orig.ident), function(m) {
      vals <- ct_meta$mScarlet_raw[ct_meta$orig.ident == m]
      if (length(vals) < 5) return(NULL)
      data.frame(orig.ident = m, lo = quantile(vals, TRIM_LO), hi = quantile(vals, TRIM_HI), stringsAsFactors = FALSE)
    })
    bounds_df <- dplyr::bind_rows(bounds[!sapply(bounds, is.null)])
    if (nrow(bounds_df) == 0) { trimmed_list[[ct]] <- ct_meta; next }
    common_lo <- median(bounds_df$lo); common_hi <- median(bounds_df$hi)
    if (common_lo >= common_hi) { trimmed_list[[ct]] <- ct_meta; next }
    trimmed_list[[ct]] <- ct_meta %>% dplyr::filter(mScarlet_raw >= common_lo & mScarlet_raw <= common_hi)
    cat(sprintf("  %-22s: [%.1f, %.1f] | %d -> %d cells\n", ct, common_lo, common_hi, nrow(ct_meta), nrow(trimmed_list[[ct]])))
  }
  meta_trimmed <- dplyr::bind_rows(trimmed_list)
  write.csv(data.frame(celltype = names(trimmed_list),
                       n_before = sapply(names(trimmed_list), function(ct) sum(meta_base$celltype == ct)),
                       n_after  = sapply(trimmed_list, nrow)),
            file.path(out_dir, paste0(RUN_DATE, "_TrimReport.csv")), row.names = FALSE)

  # Downsample every cell type to its smallest mouse x group count
  target_per_ct <- meta_trimmed %>% dplyr::count(celltype, orig.ident, age_condition) %>%
    dplyr::group_by(celltype) %>% dplyr::summarise(target_n = min(n), .groups = "drop") %>%
    dplyr::filter(target_n >= MIN_CELLS_FLOOR)
  dropped <- setdiff(unique(meta_trimmed$celltype), target_per_ct$celltype)
  if (length(dropped)) cat("Not tested (fewer than", MIN_CELLS_FLOOR, "cells in a mouse):", paste(dropped, collapse = ", "), "\n")
  write.csv(target_per_ct, file.path(out_dir, paste0(RUN_DATE, "_DownsampleTargets.csv")), row.names = FALSE)
  set.seed(SEED)
  cells_D <- lapply(target_per_ct$celltype, function(ct) {
    tgt <- target_per_ct$target_n[target_per_ct$celltype == ct]
    meta_trimmed %>% dplyr::filter(celltype == ct) %>% dplyr::group_by(orig.ident, age_condition) %>%
      dplyr::group_modify(~dplyr::slice_sample(.x, n = min(tgt, nrow(.x)))) %>% dplyr::ungroup()
  }) %>% dplyr::bind_rows()
  write.csv(cells_D %>% dplyr::select(cell_id, celltype, orig.ident, age_condition, mScarlet_raw),
            file.path(out_dir, paste0(RUN_DATE, "_CellsUsed.csv")), row.names = FALSE)
  cat(sprintf("After downsampling: %d cells\n", nrow(cells_D)))

  PB <- pseudobulk(cells_D)
  res_D <- run_deseq2(PB, out_dir, "D")
  strip_plots(res_D, "D")
}

cat("\n=== DONE ===\nDESeq2 tables:", file.path(OUT_ROOT, APPROACH_DIRS[RUN], RUN_DATE), "\n")
