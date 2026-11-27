################################################################################
### ED7d_Fig3d_ED7e_BasalDerepression.R
###
### Extended Data Fig. 7d (failed induction vs basal derepression among lost young NE genes) and Fig. 3d /
### Extended Data Fig. 7e (CollecTRI regulon enrichment in the failed and basal gene sets).
### Reads only the per-sex DESeq2 tables written by 01_DESeq2/Fig2_ED4_Pipeline.R (DESeq2_Results/); no Seurat object,
### DESeq2 is not rerun.
### Run:  Rscript ED7d_Fig3d_ED7e_BasalDerepression.R     (set BASE and COLLECTRI_RDS below first)
###
### Notation (per sex x cell type x age; "aged" = middle or old):
###   Y  = YEEYHC           young NE vs young HC  (the young NE program; anchor)
###   A  = MEEMHC / OEEOHC  aged NE vs aged HC    (NE response at the older age)
###   B  = MHCYHC / OHCYHC  aged HC vs young HC   (resting change with age)
###   young gene = BH padj < 0.05 in Y; direction = sign of its Y log2FC
###   lost       = young gene, not padj < 0.05 in A (either direction)
###   basal      = lost AND resting fraction (B log2FC / Y log2FC) >= BASAL_FRAC (0.5): the resting level already moved
###                at least half as far as the young NE response, in the same direction (no significance required on B)
###   failed     = lost AND not basal (failed induction = Up, failed repression = Down)
###
### Steps (RUN) and outputs (Paper_Figures/Figure3/BasalDerepression/<step>/):
###   table    GeneTable_<sex>_<age>.csv: every young gene x cell type with Y / A / B log2FC + padj, lost and basal calls
###   counts   Counts_D5.pdf (ED7d): mirrored stacked bar, males left, females right, Middle | Old;
###            basal (greens) vs failed (reds), dark = Up, light = Down
###   tf       CollecTRI regulon over-representation (one-sided Fisher's exact / hypergeometric) in the basal and failed
###            sets per sex x age x cell type x direction (UP / DOWN / All), vs all genes tested in YEEYHC for that
###            cell type; regulons with >= 3 targets in the universe, sets with >= TF_MIN_GENES genes;
###            effect size = log2((observed + 0.5) / (expected + 0.5)); BH across TFs within each test.
###            Across cell types: x = mean log2 fold enrichment, y = -log10(lowest FDR), size = number of cell types
###            with FDR < TF_FDR_THR.  tf/FDR0.05/Fisher/Fisher_BasalVsFailed_D5_<sex>_<age>_All.pdf
###            (Fig. 3d: females; ED7e: males)
### Each sex is analysed on its own. Existing files are skipped.
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

RUN <- c("table", "counts", "tf")
BASE          <- "PATH_TO_OUTPUT"                                # same BASE as 01_DESeq2/Fig2_ED4_Pipeline.R
DESEQ_DIR     <- file.path(BASE, "DESeq2_Results")               # <sex>_ALL_<comparison>.csv
OUT_ROOT      <- file.path(BASE, "Paper_Figures", "Figure3", "BasalDerepression")
COLLECTRI_RDS <- "PATH_TO_collectri_mouse.rds"                   # CollecTRI mouse network: data frame with source, target, mor
dir.create(OUT_ROOT, recursive=TRUE, showWarnings=FALSE)

SEXES      <- c("male", "female")
PADJ_THR   <- 0.05      # BH cutoff for young genes and for "lost"
BASAL_FRAC <- 0.5       # basal: resting change >= this fraction of the young NE log2FC (same sign)
SKIP_EXISTING <- TRUE
# -- tf step (CollecTRI) ---------------------------------------------------------------------------------------------
TF_FDR_THR    <- 0.05   # TF significance: dashed line, dot size, labels
TF_MIN_GENES  <- 5      # minimum genes per set
TF_LABEL_N    <- 8      # TFs labelled per panel

AGES <- list(Middle = c(A = "MEEMHC", B = "MHCYHC"),
             Old    = c(A = "OEEOHC", B = "OHCYHC"))
CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN","Astrocytes","Oligodendrocytes","OPC",
              "Microglia","Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell","Macrophage","Epithelial","T Cell")
D5_LABEL <- sprintf("D5 Lost + rest change >= %g x young", BASAL_FRAC)
D5_RULE  <- sprintf("lost & resting log2FC >= %g x young log2FC", BASAL_FRAC)

# -- helpers -------------------------------------------------------------------------------------------------------
LOG <- file.path(OUT_ROOT, paste0("BasalDerepression_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
say <- function(...) { m <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", ...); cat(m, "\n"); cat(m, "\n", file=LOG, append=TRUE) }

ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
  if (SKIP_EXISTING && file.exists(filename)) { cat("  exists, skipped:", basename(filename), "\n"); return(invisible(filename)) }
  ggplot2::ggsave(filename, plot, ...); cat("  saved:", basename(filename), "\n"); invisible(filename)
}
write.csv <- function(x, file = "", ...) {
  if (SKIP_EXISTING && is.character(file) && nzchar(file) && file.exists(file)) { cat("  exists, skipped:", basename(file), "\n"); return(invisible(NULL)) }
  utils::write.csv(x, file, ...)
}
save_both <- function(p, stem, w, h) {
  p <- p + labs(caption = paste0(if (!is.null(p$labels$caption)) paste0(p$labels$caption, "\n") else "",
                                 "Significance: BH padj < ", PADJ_THR, " throughout"))
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h, limitsize = FALSE)
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 300, bg = "white", limitsize = FALSE)
}
step_dir <- function(s) { d <- file.path(OUT_ROOT, s); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
ct_levels <- function(cts) c(intersect(CT_ORDER, cts), sort(setdiff(cts, CT_ORDER)))
T_ <- function(x) !is.na(x) & x                       # NA-safe TRUE
sig <- function(p) !is.na(p) & p < PADJ_THR

################################################################################
### PRE-FLIGHT + GENE TABLE (every step uses it)
################################################################################

need <- c("YEEYHC", unlist(AGES))
miss <- setdiff(paste0(rep(SEXES, each = length(need)), "_ALL_", need, ".csv"), list.files(DESEQ_DIR))
if (length(miss)) stop("Missing DESeq2 tables in ", DESEQ_DIR, ": ", paste(miss, collapse = ", "))
say("Basal derepression (D5). DESeq2: ", DESEQ_DIR, "  Output: ", OUT_ROOT)

read_comp <- function(sx, comp, tag) {
  d <- read.csv(file.path(DESEQ_DIR, paste0(sx, "_ALL_", comp, ".csv")), stringsAsFactors = FALSE)
  d <- d[, c("celltype", "gene", "log2FoldChange", "padj")]
  names(d)[3:4] <- paste0(c("lfc_", "padj_"), tag); d
}
build_table <- function(sx, age) {
  cm <- AGES[[age]]
  m <- read_comp(sx, "YEEYHC", "Y") %>%
    dplyr::full_join(read_comp(sx, cm[["A"]], "A"), by = c("celltype", "gene")) %>%
    dplyr::full_join(read_comp(sx, cm[["B"]], "B"), by = c("celltype", "gene"))
  m %>% dplyr::mutate(
    sex = sx, age = age,
    young = sig(padj_Y), dir_young = ifelse(young, ifelse(lfc_Y > 0, "Up", "Down"), NA),
    lost  = young & !sig(padj_A),
    rest_frac = ifelse(young & !is.na(lfc_B) & lfc_Y != 0, lfc_B / lfc_Y, NA),   # > 0 = same direction
    D5    = lost & T_(rest_frac >= BASAL_FRAC))
}

say("Building gene tables at padj < ", PADJ_THR)
GT <- list()
for (sx in SEXES) for (age in names(AGES)) GT[[paste(sx, age)]] <- build_table(sx, age)
ALL <- dplyr::bind_rows(GT)
ALL <- dplyr::mutate(ALL, celltype = factor(celltype, levels = ct_levels(unique(as.character(celltype)))))

SECTIONS <- list()

################################################################################
### table: gene-level table
################################################################################
SECTIONS[["table"]] <- function() {
  OUT <- step_dir("table")
  for (k in names(GT)) {
    g <- GT[[k]] %>% dplyr::filter(young)
    write.csv(g, file.path(OUT, paste0("GeneTable_", sub(" ", "_", k), ".csv")), row.names = FALSE)
  }
  cnt <- ALL %>% dplyr::filter(lost) %>% dplyr::mutate(class = ifelse(D5, "basal", "failed")) %>%
    dplyr::count(sex, age, celltype, dir_young, class, name = "n_genes")
  write.csv(cnt, file.path(OUT, "Counts_D5.csv"), row.names = FALSE)
}

################################################################################
### counts: ED7d mirrored stacked bar
################################################################################

mirror_bar <- function(d, seg_order, seg_cols, title, xlab = "Number of genes") {
  cts <- rev(levels(droplevels(d$celltype)))
  r <- d %>% dplyr::filter(n_genes > 0) %>% dplyr::mutate(segment = factor(segment, levels = seg_order)) %>%
    dplyr::group_by(sex, age, celltype) %>% dplyr::arrange(segment, .by_group = TRUE) %>%
    dplyr::mutate(e = cumsum(n_genes), s = e - n_genes,
                  xmin = ifelse(sex == "male", -e, s), xmax = ifelse(sex == "male", -s, e)) %>% dplyr::ungroup() %>%
    dplyr::mutate(y = match(as.character(celltype), cts), age = factor(age, levels = names(AGES)))
  if (!nrow(r)) return(NULL)
  lim <- r %>% dplyr::group_by(age) %>% dplyr::summarise(m = max(e) * 1.25, .groups = "drop")
  tot <- r %>% dplyr::group_by(sex, age, y) %>% dplyr::summarise(t = max(e), .groups = "drop") %>% dplyr::left_join(lim, by = "age") %>%
    dplyr::mutate(x = ifelse(sex == "male", -t - m * 0.015, t + m * 0.015), h = ifelse(sex == "male", 1, 0))
  heads <- lim %>% tidyr::crossing(sex = SEXES) %>%
    dplyr::mutate(x = ifelse(sex == "male", -m * 0.5, m * 0.5), y = length(cts) + 0.9, lab = ifelse(sex == "male", "Males", "Females"))
  ggplot() +
    geom_blank(data = lim %>% tidyr::crossing(k = c(-1, 1)) %>% dplyr::mutate(x = k * m, y = 1), aes(x, y)) +
    geom_rect(data = r, aes(xmin = xmin, xmax = xmax, ymin = y - 0.38, ymax = y + 0.38, fill = segment)) +
    geom_text(data = tot, aes(x, y, label = t, hjust = h), size = 2.5) +
    geom_text(data = heads, aes(x, y, label = lab), fontface = "bold", size = 3.4) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = seg_cols, breaks = seg_order, name = NULL, drop = FALSE) +
    scale_x_continuous(labels = abs, expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(cts), labels = cts, expand = expansion(add = c(0.6, 1.3))) +
    facet_wrap(~age, nrow = 1, scales = "free_x") +
    theme_classic(base_size = 10) +
    theme(strip.text = element_text(face = "bold", size = 12), strip.background = element_blank(),
          axis.text = element_text(color = "black"), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3), plot.title = element_text(face = "bold", hjust = 0.5),
          panel.spacing = unit(1.2, "lines")) +
    labs(title = title, x = xlab, y = NULL)
}
SECTIONS[["counts"]] <- function() {
  OUT <- step_dir("counts")
  H <- max(5, length(unique(ALL$celltype)) * 0.3 + 1.6)
  segs <- c("Failed (Up)", "Failed (Down)", "Basal (Up)", "Basal (Down)")
  cols <- setNames(c("#B2402A", "#E06666", "#4E8B4E", "#97CC8E"), segs)
  d <- ALL %>% dplyr::filter(lost) %>%
    dplyr::mutate(segment = paste0(ifelse(D5, "Basal", "Failed"), " (", dir_young, ")")) %>%
    dplyr::count(sex, age, celltype, segment, name = "n_genes")
  p <- mirror_bar(d, segs, cols, paste0(D5_LABEL, "\n", D5_RULE))
  if (!is.null(p)) save_both(p, file.path(OUT, "Counts_D5"), 11, H)
}

################################################################################
### tf: CollecTRI regulon enrichment, basal vs failed (Fig. 3d, ED7e)
################################################################################
SECTIONS[["tf"]] <- function() {
  OUT <- step_dir("tf")
  if (!file.exists(COLLECTRI_RDS)) stop("CollecTRI RDS not found: ", COLLECTRI_RDS)

  net <- readRDS(COLLECTRI_RDS)
  regulons <- split(net$target, net$source)
  PROT_RE <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
  prot <- unique(net$source[grepl(PROT_RE, net$source)]); PMAP <- setNames(rep(NA_character_, length(prot)), prot)
  if (length(prot) && requireNamespace("org.Mm.eg.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
    m <- tryCatch(suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, keys = prot, column = "SYMBOL",
                                                         keytype = "UNIPROT", multiVals = "first")), error = function(e) NULL)
    if (!is.null(m)) PMAP[names(m)] <- unname(m)
  }
  disp <- function(s) { s <- as.character(s); i <- s %in% names(PMAP); s[i] <- ifelse(is.na(PMAP[s[i]]), s[i], paste0(PMAP[s[i]], "*")); s }

  # ---- Fisher's exact test per sex x age x cell type x class x direction ------------------------------------------
  fisher_set <- function(universe, genes) {
    genes <- intersect(genes, universe); k <- length(genes); N <- length(universe)
    if (k < TF_MIN_GENES) return(NULL)
    K <- vapply(regulons, function(r) sum(r %in% universe), 0); a <- vapply(regulons, function(r) sum(r %in% genes), 0)
    ok <- K >= 3; if (!any(ok)) return(NULL)
    K <- K[ok]; a <- a[ok]
    p <- stats::phyper(a - 1, K, N - K, k, lower.tail = FALSE)          # = one-sided Fisher's exact p
    E <- K * k / N                                                      # expected targets in the set
    data.frame(source = names(K), value = log2((a + 0.5) / (E + 0.5)), p_value = p, n_in_set = a,
               n_regulon = K, n_universe = N, expected = E)
  }
  res_csv <- file.path(OUT, "TF_Fisher_AllResults.csv")
  if (SKIP_EXISTING && file.exists(res_csv)) {
    say("  reading existing TF results (delete ", basename(res_csv), " to recompute)")
    RES <- read.csv(res_csv, stringsAsFactors = FALSE)
  } else {
    acc <- list()
    for (k in names(GT)) {
      g <- GT[[k]]
      for (ct in unique(g$celltype)) {
        gc <- g[g$celltype == ct, ]
        universe <- gc$gene[!is.na(gc$lfc_Y)]                            # all genes tested in YEEYHC for this cell type
        for (kind in c("basal", "failed")) {
          sel <- if (kind == "basal") T_(gc$D5) else (gc$lost & !T_(gc$D5))
          for (dr in c("UP", "DOWN", "All")) {
            genes <- gc$gene[sel & (dr == "All" | T_(gc$dir_young == ifelse(dr == "UP", "Up", "Down")))]
            if (length(genes) < TF_MIN_GENES) next
            f <- fisher_set(universe, genes)
            if (!is.null(f)) acc[[length(acc) + 1]] <- cbind(f, data.frame(sexage = k, celltype = as.character(ct), def = "D5",
                                                                           kind = kind, direction = dr, n_genes = length(genes)))
          }
        }
      }
      say("  TF tests done: ", k)
    }
    RES <- dplyr::bind_rows(acc) %>% dplyr::group_by(sexage, celltype, def, kind, direction) %>%
      dplyr::mutate(fdr = p.adjust(p_value, "BH")) %>% dplyr::ungroup()
    utils::write.csv(RES, res_csv, row.names = FALSE)
  }

  # ---- volcanoes -----------------------------------------------------------------------------------------------
  has_repel <- requireNamespace("ggrepel", quietly = TRUE)
  thr <- TF_FDR_THR
  MOUT <- file.path(OUT, paste0("FDR", thr), "Fisher"); dir.create(MOUT, recursive = TRUE, showWarnings = FALSE)
  pool <- RES %>% dplyr::group_by(sexage, def, kind, direction, source) %>%
    dplyr::summarise(mean_value = mean(value, na.rm = TRUE), n_sig = sum(fdr < thr, na.rm = TRUE),
                     best_fdr = min(fdr, na.rm = TRUE), sig_cts = paste(celltype[fdr < thr], collapse = ", "), .groups = "drop") %>%
    dplyr::mutate(y = -log10(pmax(best_fdr, 1e-300)), label = disp(source))
  write.csv(pool, file.path(MOUT, "Fisher_Pooled_AcrossCellTypes.csv"), row.names = FALSE)
  volcano <- function(d, title, facet_var, ncol) {
    lab <- d %>% dplyr::filter(n_sig > 0) %>% dplyr::group_by(.data[[facet_var]]) %>%
      dplyr::slice_max(n_sig * 100 + y, n = TF_LABEL_N, with_ties = FALSE) %>% dplyr::ungroup()
    ggplot(d, aes(x = mean_value, y = y)) +
      geom_hline(yintercept = -log10(thr), linetype = "dashed", color = "grey55", linewidth = 0.4) +
      geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4) +
      geom_point(data = d[d$n_sig == 0, ], color = "grey80", size = 0.9) +
      geom_point(data = d[d$n_sig > 0, ], aes(size = n_sig, color = kind), alpha = 0.8) +
      { if (nrow(lab)) { if (has_repel) ggrepel::geom_text_repel(data = lab, aes(label = label), size = 2.3, max.overlaps = 30, min.segment.length = 0)
                         else geom_text(data = lab, aes(label = label), size = 2.3, vjust = -0.6) } } +
      scale_color_manual(values = c(basal = "#4E8B4E", failed = "#B2402A"), guide = "none") +
      scale_size_continuous(range = c(1.2, 4.5), name = paste0("Cell types\nFDR < ", thr)) +
      facet_wrap(stats::as.formula(paste("~", facet_var)), ncol = ncol, scales = "free_y") +
      theme_classic(base_size = 9) +
      theme(strip.text = element_text(face = "bold", size = 8), strip.background = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, color = "grey35", size = 8)) +
      labs(title = title, subtitle = paste0("Regulon enrichment (Fisher's exact, CollecTRI); TF FDR < ", thr),
           x = "Mean log2 fold enrichment across cell types", y = expression(-log[10]~"(best FDR)"))
  }
  for (k in names(GT)) for (dr in c("All", "UP", "DOWN")) {
    d2 <- pool %>% dplyr::filter(sexage == k, direction == dr) %>%
      dplyr::mutate(kind = factor(kind, levels = c("failed", "basal")),
                    set = factor(ifelse(kind == "basal", "Basal (this rule)", "Failed (rest of lost genes)"),
                                 levels = c("Failed (rest of lost genes)", "Basal (this rule)")))
    if (nrow(d2)) save_both(volcano(d2, paste0(tools::toTitleCase(k), " (", dr, "): ", D5_LABEL), "set", 2),
                            file.path(MOUT, paste0("Fisher_BasalVsFailed_D5_", sub(" ", "_", k), "_", dr)), 9, 4.6)
  }
}

################################################################################
### RUN
################################################################################
status <- character(0)
for (id in intersect(c("table", "counts", "tf"), RUN)) {
  say("STEP ", id)
  status[id] <- tryCatch({ SECTIONS[[id]](); "ok" }, error = function(e) { say("  ERROR in ", id, ": ", conditionMessage(e)); "FAIL" })
}
say("Done: ", paste(names(status), status, sep = " ", collapse = " | "), "   Outputs: ", OUT_ROOT)
