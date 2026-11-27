################################################################################
### Fig4g_NE_DREADD_ConservedCompression.R
###
### Fig. 4g: genes induced by both eNE and DREADD in young females whose response is compressed with age in both.
###   Young ARG          eNE YEEYHC and DREADD YC21YPBS padj < 0.05, same direction, in the same cell type
###   Retained response  r = log2FC old / log2FC young, in each paradigm (eNE OEEOHC, DREADD OC21OPBS)
###   Compressed in both r < COMP_R (0.5) in eNE and in DREADD (eNE middle age, MEEMHC, shown but not used)
###   Gene level         one value per gene: log2FC averaged over the cell types where it is compressed in both;
###                      ranked by min(young response eNE, DREADD) x min(drop eNE, DREADD) (drop = 1 - r, capped 0..1)
###   Fig. 4g = heatmap of the TOP_GENES (20) top-ranked induced genes: eNE young / middle / old, DREADD young / old,
###             with the cell types where each gene is compressed in both (right half of the I_ figure; the left half is
###             the median (IQR) trajectory of all common genes)
### eNE females vs DREADD (female only); APPROACHES: D = trimming + downsampling (Fig. 4g), B = covariate.
###
### Inputs:  eNE  <NE_DESEQ>/female_ALL_{YEEYHC,MEEMHC,OEEOHC}.csv     (01_DESeq2/Fig2_ED4_Pipeline.R)
###          DREADD <DESEQ_ROOT>/DESeq2_mScarletPos_Approach<X>_MergedIN/<date>/<date>_{YC21YPBS,OC21OPBS}_ALL.csv (latest)
### Outputs: <OUT_DIR>/ConservedYoungARGs_AllPatterns.csv, ConservedCompressed_Genes.csv,
###          ConservedCompressed_Superpooled_GeneLevel.csv, I_Superpooled_TrajectoryHeatmap_NE<Sex>_Approach<X>_<UP|DOWN>.pdf/.png
### Run:  Rscript Fig4g_NE_DREADD_ConservedCompression.R     (set the paths below first)
################################################################################

suppressPackageStartupMessages({
  library(conflicted); library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  conflicts_prefer(dplyr::select, .quiet = TRUE); conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicts_prefer(dplyr::mutate, .quiet = TRUE); conflicts_prefer(base::intersect, .quiet = TRUE)
  conflicts_prefer(base::setdiff, .quiet = TRUE); conflicts_prefer(base::union, .quiet = TRUE)
})

################################################################################
### SETTINGS
################################################################################

DESEQ_ROOT <- "PATH_TO_DREADD_OUTPUT"                    # same as OUT_ROOT in DREADD_DESeq2.R
NE_DESEQ   <- "PATH_TO_NE_OUTPUT/DESeq2_Results"         # eNE DESeq2 tables (01_DESeq2)
OUT_DIR    <- file.path(DESEQ_ROOT, "Paper_Figures", "Figure4", "NE_DREADD_ConservedCompression")
APPROACHES <- c("D")                                     # D for Fig. 4g; c("B", "D") for both
DEG_DIRS   <- c(B = "DESeq2_mScarletPos_ApproachB_MergedIN", D = "DESeq2_mScarletPos_ApproachD_MergedIN")
NE_SEXES   <- c("female")                                # DREADD is female only
ARG_PADJ   <- 0.05
COMP_R     <- 0.5        # compressed: old response < half of young
MAINT_R    <- 0.75       # maintained: old response >= 75% of young (pattern table only)
TOP_GENES  <- 20
CT_ORDER   <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN","Astrocytes","Oligodendrocytes","OPC",
                "Microglia","Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell")
CT_SHORT   <- c("DG ExNs"="DG","CA3 ExNs"="CA3","CA1 ExNs"="CA1","CA2 ExNs"="CA2","Sub ExNs"="Sub","IN"="IN","Astrocytes"="Astro",
                "Oligodendrocytes"="Oligo","OPC"="OPC","Microglia"="Micro","Endothelial"="Endo","Pericyte"="Peri","APC"="APC",
                "Stromal"="Stromal","Cajal-Retzius Cell"="CR")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

################################################################################
### READ TABLES
################################################################################

read_deg <- function(ap, comp) {
  f <- list.files(file.path(DESEQ_ROOT, DEG_DIRS[[ap]]), pattern = paste0("_", comp, "_ALL\\.csv$"), recursive = TRUE, full.names = TRUE)
  if (!length(f)) stop("DREADD table not found: Approach ", ap, " ", comp)
  f <- sort(f, decreasing = TRUE)[1]; cat(sprintf("  Approach %s %s: %s\n", ap, comp, f))
  read.csv(f, stringsAsFactors = FALSE) %>% dplyr::transmute(celltype, gene, lfc = log2FoldChange, pvalue, padj)
}
rd_ne <- function(sx, comp) {
  f <- file.path(NE_DESEQ, paste0(sx, "_ALL_", comp, ".csv"))
  if (!file.exists(f)) stop("eNE table not found: ", f)
  read.csv(f, stringsAsFactors = FALSE) %>% dplyr::transmute(celltype, gene, lfc = log2FoldChange, pvalue, padj)
}
ne_w <- function(sx) {
  y <- rd_ne(sx, "YEEYHC"); o <- rd_ne(sx, "OEEOHC"); m <- rd_ne(sx, "MEEMHC")
  dplyr::full_join(y %>% dplyr::rename(lfc_y = lfc, p_y = pvalue, padj_y = padj),
                   o %>% dplyr::rename(lfc_o = lfc, p_o = pvalue, padj_o = padj), by = c("celltype", "gene")) %>%
    dplyr::left_join(m %>% dplyr::select(celltype, gene, lfc_m = lfc), by = c("celltype", "gene")) %>%
    dplyr::filter(celltype %in% CT_ORDER)
}
dr_w <- function(ap) {
  dy <- read_deg(ap, "YC21YPBS"); do <- read_deg(ap, "OC21OPBS")
  dplyr::full_join(dy %>% dplyr::rename(lfc_y = lfc, p_y = pvalue, padj_y = padj),
                   do %>% dplyr::rename(lfc_o = lfc, p_o = pvalue, padj_o = padj), by = c("celltype", "gene")) %>%
    dplyr::filter(celltype %in% CT_ORDER)
}
classify <- function(r) dplyr::case_when(is.na(r) ~ NA_character_, r >= MAINT_R ~ "maintained", r >= COMP_R ~ "partial",
                                         r > 0 ~ "compressed", TRUE ~ "lost/reversed")
prep <- function(d) d %>% dplyr::mutate(arg_y = !is.na(padj_y) & padj_y < ARG_PADJ & !is.na(lfc_y),
                                        r = ifelse(arg_y & !is.na(lfc_o), lfc_o / lfc_y, NA_real_), class = classify(r))

################################################################################
### SHARED YOUNG ARGS AND PATTERNS
################################################################################

G <- list()
for (sx in NE_SEXES) {
  ne <- prep(ne_w(sx))
  for (ap in APPROACHES) {
    G[[paste(sx, ap)]] <- dplyr::inner_join(ne, prep(dr_w(ap)), by = c("celltype", "gene"), suffix = c("_ne", "_dr")) %>%
      dplyr::mutate(ne_sex = sx, approach = ap,
                    conserved_arg = arg_y_ne & arg_y_dr & sign(lfc_y_ne) == sign(lfc_y_dr),
                    direction = ifelse(conserved_arg, ifelse(lfc_y_ne > 0, "up", "down"), NA_character_),
                    pattern = dplyr::case_when(!conserved_arg ~ NA_character_,
                                               !is.na(r_ne) & !is.na(r_dr) & r_ne < COMP_R & r_dr < COMP_R ~ "compressed in both",
                                               !is.na(r_ne) & !is.na(r_dr) & r_ne >= MAINT_R & r_dr >= MAINT_R ~ "maintained in both",
                                               !is.na(r_ne) & r_ne < COMP_R ~ "compressed in NE only",
                                               !is.na(r_dr) & r_dr < COMP_R ~ "compressed in DREADD only", TRUE ~ "other"),
                    score = pmin(abs(lfc_y_ne), abs(lfc_y_dr)) * pmin(1 - r_ne, 1 - r_dr))
  }
}
G <- dplyr::bind_rows(G)
write.csv(G %>% dplyr::filter(conserved_arg), file.path(OUT_DIR, "ConservedYoungARGs_AllPatterns.csv"), row.names = FALSE)
write.csv(G %>% dplyr::filter(pattern == "compressed in both") %>% dplyr::arrange(ne_sex, approach, direction, dplyr::desc(score)),
          file.path(OUT_DIR, "ConservedCompressed_Genes.csv"), row.names = FALSE)
cat("\nShared young ARGs by pattern:\n")
print(as.data.frame(G %>% dplyr::filter(conserved_arg) %>% dplyr::count(ne_sex, approach, direction, pattern) %>%
                      tidyr::pivot_wider(names_from = pattern, values_from = n, values_fill = 0)), row.names = FALSE)

################################################################################
### SUPERPOOLED GENE LEVEL + Fig. 4g HEATMAP
################################################################################

sp_list <- list()
for (sx in NE_SEXES) for (ap in APPROACHES) {
  cc <- G %>% dplyr::filter(ne_sex == sx, approach == ap, pattern == "compressed in both")
  if (!nrow(cc)) next
  gl <- cc %>% dplyr::group_by(gene) %>%
    dplyr::summarise(n_ct = dplyr::n(), cell_types = paste(celltype, collapse = ", "), direction = names(which.max(table(direction))),
                     NE_Young = mean(lfc_y_ne * sign(lfc_y_ne)), NE_Middle = mean(lfc_m * sign(lfc_y_ne), na.rm = TRUE),
                     NE_Old = mean(lfc_o_ne * sign(lfc_y_ne)), DR_Young = mean(lfc_y_dr * sign(lfc_y_dr)),
                     DR_Old = mean(lfc_o_dr * sign(lfc_y_dr)), .groups = "drop") %>%
    # rank on response size: strong young response in BOTH paradigms x large drop in BOTH (retained response capped 0..1)
    dplyr::mutate(drop_NE = 1 - pmin(pmax(NE_Old / NE_Young, 0), 1), drop_DR = 1 - pmin(pmax(DR_Old / DR_Young, 0), 1),
                  rank_score = pmin(NE_Young, DR_Young) * pmin(drop_NE, drop_DR)) %>%
    dplyr::arrange(dplyr::desc(rank_score)) %>%
    # back to real signed log2FC for plotting / saving (repressed genes negative)
    dplyr::mutate(dplyr::across(c(NE_Young, NE_Middle, NE_Old, DR_Young, DR_Old), ~ ifelse(direction == "down", -.x, .x)))
  gl$ct_short <- sapply(strsplit(gl$cell_types, ", "), function(v) {
    v <- v[order(match(v, CT_ORDER))]; paste(ifelse(v %in% names(CT_SHORT), CT_SHORT[v], v), collapse = ", ") })
  sp_list[[paste(sx, ap)]] <- gl %>% dplyr::mutate(ne_sex = sx, approach = ap, .before = 1)
  dr_lab <- paste0("DREADD (Approach ", ap, ")")

  for (dirn in c("up", "down")) {
    dl <- ifelse(dirn == "up", "induced", "repressed")
    gd <- gl %>% dplyr::filter(direction == dirn)
    top_g <- head(gd$gene, TOP_GENES)
    if (!length(top_g)) next
    long <- dplyr::bind_rows(
      gd %>% dplyr::transmute(gene, paradigm = "NE", Young = NE_Young, Middle = NE_Middle, Old = NE_Old),
      gd %>% dplyr::transmute(gene, paradigm = dr_lab, Young = DR_Young, Middle = NA_real_, Old = DR_Old)) %>%
      tidyr::pivot_longer(c(Young, Middle, Old), names_to = "age", values_to = "value") %>% dplyr::filter(!is.na(value)) %>%
      dplyr::mutate(age = factor(age, levels = c("Young", "Middle", "Old")), paradigm = factor(paradigm, levels = c("NE", dr_lab)))
    sm <- long %>% dplyr::group_by(paradigm, age) %>%
      dplyr::summarise(med = median(value), q1 = quantile(value, 0.25), q3 = quantile(value, 0.75), .groups = "drop")
    show_lines <- if (dplyr::n_distinct(long$gene) <= 400) long else long[0, ]
    p_tr <- ggplot(sm, aes(age, med, color = paradigm, group = paradigm)) + geom_hline(yintercept = 0, color = "grey70") +
      geom_line(data = show_lines, aes(age, value, group = interaction(gene, paradigm), color = paradigm),
                alpha = 0.06, linewidth = 0.3, inherit.aes = FALSE) +
      geom_ribbon(aes(ymin = q1, ymax = q3, fill = paradigm), alpha = 0.15, color = NA) + geom_line(linewidth = 1.2) + geom_point(size = 2.4) +
      scale_x_discrete(drop = FALSE) + scale_color_manual(values = c("#E6550D", "#3182BD"), name = NULL) +
      scale_fill_manual(values = c("#E6550D", "#3182BD"), guide = "none") + theme_classic(base_size = 9) + theme(legend.position = "top") +
      labs(title = sprintf("All %d common %s genes", dplyr::n_distinct(gd$gene), dl), x = NULL, y = "log2 FC (activity vs control)")
    hm <- gd %>% head(TOP_GENES) %>%
      dplyr::transmute(gene, ct_short, `NE young` = NE_Young, `NE middle` = NE_Middle, `NE old` = NE_Old,
                       `DREADD young` = DR_Young, `DREADD old` = DR_Old) %>%
      tidyr::pivot_longer(-c(gene, ct_short), names_to = "col", values_to = "v") %>%
      dplyr::mutate(col = factor(col, levels = c("NE young", "NE middle", "NE old", "DREADD young", "DREADD old")),
                    gene = factor(gene, levels = rev(top_g)))
    lim <- max(abs(hm$v), na.rm = TRUE)
    p_hm <- ggplot(hm, aes(col, gene, fill = v)) + geom_tile(color = "white", linewidth = 0.4) +
      geom_text(aes(label = sprintf("%.1f", v)), size = 2, color = "grey15") +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-lim, lim), name = "log2 FC") +
      geom_text(data = hm %>% dplyr::distinct(gene, ct_short), aes(x = 5.65, y = gene, label = ct_short), inherit.aes = FALSE,
                size = 2.2, hjust = 0, color = "grey30") +
      annotate("text", x = 5.65, y = Inf, label = "Cell types", hjust = 0, vjust = -0.6, size = 2.4, fontface = "bold", color = "grey30") +
      scale_x_discrete(expand = expansion(add = c(0.5, 2.2))) + coord_cartesian(clip = "off") + theme_classic(base_size = 8) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.text.y = element_text(face = "italic"), axis.line = element_blank(),
            axis.ticks = element_blank()) +
      labs(title = sprintf("Top %d by rank", length(top_g)), x = NULL, y = NULL)
    p <- p_tr + p_hm + patchwork::plot_layout(widths = c(1, 1.1)) +
      patchwork::plot_annotation(title = sprintf("Super-pooled %s genes compressed in both (NE %s, DREADD Approach %s)", dl, sx, ap),
        subtitle = "Left: median (IQR) trajectory of all common genes. Right: top genes ranked by strong young response x large drop in both paradigms (cell types where each gene is compressed in both on the right)",
        theme = theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 8, color = "grey40")))
    f <- file.path(OUT_DIR, sprintf("I_Superpooled_TrajectoryHeatmap_NE%s_Approach%s_%s", tools::toTitleCase(sx), ap, toupper(dirn)))
    h <- max(5, 1.5 + 0.2 * length(top_g))
    ggsave(paste0(f, ".pdf"), p, width = 13.5, height = h)
    ggsave(paste0(f, ".png"), p, width = 13.5, height = h, dpi = 600, bg = "white")
    cat("  saved:", basename(f), "\n")
  }
}
if (length(sp_list)) {
  sp <- dplyr::bind_rows(sp_list)
  write.csv(sp, file.path(OUT_DIR, "ConservedCompressed_Superpooled_GeneLevel.csv"), row.names = FALSE)
  cat("\nTop genes (Fig. 4g):\n")
  print(as.data.frame(sp %>% dplyr::filter(direction == "up") %>% dplyr::group_by(ne_sex, approach) %>%
                        dplyr::slice_head(n = TOP_GENES) %>% dplyr::select(ne_sex, approach, gene, ct_short)), row.names = FALSE)
}

cat("\n=== DONE ===\nOutputs:", OUT_DIR, "\n")
