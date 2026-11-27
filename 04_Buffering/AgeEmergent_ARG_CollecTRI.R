################################################################################
### AgeEmergent_ARG_CollecTRI.R
###
### CollecTRI regulon over-representation among age-emergent ARGs, standalone.
### Same gene sets as ED8a (stacked bar) and ED8b-d (GO):
###   age-emergent = padj < 0.05 in MEEMHC / OEEOHC, and NOT padj < 0.05 in YEEYHC (or not tested there).
###   No fold-change threshold.
### Same test as the Figure 3 pipeline "tf" step (Fisher method):
###   per sex x age x cell type x direction (Up / Down / All): one-sided Fisher's exact test of each TF's regulon
###   (CollecTRI targets, unsigned) against all genes tested in the aged NE vs HC comparison for that cell type;
###   regulon needs >= MIN_TARGETS (5) targets in that universe, gene set >= MIN_GENES genes; BH across TFs within each test.
###   Effect size = log2 fold enrichment = log2((observed + 0.5) / (expected + 0.5)), as in Fig. 3d / ED7e.
### Across cell types (volcano): x = mean log2 fold enrichment over the cell types tested, y = -log10(lowest FDR across cell types),
###   size = number of cell types with FDR < FDR_THR.
### CollecTRI network prepared as in the buffering / Figure 6 analyses: TFs and targets named by a UniProt accession
###   removed, duplicate source-target edges removed, >= 5 targets (decoupleR default minsize).
###
### Reads DESeq2_Results only. Each sex analysed on its own.
### Output: Paper_Figures/Figure5/ED8ef_AgeEmergent_CollecTRI/
###   AgeEmergent_GeneSets.csv              every age-emergent gene (sex, age, cell type, direction)
###   Fisher_PerCellType.csv                every test (TF x cell type x sex x age x direction)
###   Fisher_AcrossCellTypes.csv            volcano table (mean log2 OR, lowest FDR, n cell types FDR < FDR_THR)
###   ED8ef_Volcano_<dir>.pdf/.png          males | females, middle (top) and old (bottom); <dir> = All, Up, Down
### Existing figures are skipped (SKIP_EXISTING); results are recomputed every run (fast).
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(conflicted)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::mutate, .quiet=TRUE); conflicts_prefer(base::intersect, .quiet=TRUE)
  conflicts_prefer(base::setdiff, .quiet=TRUE); conflicts_prefer(base::union, .quiet=TRUE)
})

################################################################################
### SETTINGS
################################################################################

BASE          <- "PATH_TO_NE_OUTPUT"
DESEQ_DIR     <- file.path(BASE, "DESeq2_Results")
COLLECTRI_RDS <- "PATH_TO_COLLECTRI/collectri_mouse.rds"
OUT           <- file.path(BASE, "Paper_Figures", "Figure5", "ED8ef_AgeEmergent_CollecTRI")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

SEXES       <- c("male", "female")
AGES        <- c(Middle = "MEEMHC", Old = "OEEOHC")
PADJ_THR    <- 0.05      # age-emergent gene definition (both comparisons)
FDR_THR     <- 0.05      # TF significance: dot size, dashed line, labels (same cutoff as the gene level)
MIN_GENES   <- 5         # minimum genes in a set to test
MIN_TARGETS <- 5         # minimum regulon targets in the universe (= ULM minsize in the buffering analyses)
DIRS_PLOT   <- c("All", "Up", "Down")   # one volcano figure per direction ("All" = Up + Down, as in the figure)
N_LABEL     <- 12        # TFs labelled per panel (lowest FDR among significant TFs)
SKIP_EXISTING <- TRUE

CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN",
              "Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell","Macrophage","Epithelial","T Cell")
COL_SIG <- "#C0442B"; COL_NS <- "grey80"

save_both <- function(p, stem, w, h) {
  if (SKIP_EXISTING && file.exists(paste0(stem, ".pdf"))) { cat("  exists, skipped:", basename(stem), "\n"); return(invisible()) }
  ggsave(paste0(stem, ".pdf"), p, width=w, height=h, limitsize=FALSE)
  ggsave(paste0(stem, ".png"), p, width=w, height=h, dpi=600, bg="white", limitsize=FALSE)
  cat("  saved:", basename(stem), "\n")
}

################################################################################
### DATA
################################################################################

read_deseq <- function(sx, comp) {
  f <- file.path(DESEQ_DIR, paste0(sx, "_ALL_", comp, ".csv"))
  if (!file.exists(f)) stop("Not found: ", f)
  read.csv(f, stringsAsFactors=FALSE)
}
if (!file.exists(COLLECTRI_RDS)) stop("CollecTRI RDS not found: ", COLLECTRI_RDS)
uniprot <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"   # as in the buffering scripts
net <- readRDS(COLLECTRI_RDS) %>%
  dplyr::filter(!grepl(uniprot, source), !grepl(uniprot, target)) %>%
  dplyr::distinct(source, target, .keep_all=TRUE)
regulons <- split(net$target, net$source)
cat(sprintf("CollecTRI (UniProt-ID TFs/targets removed): %d TFs, %d interactions\n", length(regulons), nrow(net)))
tf_label <- function(s) s

################################################################################
### GENE SETS + FISHER
################################################################################

fisher_one <- function(universe, genes) {
  genes <- intersect(genes, universe); k <- length(genes); N <- length(universe)
  if (k < MIN_GENES) return(NULL)
  dplyr::bind_rows(lapply(names(regulons), function(tf) {
    reg <- intersect(regulons[[tf]], universe); n_reg <- length(reg)
    if (n_reg < MIN_TARGETS) return(NULL)
    a <- sum(reg %in% genes); b <- n_reg - a; c <- k - a; d <- N - k - b
    p <- if (a == 0) 1 else stats::fisher.test(matrix(c(a, c, b, d), 2), alternative="greater")$p.value
    data.frame(source=tf, n_regulon=n_reg, n_in_set=a, n_set=k, n_universe=N,
               expected=n_reg * k / N, log2_fe=log2((a + 0.5) / (n_reg * k / N + 0.5)), p=p)
  }))
}

sets_out <- list(); res <- list()
for (sx in SEXES) {
  yee <- read_deseq(sx, "YEEYHC")
  for (ag in names(AGES)) {
    agd <- read_deseq(sx, AGES[[ag]])
    cts <- c(intersect(CT_ORDER, unique(agd$celltype)), sort(setdiff(unique(agd$celltype), CT_ORDER)))
    for (ct in cts) {
      y <- yee[yee$celltype == ct, ]; d <- agd[agd$celltype == ct & !is.na(agd$log2FoldChange), ]
      y_sig <- unique(y$gene[!is.na(y$padj) & y$padj < PADJ_THR])                 # significant in young (any LFC)
      up <- setdiff(unique(d$gene[!is.na(d$padj) & d$padj < PADJ_THR & d$log2FoldChange > 0]), y_sig)
      dn <- setdiff(unique(d$gene[!is.na(d$padj) & d$padj < PADJ_THR & d$log2FoldChange < 0]), y_sig)
      universe <- unique(d$gene)                                                    # genes tested in the aged comparison
      sets <- list(Up=up, Down=dn, All=union(up, dn))
      if (length(sets$All))
        sets_out[[paste(sx, ag, ct)]] <- data.frame(sex=sx, age=ag, celltype=ct, gene=c(up, dn),
                                                    direction=rep(c("Up", "Down"), c(length(up), length(dn))))
      for (dr in names(sets)) {
        f <- fisher_one(universe, sets[[dr]]); if (is.null(f) || !nrow(f)) next
        res[[paste(sx, ag, ct, dr)]] <- f %>% dplyr::mutate(sex=sx, age=ag, celltype=ct, direction=dr,
                                                           fdr=p.adjust(p, "BH"), .before=1)
      }
      cat(sprintf("  %-6s %-6s %-20s age-emergent up=%d down=%d\n", sx, ag, ct, length(up), length(dn)))
    }
  }
}
write.csv(dplyr::bind_rows(sets_out), file.path(OUT, "AgeEmergent_GeneSets.csv"), row.names=FALSE)
per_ct <- dplyr::bind_rows(res)
if (!nrow(per_ct)) stop("No gene set reached MIN_GENES")
per_ct <- per_ct %>% dplyr::mutate(label=tf_label(source))
write.csv(per_ct, file.path(OUT, "Fisher_PerCellType.csv"), row.names=FALSE)

across <- per_ct %>% dplyr::group_by(sex, age, direction, source, label) %>%
  dplyr::summarise(mean_log2_fe=mean(log2_fe), n_celltypes_tested=dplyr::n(),
                   n_sig=sum(fdr < FDR_THR), best_fdr=min(fdr),
                   sig_celltypes=paste(celltype[fdr < FDR_THR], collapse=", "), .groups="drop") %>%
  dplyr::mutate(sig=best_fdr < FDR_THR) %>% dplyr::arrange(sex, age, direction, best_fdr)
write.csv(across, file.path(OUT, "Fisher_AcrossCellTypes.csv"), row.names=FALSE)

################################################################################
### VOLCANOES (males | females; middle top, old bottom)
################################################################################

# Figures depend on FDR_THR: if they were drawn with another threshold, redraw them
stamp <- file.path(OUT, "figures_drawn_with.txt"); stamp_now <- paste0("FDR_THR=", FDR_THR)
if (SKIP_EXISTING && !identical(if (file.exists(stamp)) readLines(stamp, warn=FALSE)[1] else "", stamp_now)) {
  cat("  FDR_THR changed (or not recorded): redrawing the volcanoes\n"); SKIP_EXISTING <- FALSE }
writeLines(stamp_now, stamp)

has_repel <- requireNamespace("ggrepel", quietly=TRUE)
for (dr in DIRS_PLOT) {
  d <- across %>% dplyr::filter(direction == dr) %>%
    dplyr::mutate(y=-log10(pmax(best_fdr, 1e-300)),
                  sex_lab=factor(tools::toTitleCase(sex), levels=c("Male", "Female")),
                  age=factor(age, levels=names(AGES)))
  if (!nrow(d)) next
  lab <- d %>% dplyr::filter(sig) %>% dplyr::group_by(sex_lab, age) %>%
    dplyr::slice_min(best_fdr, n=N_LABEL, with_ties=FALSE) %>% dplyr::ungroup()
  p <- ggplot(d, aes(mean_log2_fe, y)) +
    geom_vline(xintercept=0, colour="grey75", linewidth=0.3) +
    geom_hline(yintercept=-log10(FDR_THR), linetype="dashed", colour="grey40", linewidth=0.3) +
    geom_point(data=d %>% dplyr::filter(!sig), colour=COL_NS, size=0.8) +
    geom_point(data=d %>% dplyr::filter(sig), aes(size=n_sig), colour=COL_SIG, alpha=0.9) +
    scale_size_area(max_size=4, name=paste0("Cell types\nFDR < ", FDR_THR),
                    breaks=function(l) { b <- unique(round(pretty(l))); b[b >= l[1] & b <= l[2]] }) +
    facet_grid(age ~ sex_lab, scales="free_y") +
    theme_classic(base_size=10) +
    theme(strip.background=element_blank(), strip.text=element_text(face="bold", size=11)) +
    labs(x="Mean log2 fold enrichment across cell types", y=expression(-log[10]~"(lowest FDR across cell types)"),
         title=paste0("CollecTRI regulons enriched in age-emergent ARGs (", dr, ")"))
  if (has_repel) p <- p + ggrepel::geom_text_repel(data=lab, aes(label=label), size=2.6, max.overlaps=Inf,
                                                   segment.size=0.2, min.segment.length=0, seed=1)
  save_both(p, file.path(OUT, paste0("ED8ef_Volcano_", dr)), 8, 7)
}
cat("\nDone:", OUT, "\n")
