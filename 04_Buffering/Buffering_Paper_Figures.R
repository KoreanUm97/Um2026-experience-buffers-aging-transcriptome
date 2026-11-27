################################################################################
### Buffering_Paper_Figures.R
###
### Panels for Figures 4-6 (one PDF + PNG per panel, plus source-data CSVs).
### No DESeq2 refits: reads existing outputs; only the bias-corrected buffering
### index (and, for DREADD, the bias-corrected cosine + simulated null) is
### computed here from the pseudobulk caches (minutes).
###
### FIGURE 4  NE (each sex separately; middle and old)
###   4a  Held-out cosine per cell type (style of the old figure: females = triangles,
###       males = circles; grey dots = each held-out HC mouse). Solid = passes the
###       gate (held-out cosine < 0 AND bias-corrected cosine beyond the simulated null),
###       faded = does not.                                         one panel per age
###   4b  Buffering index (bias-corrected; 1 = as far from young as HC; < 1 = closer)
###   4c  Quadrant: held-out cosine (x) vs buffering index (y), per age
###   4d  Heatmaps of buffered genes per gated cell type (z-scored expression per mouse:
###       young HC | aged HC | aged NE), T1 genes (T3 if too few)
###   ED  Bias-corrected cosine with the simulated no-effect band (gate, supplement)
###
### FIGURE 5  DREADD (female, old), Approach B and D separately, NE old female alongside
###   5a  Held-out cosine (grey band = C21-vs-C21 mouse-noise range, where available)
###   5b  Buffering index      5c  Quadrant (NE old female shown faded for reference)
###   ED  Bias-corrected cosine with simulated null
###
### FIGURE 6  NE-buffered genes and NE/DREADD moving together (old female)
###   6a  NE-buffered genes: aging LFC | NE LFC | DREADD LFC (per cell type)
###   6b  Gene-level agreement: NE effect vs DREADD effect, both oriented to NE's aging
###       direction (lower-left = both reverse aging, upper-right = both push further)
###   6c  Overlap odds ratios (Part 3) per cell type, Approach B and D
###   6d  TF direction agreement: % of NE-called TFs with the same sign in DREADD,
###       real labeling vs the relabelings
###   6e  TF buffering score NE vs DREADD (NE aging axis), shortlist labelled
###   6f  Shortlist TFs per cell type (Nr3c1, Creb1): NE vs DREADD
###
### Inputs
###   NE   Buffering_LOO_Split/pseudobulk_counts_and_meta.RDS
###        Buffering_MiddleAge_Robustness/*_Check1_BiasCorrected.csv, *_Check2_HeldOut_PerMouse.csv
###        Buffering_Candidates_Master/Step2_GeneCalls.RDS, Step5_TFActivity_PerCellType.csv
###   DREADD  DREADD_Buffering_HeldOut/pseudobulk_mScarlet_50_Approach_{B,D}.RDS
###           DREADD_Buffering_HeldOut/<date>/*_DREADD_HeldOut_PerFit.csv, *_Summary.csv
###           Buffering_DREADD_vs_NE/Approach_{B,D}/ (Part2/3/4 CSVs, units/)
### Run after Buffering_Pipeline.R steps 1-2 (robustness and master outputs).
################################################################################

suppressPackageStartupMessages({
  library(conflicted); library(DESeq2); library(dplyr); library(tidyr); library(ggplot2)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::mutate, .quiet=TRUE); conflicts_prefer(dplyr::slice_max, .quiet=TRUE)
  conflicts_prefer(base::intersect, .quiet=TRUE); conflicts_prefer(base::setdiff, .quiet=TRUE)
  conflicts_prefer(base::union, .quiet=TRUE)
})
HAS_REPEL <- requireNamespace("ggrepel", quietly=TRUE)

################################################################################
### PARAMETERS
################################################################################

NE_BASE <- "PATH_TO_NE_OUTPUT"
DR_BASE <- "PATH_TO_DREADD_OUTPUT"
NE_MASTER <- file.path(NE_BASE, "Buffering_Candidates_Master")
DR_VS_NE  <- file.path(DR_BASE, "Buffering_DREADD_vs_NE")
fig_date  <- format(Sys.Date(), "%Y-%m-%d")
OUT       <- file.path(NE_BASE, "Paper_Figures", fig_date)

RUN <- c(fig4=TRUE, fig5=TRUE, fig6=TRUE)
APPROACHES     <- c("B","D")
MIN_GENES      <- 1000          # cell types with fewer genes are left out of the plots
N_SIM          <- 50            # simulated null draws (DREADD bias-corrected cosine)
HEAT_TIER      <- "T1"          # heatmap genes; falls back to T3 if fewer than HEAT_MIN
HEAT_MIN       <- 10
HEAT_PER_DIR   <- 25            # max genes per direction (up / down with age) per heatmap
FIG6_CTS       <- c("DG ExNs","CA1 ExNs","OPC","Endothelial")   # cell types for 6a/6b
FIG6_N_GENES   <- 15            # genes per cell type in 6a
FIG6_TFS       <- c("Nr3c1","Creb1")                            # 6f panels
LABEL_TFS      <- c("Nr3c1","Creb1","Atf6","Sirt1")             # labelled in 6e
set.seed(20261022)

CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN",
              "Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell","Macrophage","T Cell")
CT_COLS <- c("DG ExNs"="#E53935","CA1 ExNs"="#CC0000","CA2 ExNs"="#8E24AA",
             "CA3 ExNs"="#5E35B1","Sub ExNs"="#80D8FF","IN"="#00897B",
             "Astrocytes"="#43A047","Oligodendrocytes"="#1E88E5","OPC"="#00ACC1",
             "Microglia"="#7CB342","Endothelial"="#FB8C00","Pericyte"="#6D4C41",
             "APC"="#546E7A","Stromal"="#27408B","Macrophage"="#C62828",
             "Cajal-Retzius Cell"="#FDD835","T Cell"="#757575")
CT_SHORT <- c("Astrocytes"="Astro","Oligodendrocytes"="Oligo","Cajal-Retzius Cell"="CR cell",
              "Endothelial"="Endo","Microglia"="Micro","Pericyte"="Peri")
short_ct <- function(x) ifelse(x %in% names(CT_SHORT), CT_SHORT[x], x)
COL_BUF <- "#1B7837"; COL_AMP <- "#762A83"; COL_OTH <- "grey80"
COL_AP  <- c(B="#D6604D", D="#4393C3")

th <- theme_classic(base_size=7) +
  theme(strip.background=element_rect(fill="grey92", color="black", linewidth=0.3),
        strip.text=element_text(face="bold", size=7), axis.line=element_line(linewidth=0.3),
        axis.ticks=element_line(linewidth=0.3), plot.title=element_text(face="bold", hjust=0.5, size=8),
        plot.subtitle=element_text(color="grey35", hjust=0.5, size=6), legend.key.size=unit(3, "mm"))

dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
mm <- function(x) x / 25.4
sv <- function(p, name, w_mm, h_mm) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width=mm(w_mm), height=mm(h_mm), limitsize=FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width=mm(w_mm), height=mm(h_mm), dpi=600, bg="white", limitsize=FALSE)
  cat("  saved:", name, "\n")
}
wsrc <- function(x, name) write.csv(x, file.path(OUT, paste0("SourceData_", name, ".csv")), row.names=FALSE)
latest <- function(dir, pat) {
  f <- list.files(dir, pattern=pat, recursive=TRUE, full.names=TRUE)
  if (!length(f)) stop("No file matching '", pat, "' under ", dir)
  f[order(file.mtime(f), decreasing=TRUE)][1]
}
ct_levels <- function(x) rev(intersect(CT_ORDER, unique(x)))

################################################################################
### SHARED: log-expression, bias-corrected cosine + buffering index
################################################################################

log_norm <- function(cnt) {
  cnt <- round(as.matrix(cnt)); cnt <- cnt[rowSums(cnt >= 10) >= 2, , drop=FALSE]
  sf  <- DESeq2::estimateSizeFactorsForMatrix(cnt)
  log2(sweep(cnt, 2, sf, "/") + 1)
}
# B = aged ctrl - young; E = treated - aged ctrl; N = treated - young (group means, log scale)
stats3 <- function(X, y, o, e) {
  mY <- rowMeans(X[, y, drop=FALSE]); mO <- rowMeans(X[, o, drop=FALSE]); mE <- rowMeans(X[, e, drop=FALSE])
  vY <- apply(X[, y, drop=FALSE], 1, var) / length(y)
  vO <- apply(X[, o, drop=FALSE], 1, var) / length(o)
  vE <- apply(X[, e, drop=FALSE], 1, var) / length(e)
  B <- mO - mY; E <- mE - mO; N <- mE - mY
  BB <- sum(B^2); BB_c <- BB - sum(vO + vY); NN_c <- sum(N^2) - sum(vE + vY)
  c(cos_raw = sum(B*E) / sqrt(BB * sum(E^2)),
    cos_dc  = (sum(B*E) + sum(vO)) / sqrt(BB * sum(E^2)),
    BI_raw  = sum(N^2) / BB,
    BI_corr = if (BB_c > 0) NN_c / BB_c else NA_real_)
}
sim_null <- function(X, y, o, e) {
  s  <- sqrt(rowMeans(cbind(apply(X[, y, drop=FALSE], 1, var), apply(X[, o, drop=FALSE], 1, var),
                            apply(X[, e, drop=FALSE], 1, var))))
  mu <- rowMeans(X)
  v  <- replicate(N_SIM, { Z <- matrix(rnorm(length(mu) * ncol(X), mu, s), nrow=length(mu), dimnames=dimnames(X))
                           stats3(Z, y, o, e)[["cos_dc"]] })
  c(null_mean=mean(v), null_sd=sd(v))
}

################################################################################
### SHARED PLOTS
################################################################################

# Held-out cosine dot plot. sm: one row per panel x celltype (cos_heldout, pass, shape_lab);
# pts: per held-out fit (cosine); band (optional): lo, hi per panel x celltype
plot_cos <- function(sm, pts, title, band=NULL, subtitle=NULL) {
  sm  <- sm %>% dplyr::mutate(gate=factor(ifelse(pass, "passes", "does not pass"), levels=c("passes","does not pass")))
  lim <- max(0.1, abs(c(sm$cos_heldout, pts$cosine, band$lo, band$hi)), na.rm=TRUE) * 1.05
  g <- ggplot(sm, aes(y=celltype)) + geom_vline(xintercept=0, linetype="dashed", color="grey40", linewidth=0.3)
  if (!is.null(band) && nrow(band))
    g <- g + geom_errorbar(data=band, aes(y=celltype, xmin=lo, xmax=hi), inherit.aes=FALSE,
                           width=0, color="grey88", linewidth=2.2, orientation="y")
  g + geom_point(data=pts, aes(x=cosine, y=celltype), inherit.aes=FALSE, color="grey55", size=0.6) +
    geom_point(aes(x=cos_heldout, fill=celltype, shape=shape_lab, alpha=gate), size=2.3, stroke=0.3, color="black") +
    scale_fill_manual(values=CT_COLS, guide="none") +
    scale_shape_manual(values=c(Females=24, Males=21), guide="none") +
    scale_alpha_manual(values=c(passes=1, `does not pass`=0.3), name="Buffering gate", drop=FALSE) +
    scale_x_continuous(limits=c(-lim, lim)) +
    facet_grid(. ~ panel) + th + theme(legend.position="bottom") +
    labs(title=title, subtitle=subtitle, x="Held-out cosine similarity (aging, stimulation)", y=NULL)
}
plot_bi <- function(d, title) {
  ggplot(d %>% dplyr::filter(!is.na(BI_corr)), aes(y=celltype)) +
    geom_vline(xintercept=1, linetype="dotted", color="grey30", linewidth=0.4) +
    geom_segment(aes(x=1, xend=BI_corr, yend=celltype), color="grey70", linewidth=0.5) +
    geom_point(aes(x=BI_corr, color=celltype, size=n_genes)) +
    scale_color_manual(values=CT_COLS, guide="none") +
    scale_size_continuous(range=c(0.8, 2.8), name="Genes tested") +
    facet_grid(. ~ panel, scales="free_x") + th + theme(legend.position="bottom") +
    labs(title=title, subtitle="< 1: stimulated closer to young than aged controls; > 1: further",
         x="Buffering index (bias-corrected)", y=NULL)
}
plot_quad <- function(d, title, ref=NULL, facet=NULL) {
  d <- d %>% dplyr::filter(!is.na(BI_corr), !is.na(cos_heldout)) %>% dplyr::mutate(lab=short_ct(as.character(celltype)))
  g <- ggplot(d, aes(cos_heldout, BI_corr)) +
    geom_hline(yintercept=1, linetype="dotted", color="grey30", linewidth=0.4) +
    geom_vline(xintercept=0, linetype="dashed", color="grey30", linewidth=0.3)
  if (!is.null(ref) && nrow(ref))
    g <- g + geom_point(data=ref, aes(cos_heldout, BI_corr), inherit.aes=FALSE, shape=24, fill="grey85",
                        color="grey60", size=1.6, stroke=0.2)
  g <- g + geom_point(aes(fill=celltype, shape=shape_lab), size=2.2, stroke=0.3, color="black") +
    scale_fill_manual(values=CT_COLS, guide="none") +
    scale_shape_manual(values=c(Females=24, Males=21), name=NULL) +
    annotate("text", x=-Inf, y=-Inf, hjust=-0.05, vjust=-0.4, size=1.9, color="grey35",
             label="Opposes aging,\ncloser to young") +
    annotate("text", x=Inf, y=Inf, hjust=1.05, vjust=1.3, size=1.9, color="grey35",
             label="Follows aging,\nfurther from young") +
    annotate("text", x=-Inf, y=Inf, hjust=-0.05, vjust=1.3, size=1.9, color="grey35",
             label="Opposes aging,\nfurther from young") +
    annotate("text", x=Inf, y=-Inf, hjust=1.05, vjust=-0.4, size=1.9, color="grey35",
             label="Follows aging,\ncloser to young") +
    th + theme(legend.position="bottom") +
    labs(title=title, x="Held-out cosine similarity", y="Buffering index (bias-corrected)")
  g <- g + if (HAS_REPEL) ggrepel::geom_text_repel(aes(label=lab), size=1.8, segment.size=0.2,
                                                   max.overlaps=30, min.segment.length=0.1)
           else geom_text(aes(label=lab), size=1.8, vjust=-0.9)
  if (!is.null(facet)) g <- g + facet_wrap(as.formula(paste("~", facet)))
  g
}
plot_dc <- function(d, title) {
  ggplot(d, aes(y=celltype)) +
    geom_vline(xintercept=0, linetype="dashed", color="grey40", linewidth=0.3) +
    geom_errorbar(aes(xmin=null_mean - 2*null_sd, xmax=null_mean + 2*null_sd), width=0,
                  color="grey88", linewidth=2.2, orientation="y") +
    geom_point(aes(x=cos_raw), shape=5, size=1.4, color="grey50", stroke=0.3) +
    geom_point(aes(x=cos_dc, fill=celltype, shape=shape_lab), size=2.2, stroke=0.3, color="black") +
    scale_fill_manual(values=CT_COLS, guide="none") +
    scale_shape_manual(values=c(Females=24, Males=21), guide="none") +
    facet_grid(. ~ panel) + th +
    labs(title=title, subtitle="Open diamond = raw coupled | filled = bias-corrected | grey = no-effect null (+/- 2 SD)",
         x="Cosine similarity (all mice)", y=NULL)
}

################################################################################
### NE DATA (needed by figures 4, 5 and 6)
################################################################################

cat("Loading NE inputs\n")
pb <- readRDS(file.path(NE_BASE, "Buffering_LOO_Split", "pseudobulk_counts_and_meta.RDS"))
pb_all <- pb$pb_all; pb_meta <- pb$pb_meta
age_vals <- unique(pb_meta$age_norm)
AGES <- c(middle = age_vals[!age_vals %in% c("young","old")][1], old = "old")
AGE_LAB <- c(middle="Middle age", old="Old age")
SEX_LAB <- c(female="Females", male="Males")

rob_dir <- file.path(NE_BASE, "Buffering_MiddleAge_Robustness")
c1 <- read.csv(latest(rob_dir, "_Check1_BiasCorrected\\.csv$"), stringsAsFactors=FALSE) %>% dplyr::mutate(age=tolower(age))
c2 <- read.csv(latest(rob_dir, "_Check2_HeldOut_PerMouse\\.csv$"), stringsAsFactors=FALSE) %>% dplyr::mutate(age=tolower(age))

ne_groups <- function(ag, sx, ct) {
  m <- pb_meta[pb_meta$celltype == ct & pb_meta$sex_norm == sx, ]
  list(m=m, y=m$sample[m$age_norm=="young" & m$condition=="HC"],
       o=m$sample[m$age_norm==AGES[[ag]] & m$condition=="HC"],
       e=m$sample[m$age_norm==AGES[[ag]] & m$condition=="EE"])
}
bi_ne <- dplyr::bind_rows(lapply(seq_len(nrow(c1)), function(i) {
  g <- ne_groups(c1$age[i], c1$sex[i], c1$celltype[i])
  if (min(length(g$y), length(g$o), length(g$e)) < 2) return(NULL)
  s <- stats3(log_norm(pb_all[, c(g$y, g$o, g$e), drop=FALSE]), g$y, g$o, g$e)
  data.frame(age=c1$age[i], sex=c1$sex[i], celltype=c1$celltype[i], BI_raw=s[["BI_raw"]], BI_corr=s[["BI_corr"]])
}))
ho_ne <- c2 %>% dplyr::group_by(age, sex, celltype) %>%
  dplyr::summarise(cos_heldout=mean(cosine, na.rm=TRUE), .groups="drop")
ne_ct <- c1 %>% dplyr::select(age, sex, celltype, n_genes, cos_raw, cos_dc, null_mean, null_sd, beyond_null) %>%
  dplyr::left_join(ho_ne, by=c("age","sex","celltype")) %>%
  dplyr::left_join(bi_ne, by=c("age","sex","celltype")) %>%
  dplyr::filter(n_genes >= MIN_GENES) %>%
  dplyr::mutate(pass = !is.na(cos_heldout) & cos_heldout < 0 & beyond_null %in% TRUE,
                shape_lab = SEX_LAB[sex])
calls <- readRDS(file.path(NE_MASTER, "Step2_GeneCalls.RDS"))
ne_act <- read.csv(file.path(NE_MASTER, "Step5_TFActivity_PerCellType.csv"), stringsAsFactors=FALSE) %>%
  dplyr::filter(sex == "female", age == "old") %>%
  dplyr::transmute(celltype, source, NE_S=S, NE_call=dplyr::case_when(buf ~ "buffered", amp ~ "amplified", TRUE ~ "-"))

################################################################################
### FIGURE 4: NE
################################################################################

if (RUN[["fig4"]]) {
  cat("\n== FIGURE 4 ==\n")
  wsrc(ne_ct, "Fig4_NE_CellType_Metrics")
  for (ag in c("middle","old")) {
    sm  <- ne_ct %>% dplyr::filter(age == ag) %>%
      dplyr::mutate(panel=factor(SEX_LAB[sex], levels=SEX_LAB), celltype=factor(celltype, levels=ct_levels(celltype)))
    pts <- c2 %>% dplyr::filter(age == ag) %>% dplyr::semi_join(sm, by=c("sex","celltype")) %>%
      dplyr::mutate(panel=factor(SEX_LAB[sex], levels=SEX_LAB), celltype=factor(celltype, levels=levels(sm$celltype)))
    h <- 25 + 5 * nlevels(sm$celltype)
    sv(plot_cos(sm, pts, paste("Direction of buffering:", AGE_LAB[[ag]]),
                subtitle="Large symbol = mean of held-out fits; grey = each held-out HC mouse"),
       paste0("Fig4a_HeldOutCosine_", ag), 110, h)
    sv(plot_bi(sm, paste("Buffering index:", AGE_LAB[[ag]])), paste0("Fig4b_BufferingIndex_", ag), 110, h)
    sv(plot_dc(sm, paste("Bias-corrected cosine:", AGE_LAB[[ag]])), paste0("ED4_BiasCorrectedCosine_", ag), 110, h)
  }
  q <- ne_ct %>% dplyr::mutate(age_lab=factor(AGE_LAB[age], levels=AGE_LAB), celltype=factor(celltype, levels=CT_ORDER))
  sv(plot_quad(q, "Cosine vs buffering index (NE)", facet="age_lab"), "Fig4c_Quadrant_NE", 183, 90)
  for (sx in c("female","male"))
    sv(plot_quad(q %>% dplyr::filter(sex == sx), paste("Cosine vs buffering index (NE),", SEX_LAB[[sx]]), facet="age_lab"),
       paste0("Fig4c_Quadrant_NE_", sx), 183, 90)

  # 4d heatmaps of buffered genes, gated cell types
  units <- ne_ct %>% dplyr::filter(pass) %>% dplyr::select(age, sex, celltype) %>%
    dplyr::semi_join(calls %>% dplyr::distinct(age, sex, celltype), by=c("age","sex","celltype")) %>%
    dplyr::arrange(sex, age, factor(celltype, levels=CT_ORDER))
  heat_dir <- file.path(OUT, "Fig4d_Heatmaps"); dir.create(heat_dir, showWarnings=FALSE)
  heat_src <- list(); pages <- list()
  for (i in seq_len(nrow(units))) {
    ag <- units$age[i]; sx <- units$sex[i]; ct <- units$celltype[i]
    cl <- calls %>% dplyr::filter(age == ag, sex == sx, celltype == ct)
    tier <- HEAT_TIER; sel <- cl[cl[[tier]] %in% TRUE, ]
    if (nrow(sel) < HEAT_MIN) { tier <- "T3"; sel <- cl[cl$T3 %in% TRUE, ] }
    if (nrow(sel) < 3) next
    sel <- sel %>% dplyr::mutate(dir=ifelse(aging_lfc > 0, "Up with age", "Down with age")) %>%
      dplyr::group_by(dir) %>% dplyr::slice_max(abs(aging_lfc), n=HEAT_PER_DIR, with_ties=FALSE) %>% dplyr::ungroup()
    g <- ne_groups(ag, sx, ct)
    X <- log_norm(pb_all[, c(g$y, g$o, g$e), drop=FALSE])
    genes <- intersect(sel$gene, rownames(X)); if (length(genes) < 3) next
    Z <- t(scale(t(X[genes, , drop=FALSE]))); Z[is.na(Z)] <- 0
    grp <- setNames(c(rep("Young HC", length(g$y)), rep(paste(tools::toTitleCase(ag), "HC"), length(g$o)),
                      rep(paste(tools::toTitleCase(ag), "NE"), length(g$e))), c(g$y, g$o, g$e))
    ord <- sel %>% dplyr::filter(gene %in% genes) %>% dplyr::arrange(dir, aging_lfc)
    long <- as.data.frame(as.table(Z)) %>% setNames(c("gene","sample","z")) %>%
      dplyr::mutate(gene=as.character(gene), sample=as.character(sample)) %>%
      dplyr::left_join(ord %>% dplyr::select(gene, dir, aging_lfc, ne_lfc), by="gene") %>%
      dplyr::mutate(group=factor(grp[sample], levels=unique(grp)),
                    gene=factor(gene, levels=ord$gene), dir=factor(dir, levels=c("Up with age","Down with age")),
                    z=pmax(pmin(z, 2), -2))
    heat_src[[length(heat_src)+1]] <- long %>% dplyr::mutate(age=ag, sex=sx, celltype=ct, tier=tier)
    p <- ggplot(long, aes(sample, gene, fill=z)) + geom_tile(color="white", linewidth=0.1) +
      facet_grid(dir ~ group, scales="free", space="free") +
      scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", limits=c(-2, 2), name="z-score") +
      th + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=4.5),
                 axis.line=element_blank(), panel.spacing=unit(0.6, "mm")) +
      labs(title=sprintf("%s, %s, %s", ct, SEX_LAB[[sx]], AGE_LAB[[ag]]),
           subtitle=sprintf("NE-buffered genes (%s; top %d per direction by aging LFC); each column = one mouse",
                            tier, HEAT_PER_DIR), x=NULL, y=NULL)
    hh <- 30 + 2.4 * length(genes)
    nm <- paste0("Fig4d_", ag, "_", sx, "_", gsub("[^A-Za-z0-9]+", "_", ct))
    ggsave(file.path(heat_dir, paste0(nm, ".pdf")), p, width=mm(70), height=mm(hh), limitsize=FALSE)
    ggsave(file.path(heat_dir, paste0(nm, ".png")), p, width=mm(70), height=mm(hh), dpi=600, bg="white", limitsize=FALSE)
    pages[[length(pages)+1]] <- p
  }
  if (length(pages)) {
    pdf(file.path(OUT, "Fig4d_Heatmaps_AllGatedCellTypes.pdf"), width=mm(70), height=mm(180))
    for (p in pages) print(p); dev.off()
    wsrc(dplyr::bind_rows(heat_src), "Fig4d_Heatmaps")
    cat(sprintf("  heatmaps: %d gated cell types\n", length(pages)))
  }
}

################################################################################
### DREADD DATA (figures 5 and 6)
################################################################################

dr_ct <- NULL; dr_pts <- NULL; dr_band <- NULL
if (RUN[["fig5"]]) {
  cat("\nLoading DREADD inputs\n")
  ho_dir  <- file.path(DR_BASE, "DREADD_Buffering_HeldOut")
  perfit  <- read.csv(latest(ho_dir, "_DREADD_HeldOut_PerFit\\.csv$"), stringsAsFactors=FALSE)
  hsum    <- read.csv(latest(ho_dir, "_DREADD_HeldOut_Summary\\.csv$"), stringsAsFactors=FALSE)
  rows <- list()
  for (ap in APPROACHES) {
    PB <- readRDS(file.path(ho_dir, paste0("pseudobulk_mScarlet_50_Approach_", ap, ".RDS")))
    for (ct in intersect(CT_ORDER, unique(PB$meta$celltype))) {
      d <- PB$meta[PB$meta$celltype == ct, ]
      y <- d$sample[d$age=="young" & d$condition=="PBS"]; o <- d$sample[d$age=="old" & d$condition=="PBS"]
      e <- d$sample[d$age=="old" & d$condition=="C21"]
      if (min(length(y), length(o), length(e)) < 2) next
      X <- log_norm(PB$counts[, c(y, o, e), drop=FALSE]); if (nrow(X) < MIN_GENES) next
      rows[[length(rows)+1]] <- data.frame(approach=ap, celltype=ct, n_genes=nrow(X),
                                           t(stats3(X, y, o, e)), t(sim_null(X, y, o, e)))
      cat(sprintf("  DREADD %s %s\n", ap, ct))
    }
  }
  st <- paste0("mScarlet_50_Approach_", APPROACHES)
  dr_ct <- dplyr::bind_rows(rows) %>%
    dplyr::left_join(hsum %>% dplyr::filter(state %in% paste0("mScarlet_50_Approach_", APPROACHES)) %>%
                       dplyr::transmute(approach, celltype, cos_heldout=cosine_heldout, mouse_null_maxabs), by=c("approach","celltype")) %>%
    dplyr::mutate(beyond_null = cos_dc < null_mean - 2*null_sd,
                  pass = !is.na(cos_heldout) & cos_heldout < 0 & beyond_null, shape_lab="Females")
  dr_pts  <- perfit %>% dplyr::filter(state %in% st, !is.na(cosine)) %>%
    dplyr::mutate(approach=sub(".*Approach_", "", state))
  dr_band <- dr_ct %>% dplyr::filter(is.finite(mouse_null_maxabs)) %>%
    dplyr::transmute(approach, celltype, lo=-mouse_null_maxabs, hi=mouse_null_maxabs)
}

################################################################################
### FIGURE 5: DREADD (B and D separately; NE old female alongside)
################################################################################

if (RUN[["fig5"]]) {
  cat("\n== FIGURE 5 ==\n")
  wsrc(dr_ct, "Fig5_DREADD_CellType_Metrics")
  ne_ref <- ne_ct %>% dplyr::filter(age == "old", sex == "female")
  ne_pts <- c2 %>% dplyr::filter(age == "old", sex == "female") %>% dplyr::semi_join(ne_ref, by="celltype")
  for (ap in APPROACHES) {
    pn_dr <- paste0("DREADD (Approach ", ap, ")"); pn_ne <- "NE (old females)"
    d  <- dr_ct %>% dplyr::filter(approach == ap)
    cts <- ct_levels(c(d$celltype, ne_ref$celltype))
    sm <- dplyr::bind_rows(ne_ref %>% dplyr::mutate(panel=pn_ne), d %>% dplyr::mutate(panel=pn_dr)) %>%
      dplyr::mutate(panel=factor(panel, levels=c(pn_ne, pn_dr)), celltype=factor(celltype, levels=cts))
    pts <- dplyr::bind_rows(ne_pts %>% dplyr::transmute(celltype, cosine, panel=pn_ne),
                            dr_pts %>% dplyr::filter(approach == ap) %>% dplyr::transmute(celltype, cosine, panel=pn_dr)) %>%
      dplyr::filter(celltype %in% cts) %>%
      dplyr::mutate(panel=factor(panel, levels=levels(sm$panel)), celltype=factor(celltype, levels=cts))
    bd <- dr_band %>% dplyr::filter(approach == ap, celltype %in% cts) %>%
      dplyr::mutate(panel=factor(pn_dr, levels=levels(sm$panel)), celltype=factor(celltype, levels=cts))
    h <- 25 + 5 * length(cts)
    sv(plot_cos(sm, pts, "Direction of buffering: NE vs DREADD (old females)", band=bd,
                subtitle="DREADD grey band = C21-vs-C21 mouse-noise range (not available for Approach B neurons)"),
       paste0("Fig5a_HeldOutCosine_Approach", ap), 110, h)
    sv(plot_bi(sm, "Buffering index: NE vs DREADD (old females)"), paste0("Fig5b_BufferingIndex_Approach", ap), 110, h)
    sv(plot_dc(sm, "Bias-corrected cosine: NE vs DREADD"), paste0("ED5_BiasCorrectedCosine_Approach", ap), 110, h)
    sv(plot_quad(d %>% dplyr::mutate(celltype=factor(celltype, levels=CT_ORDER)),
                 paste0("Cosine vs buffering index: DREADD (Approach ", ap, ")"), ref=ne_ref) +
         labs(subtitle="Grey triangles = NE, old females"),
       paste0("Fig5c_Quadrant_Approach", ap), 89, 85)
  }
}

################################################################################
### FIGURE 6: NE-buffered genes, NE and DREADD move together
################################################################################

if (RUN[["fig6"]]) {
  cat("\n== FIGURE 6 ==\n")
  ne_g <- calls %>% dplyr::filter(sex == "female", age == "old") %>%
    dplyr::transmute(celltype, gene, ne_aging_lfc=aging_lfc, ne_lfc, ne_agechg = n_aging >= 2, ne_T3=T3, ne_A3=A3)
  for (ap in APPROACHES) {
    adir <- file.path(DR_VS_NE, paste0("Approach_", ap))
    dr_g <- read.csv(file.path(adir, "Part2_GeneCalls_RealLabeling.csv"), stringsAsFactors=FALSE) %>%
      dplyr::select(celltype, gene, c21_lfc)
    j <- ne_g %>% dplyr::inner_join(dr_g, by=c("celltype","gene"))

    # 6a: NE-buffered genes, aging | NE | DREADD
    a <- j %>% dplyr::filter(celltype %in% FIG6_CTS, ne_T3) %>% dplyr::group_by(celltype) %>%
      dplyr::slice_max(abs(ne_lfc), n=FIG6_N_GENES, with_ties=FALSE) %>% dplyr::ungroup() %>%
      dplyr::arrange(celltype, sign(ne_aging_lfc), ne_aging_lfc) %>%
      dplyr::mutate(row=paste(celltype, gene, sep="::"))
    if (nrow(a)) {
      la <- a %>% tidyr::pivot_longer(c(ne_aging_lfc, ne_lfc, c21_lfc), names_to="col", values_to="lfc") %>%
        dplyr::mutate(col=factor(dplyr::recode(col, ne_aging_lfc="Aging", ne_lfc="NE", c21_lfc="DREADD"),
                                 levels=c("Aging","NE","DREADD")),
                      row=factor(row, levels=rev(unique(a$row))), celltype=factor(celltype, levels=FIG6_CTS))
      p <- ggplot(la, aes(col, row, fill=pmax(pmin(lfc, 1), -1))) + geom_tile(color="white", linewidth=0.2) +
        scale_y_discrete(labels=function(x) sub(".*::", "", x)) +
        facet_grid(celltype ~ ., scales="free_y", space="free_y") +
        scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", limits=c(-1, 1), name="log2 FC") +
        th + theme(axis.text.y=element_text(size=4.5), axis.line=element_blank(), strip.text.y=element_text(angle=0)) +
        labs(title="NE-buffered genes (old females)",
             subtitle=paste0("Aging = old vs young HC; NE = old NE vs old HC; DREADD = old C21 vs old PBS (Approach ", ap, ")"),
             x=NULL, y=NULL)
      sv(p, paste0("Fig6a_NEbufferedGenes_Approach", ap), 75, 25 + 2.6 * nrow(a))
      wsrc(a, paste0("Fig6a_Approach", ap))
    }

    # 6b: gene-level agreement, oriented to NE's aging direction
    b <- j %>% dplyr::filter(celltype %in% FIG6_CTS, ne_agechg) %>%
      dplyr::mutate(s=sign(ne_aging_lfc), x=ne_lfc * s, y=c21_lfc * s,
                    cls=factor(dplyr::case_when(ne_T3 ~ "NE-buffered", ne_A3 ~ "NE-amplified", TRUE ~ "Other aging genes"),
                               levels=c("Other aging genes","NE-amplified","NE-buffered")),
                    celltype=factor(celltype, levels=FIG6_CTS)) %>% dplyr::arrange(cls)
    if (nrow(b)) {
      ann <- b %>% dplyr::group_by(celltype) %>%
        dplyr::summarise(rho=suppressWarnings(cor(x, y, method="spearman")),
                         pct_buf=100 * mean(y[cls=="NE-buffered"] < 0), pct_amp=100 * mean(y[cls=="NE-amplified"] > 0),
                         .groups="drop") %>%
        dplyr::mutate(lab=sprintf("rho = %.2f\nbuffered reversed: %.0f%%\namplified pushed: %.0f%%", rho, pct_buf, pct_amp))
      lx <- quantile(abs(b$x), 0.99, na.rm=TRUE); ly <- quantile(abs(b$y), 0.99, na.rm=TRUE)
      p <- ggplot(b, aes(x, y, color=cls)) +
        geom_hline(yintercept=0, color="grey50", linewidth=0.3) + geom_vline(xintercept=0, color="grey50", linewidth=0.3) +
        geom_point(size=0.35, alpha=0.6) +
        geom_text(data=ann, aes(x=-Inf, y=Inf, label=lab), inherit.aes=FALSE, hjust=-0.05, vjust=1.1, size=1.8) +
        scale_color_manual(values=c(`Other aging genes`=COL_OTH, `NE-amplified`=COL_AMP, `NE-buffered`=COL_BUF), name=NULL) +
        coord_cartesian(xlim=c(-lx, lx), ylim=c(-ly, ly)) + facet_wrap(~celltype, nrow=1) +
        th + theme(legend.position="bottom") + guides(color=guide_legend(override.aes=list(size=1.5, alpha=1))) +
        labs(title=paste0("NE and DREADD move aging genes together (Approach ", ap, ")"),
             subtitle="Both effects oriented to NE's aging direction: < 0 = reverses aging, > 0 = pushes further",
             x="NE effect (log2 FC, aging-oriented)", y="DREADD effect (log2 FC, aging-oriented)")
      sv(p, paste0("Fig6b_GeneAgreement_Approach", ap), 183, 70)
      wsrc(ann, paste0("Fig6b_Stats_Approach", ap))
    }
  }

  # 6c: Part 3 odds ratios, B and D
  p3 <- dplyr::bind_rows(lapply(APPROACHES, function(ap)
    read.csv(file.path(DR_VS_NE, paste0("Approach_", ap), "Part3_Overlap_withNE_PerCellType.csv"), stringsAsFactors=FALSE) %>%
      dplyr::mutate(approach=ap))) %>%
    dplyr::filter(n_ne_set >= 20, !grepl("\\(", test)) %>%
    dplyr::mutate(test=ifelse(grepl("buffered", test), "NE-buffered genes reversed by DREADD",
                              "NE-amplified genes pushed further by DREADD"),
                  test=factor(test, levels=c("NE-buffered genes reversed by DREADD","NE-amplified genes pushed further by DREADD")),
                  celltype=factor(celltype, levels=ct_levels(celltype)),
                  lab=sprintf("%d/%d", rank, n_labelings))
  wsrc(p3, "Fig6c_Overlap")
  p <- ggplot(p3, aes(odds_ratio, celltype, color=approach)) +
    geom_vline(xintercept=1, linetype="dashed", color="grey40", linewidth=0.3) +
    geom_point(position=position_dodge(width=0.6), size=1.8) +
    geom_text(aes(label=lab), position=position_dodge(width=0.6), hjust=-0.35, size=1.7, show.legend=FALSE) +
    scale_x_log10() + scale_color_manual(values=COL_AP, name="Approach") +
    facet_wrap(~test) + th + theme(legend.position="bottom") +
    labs(title="DREADD moves NE's genes the way NE does",
         subtitle="Odds ratio vs other NE aging genes; label = rank of the real labeling among relabelings",
         x="Odds ratio (log scale)", y=NULL)
  sv(p, "Fig6c_Overlap_OddsRatio", 150, 25 + 5 * nlevels(p3$celltype))

  # 6d: TF direction agreement, real vs relabelings
  ne_called <- ne_act %>% dplyr::filter(NE_call != "-")
  ag <- dplyr::bind_rows(lapply(APPROACHES, function(ap) {
    U <- lapply(list.files(file.path(DR_VS_NE, paste0("Approach_", ap), "units"), pattern="\\.RDS$", full.names=TRUE), readRDS)
    dplyr::bind_rows(lapply(U, function(u) if (!is.null(u$tf) && nrow(u$tf)) dplyr::mutate(u$tf, celltype=u$celltype))) %>%
      dplyr::inner_join(ne_called, by=c("celltype","source")) %>%
      dplyr::group_by(celltype, key, real) %>%
      dplyr::summarise(pct=100 * mean(sign(S_NEax) == sign(NE_S), na.rm=TRUE), n_tf=dplyr::n(), .groups="drop") %>%
      dplyr::mutate(approach=ap)
  })) %>% dplyr::filter(n_tf >= 20) %>% dplyr::mutate(celltype=factor(celltype, levels=ct_levels(celltype)))
  wsrc(ag, "Fig6d_TFAgreement")
  p <- ggplot(ag, aes(pct, celltype)) +
    geom_vline(xintercept=50, linetype="dashed", color="grey40", linewidth=0.3) +
    geom_point(data=ag %>% dplyr::filter(!real), color="grey70", size=0.7,
               position=position_jitter(height=0.15, width=0, seed=1)) +
    geom_point(data=ag %>% dplyr::filter(real), aes(fill=celltype), shape=23, size=2.4, stroke=0.3) +
    scale_fill_manual(values=CT_COLS, guide="none") + facet_wrap(~paste("Approach", approach)) + th +
    labs(title="TF activity: DREADD agrees with NE",
         subtitle="% of NE-called TFs with the same sign in DREADD; diamond = real, grey = relabeled mice",
         x="% same direction as NE", y=NULL)
  sv(p, "Fig6d_TFAgreement", 120, 25 + 5 * nlevels(ag$celltype))

  # 6e / 6f: TF scores NE vs DREADD
  tf <- dplyr::bind_rows(lapply(APPROACHES, function(ap)
    read.csv(file.path(DR_VS_NE, paste0("Approach_", ap), "Part4_TFActivity_PerCellType_withNE.csv"), stringsAsFactors=FALSE) %>%
      dplyr::mutate(approach=ap))) %>% dplyr::filter(!is.na(NE_S), NE_call != "-")
  wsrc(tf %>% dplyr::select(approach, celltype, source, NE_call, NE_S, S_NEax, rank_buf_NEax, rank_amp_NEax, DR_call_NEax),
       "Fig6ef_TFScores")
  hl <- tf %>% dplyr::filter(source %in% LABEL_TFS) %>% dplyr::mutate(lab=paste(source, short_ct(celltype)))
  pe <- ggplot(tf, aes(NE_S, S_NEax)) +
    geom_hline(yintercept=0, color="grey50", linewidth=0.3) + geom_vline(xintercept=0, color="grey50", linewidth=0.3) +
    geom_point(aes(color=NE_call), size=0.5, alpha=0.45) +
    geom_point(data=hl, aes(fill=source), shape=21, size=1.8, stroke=0.25) +
    scale_color_manual(values=c(buffered=COL_BUF, amplified=COL_AMP), name="NE call") +
    scale_fill_brewer(palette="Set1", name=NULL) +
    facet_wrap(~paste("Approach", approach), scales="free_y") + th + theme(legend.position="bottom") +
    labs(title="TF buffering scores: NE vs DREADD (old females)",
         subtitle="Each point = TF x cell type called by NE; both scored against NE's aging direction (> 0 = reverses aging)",
         x="NE buffering score", y="DREADD buffering score")
  pe <- pe + if (HAS_REPEL) ggrepel::geom_text_repel(data=hl %>% dplyr::filter(source %in% c("Nr3c1","Creb1")),
                                                     aes(label=lab), size=1.6, segment.size=0.2, max.overlaps=40)
             else geom_text(data=hl %>% dplyr::filter(source %in% c("Nr3c1","Creb1")), aes(label=lab), size=1.6, vjust=-0.8)
  sv(pe, "Fig6e_TFScores_NEvsDREADD", 183, 85)

  for (tfn in FIG6_TFS) {
    d <- tf %>% dplyr::filter(source == tfn)
    if (!nrow(d)) next
    long <- dplyr::bind_rows(d %>% dplyr::distinct(celltype, NE_S) %>% dplyr::transmute(celltype, who="NE", S=NE_S, z=NE_S / max(abs(NE_S))),
                             d %>% dplyr::group_by(approach) %>%
                               dplyr::transmute(celltype, who=paste("DREADD", approach), S=S_NEax, z=S_NEax / max(abs(S_NEax))) %>%
                               dplyr::ungroup() %>% dplyr::select(-approach)) %>%
      dplyr::mutate(celltype=factor(celltype, levels=ct_levels(celltype)),
                    who=factor(who, levels=c("NE", paste("DREADD", APPROACHES))))
    p <- ggplot(long, aes(z, celltype, fill=who)) +
      geom_col(position=position_dodge(width=0.75), width=0.7) +
      geom_vline(xintercept=0, color="grey30", linewidth=0.3) +
      scale_fill_manual(values=c(NE="grey30", setNames(COL_AP, paste("DREADD", names(COL_AP)))), name=NULL) +
      scale_x_continuous(limits=c(-1, 1)) + th + theme(legend.position="bottom") +
      labs(title=paste0(tfn, ": NE vs DREADD"),
           subtitle="Buffering score scaled to each dataset's maximum; < 0 = amplified, > 0 = buffered",
           x="Scaled TF buffering score", y=NULL)
    sv(p, paste0("Fig6f_", tfn, "_NEvsDREADD"), 89, 25 + 6 * nlevels(long$celltype))
  }
}

cat("\n=== DONE ===\nOutputs in:", OUT, "\n")
