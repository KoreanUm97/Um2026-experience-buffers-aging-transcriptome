################################################################################
### Fig6_ED11_TF_Figures.R
###
### Panels: Fig 6a = bubble_tp_fdr/both_<sex>_old; Fig 6b = dreadd_recap/score_OwnAxis_buffered_D;
###         ED11a = dreadd_recap/score_OwnAxis_buffered_B; ED11b = dreadd_recap/score_NEaxis_{D,B}
###
### Figure 6 TF figures built only from the original quantities:
###   aging t       = ULM activity t-value, old HC vs young HC      (mean of 3 held-out fits)
###   stimulation t = ULM activity t-value, old eNE vs old HC       (mean of 3 held-out fits)
###   buffering score = -(aging t x stimulation t), mean of 3 fits  (> 0 buffered, < 0 amplified)
###   buffered / amplified calls = the Figure 6 calls (buf / amp)
### Pairs (one TF in one cell type) are ranked by their buffering score (its size for amplified TFs).
### TFs (across cell types) are ranked by RANK_BY; default "sum" = buffering score summed over the cell types where called.
###
### Input : Buffering_Candidates_Master/Step5_TFActivity_PerCellType.csv (Buffering_Pipeline.R step 2; nothing is recomputed)
###         TF_Significance/p/TFSig_PerFit.csv (TF_Buffering_Significance.R; stars on the 6a tiles)
###         Buffering_DREADD_vs_NE/Approach_{B,D}/units/*.RDS (Buffering_DREADD_vs_NE.R; dreadd_recap)
### Output: Buffering_Candidates_Master/Figure6_TF_Score/<design>/
###
### Designs (choose with RUN). "pair" = one TF in one cell type, e.g. "Creb1 (DG)".
###   bubble_rank      : one dot per TF, ranked by its top buffering score (x); dot size = number of cell types in which
###                      it is called. Buffered only, amplified only, and both together (buffered right, amplified left).
###   bubble_t         : bubble_rank plus tiles with the aging t and stimulation t of the cell type that gives the top score
###                      (the two values the plotted score is built from), and that cell type. Buffered, amplified, both.
###   bubble_tp        : bubble_t with the ULM significance of each activity shown as stars on the tiles
###                      (median p, or FDR, over the 3 held-out fits; * < 0.05, ** < 0.01, *** < 0.001).
###                      Descriptive only: the calls and ranking are unchanged. Needs PERFIT_CSV.
###   split_heat       : top TFs x cell types, only called cells filled, number = score. Buffered and amplified separately.
###   bar_t            : top buffered (right) and amplified (left) pairs as bars, with aging t and stimulation t tiles.
###   lollipop_t       : top pairs as a lollipop, with aging t and stimulation t tiles.
###   circles_t        : top pairs, score printed in a circle, with aging t and stimulation t tiles.
###   percell_table    : each cell type side by side; top TFs with score, aging t, stimulation t.
###   percell_lollipop : each cell type side by side; lollipop + aging t and stimulation t tiles.
###   top_score_heat   : top TFs x all cell types, buffering score in every cell (blue buffered, red amplified); box = called.
###   top_t_heat       : top TFs x all cell types, aging (A) and stimulation (S) t-values in every cell; box = called.
###   dreadd_recap     : old females only. NE stimulation t vs DREADD stimulation t for every TF x cell type, neurons and glia
###                      side by side; top NE-buffered / -amplified TFs per cell type coloured; bars = % changed in the same
###                      direction as eNE. Needs the Buffering_DREADD_vs_NE units cache (DR_BASE).
### Each design is drawn per sex x age. split_heat, bar_t and lollipop_t are also drawn as one figure with all four
### groups (sexes side by side, middle above old), using the same scales in every panel.
### Existing figures are skipped; delete a figure to redraw it.
################################################################################

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
if (!requireNamespace("patchwork", quietly=TRUE)) stop("Please install patchwork: install.packages('patchwork')")
library(patchwork)

################################################################################
### SETTINGS
################################################################################

base_path <- "PATH_TO_NE_OUTPUT"
MASTER <- file.path(base_path, "Buffering_Candidates_Master")
OUT    <- file.path(MASTER, "Figure6_TF_Score")

RUN <- c("bubble_t", "bubble_tp", "bubble_rank", "split_heat", "bar_t", "lollipop_t", "circles_t", "percell_table", "percell_lollipop",
         "top_score_heat", "top_t_heat", "dreadd_recap")
N_TOP     <- 15     # rows per figure (TFs or pairs) for single-group figures
N_BAR     <- 10     # bar_t: buffered and amplified pairs each, single-group figure
N_GROUPS  <- 15     # rows per panel (per direction) in the all-groups figures
N_PER_CT  <- 5      # percell_*: TFs per cell type
N_BUBBLE  <- 20     # bubble_rank: TFs per direction (single-group figures)
PERFIT_CSV <- file.path(MASTER, "TF_Significance", "p", "TFSig_PerFit.csv")   # per-fit ULM p-values (from TF_Buffering_Significance.R)
P_TYPE     <- "fdr"   # stars on the bubble_tp tiles: "fdr" = BH FDR (within cell type and contrast), "p" = raw ULM p-value
                      # output goes to bubble_tp_fdr/ or bubble_tp_p/, so both versions can be kept side by side
RANK_BY   <- "sum"   # how TFs are ranked in the TF-level designs (bubble_rank, bubble_t, split_heat, top_score_heat, top_t_heat):
                      #   "sum"   = buffering score summed over the cell types where the TF is called (strength + recurrence)
                      #   "top"   = highest buffering score in any one cell type
                      #   "count" = number of cell types called (ties broken by summed score)
N_BUBBLE_GROUPS <- 15   # bubble_rank / bubble_t: TFs per direction in the all-groups figures (also per side in "both")

# dreadd_recap: do NE-buffered / NE-amplified TFs (old females) change the same way with DREADD activation?
DR_BASE  <- "PATH_TO_DREADD_OUTPUT"
DR_APPROACHES <- c("B", "D")
DR_TOP_PER_CT <- 5               # top NE-called TFs per cell type and direction (by buffering score) that are coloured (ED11b: top 5)
DR_N_LABEL    <- 5               # labelled TFs per direction per panel
DR_SCORE_CLIP <- 60              # buffering-score scatters: axis values clipped at +/- this (points beyond sit on the edge)
DR_N_LABEL_BUF <- 12             # score_OwnAxis_buffered: labelled buffered TFs per panel (neurons / glia)
DR_TOP_PER_CT_BUF <- 10          # score_OwnAxis_buffered: top eNE-buffered TFs per cell type that are coloured
DR_NEURONAL   <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN")

CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN","Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell")
CT_SHORT <- c("DG ExNs"="DG","CA3 ExNs"="CA3","CA1 ExNs"="CA1","CA2 ExNs"="CA2","Sub ExNs"="Sub","IN"="IN",
              "Astrocytes"="Astro","Oligodendrocytes"="Oligo","OPC"="OPC","Microglia"="Micro","Endothelial"="Endo",
              "Pericyte"="Peri","APC"="APC","Stromal"="Stromal","Cajal-Retzius Cell"="CR")
GROUP_ORDER <- c("female middle", "male middle", "female old", "male old")   # 2 x 2: sexes side by side, middle on top
COL_BUF <- "#2B6CA3"; COL_AMP <- "#B2402A"; T_LOW <- "#2166AC"; T_HIGH <- "#B2182B"

################################################################################
### DATA
################################################################################

f5 <- file.path(MASTER, "Step5_TFActivity_PerCellType.csv")
if (!file.exists(f5)) stop("Not found: ", f5)
raw <- read.csv(f5, stringsAsFactors=FALSE)
ct_lev <- c(intersect(CT_ORDER, unique(raw$celltype)), setdiff(unique(raw$celltype), CT_ORDER))
short_of <- function(x) ifelse(x %in% names(CT_SHORT), CT_SHORT[x], x)
d <- raw %>% dplyr::mutate(
  ct    = factor(short_of(celltype), levels=unname(short_of(ct_lev))),
  pair  = paste0(source, " (", ct, ")"),
  score = abs(S),
  dir   = dplyr::case_when(buf %in% TRUE ~ "Buffered", amp %in% TRUE ~ "Amplified", TRUE ~ NA_character_),
  group = paste(sex, age))
# ULM p-values for each activity: median over the 3 held-out fits (real labels)
HAVE_P <- file.exists(PERFIT_CSV)
if (HAVE_P) {
  pf <- read.csv(PERFIT_CSV, stringsAsFactors=FALSE)
  pcolB <- if (P_TYPE == "fdr") "fdr_B" else "p_B"; pcolE <- if (P_TYPE == "fdr") "fdr_E" else "p_E"
  pm <- pf %>% dplyr::group_by(age, sex, celltype, source) %>%
    dplyr::summarise(p_aging = stats::median(.data[[pcolB]], na.rm=TRUE), p_stim = stats::median(.data[[pcolE]], na.rm=TRUE),
                     .groups="drop")
  d <- d %>% dplyr::left_join(pm, by=c("age", "sex", "celltype", "source"))
  cat(sprintf("p-values joined for %.0f%% of rows (%s)\n", 100 * mean(!is.na(d$p_aging)), P_TYPE))
} else {
  d <- d %>% dplyr::mutate(p_aging = NA_real_, p_stim = NA_real_)
  if ("bubble_tp" %in% RUN) cat("NOTE: ", PERFIT_CSV, " not found; bubble_tp is skipped and p columns are empty\n")
}
P_CAPTION <- paste0("Stars: ULM ", if (P_TYPE == "fdr") "FDR" else "p-value", " of each activity, median of 3 held-out fits; * < 0.05, ** < 0.01, *** < 0.001")
stars <- function(p) dplyr::case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
called <- d %>% dplyr::filter(!is.na(dir))
SLIM <- max(abs(called$S));                                   # shared scales across all figures
TLIM <- max(abs(c(called$aging_act, called$ne_act)))
groups <- intersect(GROUP_ORDER, unique(d$group))
cat(sprintf("Loaded %d rows: %d buffered, %d amplified calls\n", nrow(d), sum(called$dir == "Buffered"), sum(called$dir == "Amplified")))
cat("Making:", paste(RUN, collapse=", "), "\n")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
write.csv(called %>% dplyr::group_by(dir, sex, age) %>% dplyr::arrange(dplyr::desc(score), .by_group=TRUE) %>%
            dplyr::mutate(rank = dplyr::row_number()) %>% dplyr::ungroup() %>%
            dplyr::select(direction=dir, sex, age, rank, source, celltype, buffering_score=S, aging_t=aging_act, stimulation_t=ne_act,
                          aging_p=p_aging, stimulation_p=p_stim),
          file.path(OUT, "TF_ScoreRanking.csv"), row.names=FALSE)

################################################################################
### HELPERS
################################################################################

fmt <- function(x, digits=1) {                                   # digits may be a vector (one per value)
  s <- mapply(function(v, dg) formatC(v, format="f", digits=dg), x, rep_len(digits, length(x)))
  sub("^-(0(\\.0+)?)$", "\\1", unname(s))
}
nice_group <- function(g) { p <- strsplit(g, " ")[[1]]; paste0(tools::toTitleCase(p[1]), ", ", p[2]) }
save_fig <- function(p, design, stem, w, h) {
  dir <- file.path(OUT, design); dir.create(dir, recursive=TRUE, showWarnings=FALSE)
  f <- file.path(dir, stem)
  if (file.exists(paste0(f, ".pdf"))) { cat("  exists, skipped:", design, "/", stem, "\n"); return(invisible()) }
  ggsave(paste0(f, ".pdf"), p, width=w, height=h, limitsize=FALSE)
  ggsave(paste0(f, ".png"), p, width=w, height=h, dpi=300, bg="white", limitsize=FALSE)
  cat("  saved:", design, "/", stem, "\n")
}
th_cl <- function(base=12) theme_classic(base_size=base) + theme(plot.title=element_text(face="bold"))
th_mn <- function(base=12) theme_minimal(base_size=base) + theme(panel.grid=element_blank(), plot.title=element_text(face="bold"))
t_fill <- function() scale_fill_gradient2(low=T_LOW, mid="white", high=T_HIGH, limits=c(-TLIM, TLIM),
                                          name="t-value\nred = up\nblue = down")
dcol <- function(dirn) if (dirn == "Buffered") COL_BUF else COL_AMP
top_pairs <- function(g, dirn, n) called %>% dplyr::filter(group == g, dir == dirn) %>%
  dplyr::slice_max(score, n=n, with_ties=FALSE)
# t-value tiles for a set of pairs; `lev` = row order (bottom to top)
t_tiles <- function(x, lev, base=12) {
  tl <- x %>% dplyr::mutate(pair = factor(pair, lev)) %>% dplyr::select(pair, Aging=aging_act, Stimulation=ne_act) %>%
    tidyr::pivot_longer(-pair, names_to="k", values_to="t")
  ggplot(tl, aes(k, pair, fill=t)) + geom_tile(colour="white", linewidth=0.7) +
    geom_text(aes(label=fmt(t)), size=base / 4.2) + t_fill() +
    scale_x_discrete(position="top", labels=c(Aging="Aging\nt", Stimulation="Stim.\nt")) + labs(x=NULL, y=NULL) +
    th_mn(base) + theme(axis.text.y=element_blank())
}
no_data <- function(msg) ggplot() + annotate("text", x=0, y=0, label=msg, colour="grey50") + theme_void()

################################################################################
### DESIGNS (each returns a plot for one group, or NULL)
################################################################################

# bubble_rank: one dot per TF; x = top buffering score (no averaging), size = number of cell types called
tf_summary <- function(g, dirn) called %>% dplyr::filter(group == g, dir == dirn) %>% dplyr::group_by(source) %>%
  dplyr::summarise(n_ct = dplyr::n(), top_score = max(score), sum_score = sum(score),
                   top_ct = as.character(ct[which.max(score)]),
                   top_aging_t = aging_act[which.max(score)], top_stim_t = ne_act[which.max(score)],
                   top_p_aging = p_aging[which.max(score)], top_p_stim = p_stim[which.max(score)],
                   .groups="drop") %>%
  dplyr::mutate(dir = dirn, rank_value = switch(RANK_BY, sum = sum_score, top = top_score, count = n_ct),
                rank_tie = sum_score) %>%
  dplyr::arrange(dplyr::desc(rank_value), dplyr::desc(rank_tie))
top_tfs <- function(x, n) head(x, n)              # tf_summary output is already sorted by the ranking
RANK_LAB  <- c(sum = "Buffering score, summed over cell types", top = "Top buffering score", count = "Cell types called")[[RANK_BY]]
RANK_HEAD <- c(sum = "Summed score", top = "Top score", count = "Cell types")[[RANK_BY]]
RLIM <- max(called %>% dplyr::group_by(group, dir, source) %>%
              dplyr::summarise(v = switch(RANK_BY, sum = sum(score), top = max(score), count = dplyr::n()), .groups="drop") %>% dplyr::pull(v))
NCT_MAX <- max(called %>% dplyr::count(group, dir, source) %>% dplyr::pull(n))
bubble_size <- function() scale_size_area(max_size=9, limits=c(0, NCT_MAX),
                                          breaks=unique(pmax(1, round(seq(1, NCT_MAX, length.out=4)))), name="Cell types\ncalled")
plot_bubble_rank <- function(g, dirn, n=N_BUBBLE, base=12, title=paste0(dirn, ": ", nice_group(g))) {
  x <- tf_summary(g, dirn); if (!nrow(x)) return(NULL)
  x <- top_tfs(x, n) %>% dplyr::mutate(source = factor(source, rev(source)))
  ggplot(x, aes(rank_value, source)) + geom_segment(aes(x=0, xend=rank_value, yend=source), colour="grey85") +
    geom_point(aes(size=n_ct), colour=dcol(dirn), alpha=0.85) + bubble_size() +
    scale_x_continuous(limits=c(0, RLIM), expand=expansion(mult=c(0, 0.04))) +
    labs(title=title, x=RANK_LAB, y=NULL) + th_cl(base)
}
plot_bubble_both <- function(g, n=ceiling(N_BUBBLE / 2), base=12, title=nice_group(g)) {
  x <- dplyr::bind_rows(top_tfs(tf_summary(g, "Buffered"), n), top_tfs(tf_summary(g, "Amplified"), n))
  if (!nrow(x)) return(NULL)
  x <- x %>% dplyr::mutate(xs = ifelse(dir == "Buffered", rank_value, -rank_value), key = paste(source, dir, sep="__")) %>%
    dplyr::arrange(xs, ifelse(dir == "Buffered", rank_tie, -rank_tie)) %>% dplyr::mutate(key = factor(key, key))
  ggplot(x, aes(xs, key)) + geom_vline(xintercept=0, linewidth=0.3) +
    geom_segment(aes(x=0, xend=xs, yend=key), colour="grey85") +
    geom_point(aes(size=n_ct, colour=dir), alpha=0.85) + bubble_size() +
    scale_colour_manual(values=c(Buffered=COL_BUF, Amplified=COL_AMP), name=NULL) +
    scale_y_discrete(labels=function(v) sub("__.*$", "", v)) + scale_x_continuous(limits=c(-RLIM, RLIM)) +
    labs(title=title, x=paste0(RANK_LAB, "\n(left: amplified   right: buffered)"), y=NULL) + th_cl(base)
}

# bubble_t: bubble ranking + t-value tiles for the cell type giving the top score
plot_bubble_t <- function(g, dirn, n=N_BUBBLE, base=12, title=NULL, collect=TRUE, show_p=FALSE) {
  if (dirn == "Both") {
    k <- ceiling(n / 2)
    x <- dplyr::bind_rows(top_tfs(tf_summary(g, "Buffered"), k), top_tfs(tf_summary(g, "Amplified"), k))
  } else x <- top_tfs(tf_summary(g, dirn), n)
  if (!nrow(x)) return(NULL)
  if (is.null(title)) title <- if (dirn == "Both") nice_group(g) else paste0(dirn, ": ", nice_group(g))
  x <- x %>% dplyr::mutate(xs = ifelse(dir == "Buffered", rank_value, -rank_value), key = paste(source, dir, sep="__"),
                           ord = if (dirn == "Both") xs else rank_value, tie = ifelse(dir == "Buffered", rank_tie, -rank_tie)) %>%
    dplyr::arrange(ord, tie) %>% dplyr::mutate(key = factor(key, key))
  xvar <- if (dirn == "Both") "xs" else "rank_value"
  p1 <- ggplot(x, aes(.data[[xvar]], key)) + geom_segment(aes(x=0, xend=.data[[xvar]], yend=key), colour="grey85") +
    geom_point(aes(size=n_ct, colour=dir), alpha=0.85) + bubble_size() +
    scale_colour_manual(values=c(Buffered=COL_BUF, Amplified=COL_AMP), name=NULL, guide=if (dirn == "Both") "legend" else "none") +
    scale_y_discrete(labels=function(v) sub("__.*$", "", v)) +
    labs(title=title, y=NULL,
         x=if (dirn == "Both") paste0(RANK_LAB, "\n(left: amplified   right: buffered)") else RANK_LAB) +
    th_cl(base)
  p1 <- if (dirn == "Both") p1 + geom_vline(xintercept=0, linewidth=0.3) + scale_x_continuous(limits=c(-RLIM, RLIM)) else
    p1 + scale_x_continuous(limits=c(0, RLIM), expand=expansion(mult=c(0, 0.04)))
  tl <- x %>% dplyr::select(key, top_ct, Aging=top_aging_t, Stimulation=top_stim_t) %>%
    tidyr::pivot_longer(c(Aging, Stimulation), names_to="k", values_to="t") %>%
    dplyr::left_join(x %>% dplyr::select(key, Aging=top_p_aging, Stimulation=top_p_stim) %>%
                       tidyr::pivot_longer(c(Aging, Stimulation), names_to="k", values_to="p"), by=c("key", "k")) %>%
    dplyr::mutate(lab = if (show_p) paste0(fmt(t), stars(p)) else fmt(t))
  p2 <- ggplot(tl, aes(k, key, fill=t)) + geom_tile(colour="white", linewidth=0.7) + geom_text(aes(label=lab), size=base / 4.2) +
    geom_text(data=x, aes(x=2.7, y=key, label=top_ct), inherit.aes=FALSE, hjust=0, size=base / 4.2, colour="grey30") +
    annotate("text", x=2.7, y=nrow(x) + 0.9, label="Cell type", hjust=0, size=base / 4.2, fontface="bold", colour="grey30") +
    t_fill() + scale_x_discrete(position="top", labels=c(Aging="Aging\nt", Stimulation="Stim.\nt"), expand=expansion(add=c(0.6, 1.6))) +
    coord_cartesian(clip="off") + labs(x=NULL, y=NULL) + th_mn(base) + theme(axis.text.y=element_blank())
  (p1 | p2) + plot_layout(widths=c(2.4, 1.35), guides=if (collect) "collect" else "auto")
}

# split_heat: top TFs x cell types, only called cells filled
plot_split_heat <- function(g, dirn, n=N_TOP, base=12, title=paste0(dirn, ": ", nice_group(g))) {
  x <- called %>% dplyr::filter(group == g, dir == dirn); if (!nrow(x)) return(NULL)
  top <- top_tfs(tf_summary(g, dirn), n)
  x <- x %>% dplyr::filter(source %in% top$source) %>% dplyr::mutate(source = factor(source, rev(top$source)))
  cts <- levels(droplevels(d$ct[d$group == g]))
  grid <- expand.grid(source = levels(x$source), ct = cts) %>%
    dplyr::mutate(source = factor(source, levels(x$source)), ct = factor(ct, cts))
  x$ct <- factor(as.character(x$ct), cts)
  ggplot(x, aes(ct, source)) + geom_tile(data=grid, fill="grey96", colour="white", linewidth=0.8) +
    geom_tile(aes(fill=score), colour="white", linewidth=0.8) + geom_text(aes(label=fmt(score, 0)), size=base / 3.8) +
    scale_fill_gradient(low="white", high=dcol(dirn), limits=c(0, SLIM), name="Buffering\nscore") +
    scale_x_discrete(position="top", drop=FALSE) + labs(title=title, x=NULL, y=NULL) +
    th_mn(base) + theme(axis.text.x=element_text(face="bold"))
}

# bar_t: buffered right, amplified left, with t tiles
plot_bar_t <- function(g, n=N_BAR, base=12, title=nice_group(g)) {
  cc <- dplyr::bind_rows(top_pairs(g, "Buffered", n), top_pairs(g, "Amplified", n)); if (!nrow(cc)) return(NULL)
  cc <- cc %>% dplyr::arrange(S); lev <- unique(cc$pair); cc$pair <- factor(cc$pair, lev)
  p1 <- ggplot(cc, aes(S, pair, fill=dir)) + geom_col(width=0.75) + geom_vline(xintercept=0, linewidth=0.3) +
    scale_fill_manual(values=c(Buffered=COL_BUF, Amplified=COL_AMP), guide="none") + scale_x_continuous(limits=c(-SLIM, SLIM)) +
    labs(title=title, x="Buffering score\n(left: amplified   right: buffered)", y=NULL) + th_cl(base)
  (p1 | t_tiles(cc, lev, base)) + plot_layout(widths=c(3, 1.2))
}

# lollipop_t
plot_lollipop_t <- function(g, dirn, n=N_TOP, base=12, title=paste0(dirn, ": ", nice_group(g))) {
  x <- top_pairs(g, dirn, n); if (!nrow(x)) return(NULL)
  x <- x %>% dplyr::arrange(score); lev <- x$pair; x$pair <- factor(x$pair, lev)
  p1 <- ggplot(x, aes(score, pair)) + geom_segment(aes(x=0, xend=score, yend=pair), colour="grey75", linewidth=1) +
    geom_point(colour=dcol(dirn), size=3.6) + scale_x_continuous(limits=c(0, SLIM), expand=expansion(mult=c(0, 0.04))) +
    labs(title=title, x=if (dirn == "Buffered") "Buffering score" else "Buffering score (size)", y=NULL) + th_cl(base)
  (p1 | t_tiles(x, lev, base)) + plot_layout(widths=c(2.6, 1.2))
}

# circles_t
plot_circles_t <- function(g, dirn, n=N_TOP, base=12) {
  x <- top_pairs(g, dirn, n); if (!nrow(x)) return(NULL)
  x <- x %>% dplyr::arrange(score); lev <- x$pair; x$pair <- factor(x$pair, lev)
  p1 <- ggplot(x, aes(1, pair)) + geom_point(aes(size=score), colour=dcol(dirn)) +
    geom_text(aes(label=fmt(score, 0)), colour="white", size=base / 3.8, fontface="bold") +
    scale_size_area(max_size=10, limits=c(0, SLIM), guide="none") +
    scale_x_continuous(breaks=1, labels="Score", position="top") +
    labs(title=paste0(dirn, ": ", nice_group(g)), x=NULL, y=NULL) + th_cl(base) +
    theme(axis.line=element_blank(), axis.ticks=element_blank())
  (p1 | t_tiles(x, lev, base)) + plot_layout(widths=c(1.1, 1.5))
}

# percell_table: each cell type side by side; score, aging t, stimulation t
plot_percell_table <- function(g, dirn, n=N_PER_CT, base=11) {
  x <- called %>% dplyr::filter(group == g, dir == dirn); if (!nrow(x)) return(NULL)
  x <- x %>% dplyr::group_by(ct) %>% dplyr::slice_max(score, n=n, with_ties=FALSE) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>% dplyr::ungroup() %>% dplyr::mutate(ct = droplevels(ct))
  l <- x %>% dplyr::select(ct, rank, source, Score=score, Aging=aging_act, Stim.=ne_act) %>%
    tidyr::pivot_longer(c(Score, Aging, Stim.), names_to="k", values_to="v") %>%
    dplyr::mutate(xk = match(k, c("Score", "Aging", "Stim.")), tv = ifelse(k == "Score", NA, v))
  ggplot(l, aes(xk, -rank)) + geom_tile(aes(fill=tv), colour="white", linewidth=0.8) +
    geom_tile(data=l %>% dplyr::filter(k == "Score"), fill="grey90", colour="white", linewidth=0.8) +
    geom_text(aes(label=fmt(v, ifelse(k == "Score", 0, 1)), fontface=ifelse(k == "Score", "bold", "plain")), size=base / 3.8) +
    geom_text(data=x, aes(x=0.4, y=-rank, label=source), inherit.aes=FALSE, hjust=1, size=base / 3.6) +
    facet_wrap(~ ct, nrow=1) +
    scale_fill_gradient2(low=T_LOW, mid="white", high=T_HIGH, limits=c(-TLIM, TLIM), na.value="grey90", name="t-value") +
    scale_x_continuous(breaks=1:3, labels=c("Score", "Aging", "Stim."), position="top", limits=c(-1.6, 3.5)) +
    coord_cartesian(clip="off") + labs(title=paste0(dirn, ": ", nice_group(g), ", top ", n, " per cell type"), x=NULL, y=NULL) +
    th_mn(base) + theme(axis.text.y=element_blank(), axis.text.x=element_text(size=base - 3),
                        strip.text=element_text(face="bold", size=base + 1), panel.spacing=unit(0.6, "lines"))
}

# percell_lollipop: each cell type side by side; lollipop + t tiles
plot_percell_lollipop <- function(g, dirn, n=N_PER_CT, base=10) {
  x <- called %>% dplyr::filter(group == g, dir == dirn); if (!nrow(x)) return(NULL)
  cts <- levels(droplevels(x$ct))
  one <- function(cn) {
    y <- x %>% dplyr::filter(ct == cn) %>% dplyr::slice_max(score, n=n, with_ties=FALSE) %>% dplyr::arrange(score)
    lev <- y$pair; y$pair <- factor(y$pair, lev)
    a <- ggplot(y, aes(score, pair)) + geom_segment(aes(x=0, xend=score, yend=pair), colour="grey75", linewidth=0.9) +
      geom_point(colour=dcol(dirn), size=2.8) + scale_x_continuous(limits=c(0, SLIM)) +
      scale_y_discrete(labels=function(v) sub(" \\(.*\\)$", "", v)) +
      labs(title=cn, x="Score", y=NULL) + th_cl(base) + theme(plot.title=element_text(hjust=0.5))
    (a | t_tiles(y, lev, base)) + plot_layout(widths=c(1.6, 1))
  }
  wrap_plots(lapply(cts, one), ncol=min(4, length(cts))) + plot_layout(guides="collect") +
    plot_annotation(title=paste0(dirn, ": ", nice_group(g), ", top ", n, " per cell type"),
                    theme=theme(plot.title=element_text(face="bold", size=14)))
}

# top_score_heat: top TFs (by best score in this direction) x all cell types, buffering score everywhere
plot_top_score_heat <- function(g, dirn, n=N_TOP, base=12) {
  x <- called %>% dplyr::filter(group == g, dir == dirn); if (!nrow(x)) return(NULL)
  best <- top_tfs(tf_summary(g, dirn), n)
  lev <- rev(best$source)
  dd <- d %>% dplyr::filter(group == g, source %in% best$source) %>% dplyr::mutate(source = factor(source, lev), ct = droplevels(ct))
  bl <- best %>% dplyr::mutate(source = factor(source, lev),
                               lab = if (RANK_BY == "top") paste0(fmt(rank_value, 1), " (", top_ct, ")") else fmt(rank_value, if (RANK_BY == "count") 0 else 1))
  nct <- nlevels(dd$ct)
  ggplot(dd, aes(ct, source)) + geom_tile(aes(fill=S), colour="white", linewidth=0.8) +
    geom_tile(data=dd %>% dplyr::filter(dir %in% dirn), fill=NA, colour="black", linewidth=0.7, width=0.92, height=0.88) +
    geom_text(aes(label=fmt(S, 0)), size=base / 3.8) +
    geom_text(data=bl, aes(x=nct + 0.75, y=source, label=lab), hjust=0, size=base / 3.6, inherit.aes=FALSE) +
    annotate("text", x=nct + 0.75, y=length(lev) + 0.85, label=RANK_HEAD, hjust=0, fontface="bold", size=base / 3.6) +
    scale_fill_gradient2(low=COL_AMP, mid="white", high=COL_BUF, limits=c(-SLIM, SLIM),
                         name="Buffering score\nblue = buffered\nred = amplified") +
    scale_x_discrete(position="top", expand=expansion(add=c(0.5, 2.4))) + coord_cartesian(clip="off") +
    labs(title=paste0("Top ", tolower(dirn), " TFs: ", nice_group(g)),
         subtitle="Buffering score in every cell type; box = called", x=NULL, y=NULL) +
    th_mn(base) + theme(axis.text.x=element_text(face="bold"))
}

# top_t_heat: same rows, aging (A) and stimulation (S) t-values in every cell type
plot_top_t_heat <- function(g, dirn, n=N_TOP, base=12) {
  x <- called %>% dplyr::filter(group == g, dir == dirn); if (!nrow(x)) return(NULL)
  best <- top_tfs(tf_summary(g, dirn), n)
  lev <- rev(best$source)
  dd <- d %>% dplyr::filter(group == g, source %in% best$source) %>% dplyr::mutate(source = factor(source, lev), ct = droplevels(ct))
  cts <- levels(dd$ct); step <- 2.4
  tl <- dd %>% dplyr::select(source, ct, dir, A=aging_act, S=ne_act) %>% tidyr::pivot_longer(c(A, S), names_to="k", values_to="t") %>%
    dplyr::mutate(xx = as.numeric(ct) * step + ifelse(k == "A", -0.5, 0.5))
  box <- dd %>% dplyr::filter(dir %in% dirn) %>% dplyr::distinct(source, ct) %>% dplyr::mutate(xx = as.numeric(ct) * step)
  ggplot(tl, aes(xx, source)) + geom_tile(aes(fill=t), colour="white", linewidth=0.6, width=1) +
    geom_tile(data=box, fill=NA, colour="black", linewidth=0.7, width=2, height=0.9) +
    geom_text(aes(label=fmt(t)), size=base / 4.6) + t_fill() +
    scale_x_continuous(breaks=seq_along(cts) * step, labels=cts, position="top",
                       sec.axis=dup_axis(breaks=c(rbind(seq_along(cts) * step - 0.5, seq_along(cts) * step + 0.5)),
                                         labels=rep(c("A", "S"), length(cts)))) +
    labs(title=paste0("Top ", tolower(dirn), " TFs: ", nice_group(g)),
         subtitle="Aging (A) and stimulation (S) t-values in every cell type; box = called", x=NULL, y=NULL) +
    th_mn(base) + theme(axis.text.x.top=element_text(face="bold"), axis.text.x.bottom=element_text(size=base - 3, colour="grey40"))
}

################################################################################
### DRAW
################################################################################

n_ct_of <- function(g) dplyr::n_distinct(d$ct[d$group == g])
for (g in groups) {
  gs <- gsub(" ", "_", g)
  for (dirn in c("Buffered", "Amplified")) {
    low <- tolower(dirn)
    if ("bubble_rank" %in% RUN)      { p <- plot_bubble_rank(g, dirn);      if (!is.null(p)) save_fig(p, "bubble_rank", paste0(low, "_", gs), 6.5, 0.3 * N_BUBBLE + 1.6) }
    if ("split_heat" %in% RUN)       { p <- plot_split_heat(g, dirn);       if (!is.null(p)) save_fig(p, "split_heat", paste0(low, "_", gs), 0.55 * n_ct_of(g) + 3, 0.32 * N_TOP + 1.6) }
    if ("lollipop_t" %in% RUN)       { p <- plot_lollipop_t(g, dirn);       if (!is.null(p)) save_fig(p, "lollipop_t", paste0(low, "_", gs), 7, 0.36 * N_TOP + 1.6) }
    if ("circles_t" %in% RUN)        { p <- plot_circles_t(g, dirn);        if (!is.null(p)) save_fig(p, "circles_t", paste0(low, "_", gs), 5.6, 0.4 * N_TOP + 1.4) }
    if ("percell_table" %in% RUN)    { p <- plot_percell_table(g, dirn);    if (!is.null(p)) save_fig(p, "percell_table", paste0(low, "_", gs), 2.3 * n_ct_of(g) + 1.5, 0.45 * N_PER_CT + 1.6) }
    if ("percell_lollipop" %in% RUN) { p <- plot_percell_lollipop(g, dirn); if (!is.null(p)) { k <- min(4, n_ct_of(g)); save_fig(p, "percell_lollipop", paste0(low, "_", gs), 3.8 * k + 1.2, ceiling(n_ct_of(g) / k) * (0.45 * N_PER_CT + 1.4) + 0.8) } }
    if ("top_score_heat" %in% RUN)   { p <- plot_top_score_heat(g, dirn);   if (!is.null(p)) save_fig(p, "top_score_heat", paste0(low, "_", gs), 0.6 * n_ct_of(g) + 4.6, 0.36 * N_TOP + 1.8) }
    if ("top_t_heat" %in% RUN)       { p <- plot_top_t_heat(g, dirn);       if (!is.null(p)) save_fig(p, "top_t_heat", paste0(low, "_", gs), 1.15 * n_ct_of(g) + 3.4, 0.36 * N_TOP + 1.8) }
  }
  if ("bubble_t" %in% RUN) for (dn in c("Buffered", "Amplified", "Both")) {
    p <- plot_bubble_t(g, dn); if (!is.null(p)) save_fig(p, "bubble_t", paste0(tolower(dn), "_", gs), 9, 0.32 * N_BUBBLE + 1.8) }
  if ("bubble_tp" %in% RUN && HAVE_P) for (dn in c("Buffered", "Amplified", "Both")) {
    p <- plot_bubble_t(g, dn, show_p=TRUE)
    if (!is.null(p)) save_fig(p + plot_annotation(caption=P_CAPTION), paste0("bubble_tp_", P_TYPE), paste0(tolower(dn), "_", gs), 9, 0.32 * N_BUBBLE + 2.1) }
  if ("bar_t" %in% RUN) { p <- plot_bar_t(g); if (!is.null(p)) save_fig(p, "bar_t", gs, 8, 0.36 * 2 * N_BAR + 1.8) }
  if ("bubble_rank" %in% RUN) { p <- plot_bubble_both(g); if (!is.null(p)) save_fig(p, "bubble_rank", paste0("both_", gs), 7, 0.3 * N_BUBBLE + 1.8) }
}

# all four groups in one figure (sexes side by side, middle above old; same scales)
grid4 <- function(fun, ...) {
  ps <- lapply(GROUP_ORDER, function(g) if (g %in% groups) { p <- fun(g, ...); if (is.null(p)) no_data(paste(nice_group(g), ": none")) else p } else no_data(paste(nice_group(g), ": no data")))
  wrap_plots(ps, ncol=2) + plot_layout(guides="collect")
}
if ("bubble_t" %in% RUN) for (dn in c("Buffered", "Amplified", "Both")) {
  p <- grid4(function(g, d2) plot_bubble_t(g, d2, n=if (d2 == "Both") 2 * N_BUBBLE_GROUPS else N_BUBBLE_GROUPS, base=11,
                                           title=nice_group(g), collect=FALSE), dn) +
    plot_annotation(title=paste0(if (dn == "Both") "Buffered and amplified" else dn, " TFs across sex and age"),
                    subtitle=paste0("Ranked by ", tolower(RANK_LAB), "; dot size = cell types called; t-values from the cell type with the highest score"),
                    theme=theme(plot.title=element_text(face="bold", size=16)))
  save_fig(p, "bubble_t", paste0("AllGroups_", tolower(dn)), 16, 2 * (0.32 * (N_BUBBLE_GROUPS * if (dn == "Both") 2 else 1) + 1.6))
}
if ("bubble_tp" %in% RUN && HAVE_P) for (dn in c("Buffered", "Amplified", "Both")) {
  p <- grid4(function(g, d2) plot_bubble_t(g, d2, n=if (d2 == "Both") 2 * N_BUBBLE_GROUPS else N_BUBBLE_GROUPS, base=11,
                                           title=nice_group(g), collect=FALSE, show_p=TRUE), dn) +
    plot_annotation(title=paste0(if (dn == "Both") "Buffered and amplified" else dn, " TFs across sex and age"),
                    subtitle=paste0("Ranked by ", tolower(RANK_LAB), "; dot size = cell types called; t-values from the cell type with the highest score"),
                    caption=P_CAPTION, theme=theme(plot.title=element_text(face="bold", size=16)))
  save_fig(p, paste0("bubble_tp_", P_TYPE), paste0("AllGroups_", tolower(dn)), 16, 2 * (0.32 * (N_BUBBLE_GROUPS * if (dn == "Both") 2 else 1) + 1.6) + 0.3)
}
if ("bubble_rank" %in% RUN) {
  for (dirn in c("Buffered", "Amplified")) {
    p <- grid4(function(g, dn) plot_bubble_rank(g, dn, n=N_BUBBLE_GROUPS, base=11, title=nice_group(g)), dirn) +
      plot_annotation(title=paste0(dirn, " TFs across sex and age"), subtitle=paste0("Ranked by ", tolower(RANK_LAB), "; dot size = cell types called"),
                      theme=theme(plot.title=element_text(face="bold", size=16)))
    save_fig(p, "bubble_rank", paste0("AllGroups_", tolower(dirn)), 11, 2 * (0.3 * N_BUBBLE_GROUPS + 1.4))
  }
  p <- grid4(function(g) plot_bubble_both(g, n=N_BUBBLE_GROUPS, base=11)) +
    plot_annotation(title="Buffered and amplified TFs across sex and age", subtitle=paste0("Ranked by ", tolower(RANK_LAB), "; dot size = cell types called"),
                    theme=theme(plot.title=element_text(face="bold", size=16)))
  save_fig(p, "bubble_rank", "AllGroups_both", 12, 2 * (0.3 * 2 * N_BUBBLE_GROUPS + 1.6))
}
if ("split_heat" %in% RUN) for (dirn in c("Buffered", "Amplified")) {
  p <- grid4(function(g, dn) plot_split_heat(g, dn, n=N_GROUPS, base=11, title=nice_group(g)), dirn) +
    plot_annotation(title=paste0(dirn, " TFs across sex and age"), theme=theme(plot.title=element_text(face="bold", size=16)))
  save_fig(p, "split_heat", paste0("AllGroups_", tolower(dirn)), 2 * (0.5 * max(sapply(groups, n_ct_of)) + 2.6), 2 * (0.32 * N_GROUPS + 1.4))
}
if ("lollipop_t" %in% RUN) for (dirn in c("Buffered", "Amplified")) {
  p <- grid4(function(g, dn) plot_lollipop_t(g, dn, n=N_GROUPS, base=11, title=nice_group(g)), dirn) +
    plot_annotation(title=paste0(dirn, " TFs across sex and age"), theme=theme(plot.title=element_text(face="bold", size=16)))
  save_fig(p, "lollipop_t", paste0("AllGroups_", tolower(dirn)), 13, 2 * (0.36 * N_GROUPS + 1.6))
}
if ("bar_t" %in% RUN) {
  p <- grid4(function(g) plot_bar_t(g, n=N_GROUPS, base=11)) +
    plot_annotation(title="TF buffering across sex and age",
                    subtitle=paste0("Top ", N_GROUPS, " buffered (right) and amplified (left) TF / cell-type pairs; same scales in every panel"),
                    theme=theme(plot.title=element_text(face="bold", size=16)))
  save_fig(p, "bar_t", "AllGroups", 15, 2 * (0.36 * 2 * N_GROUPS + 1.8))
}
################################################################################
### dreadd_recap: NE-defined TFs (old females) in DREADD-activated old females
###   x = NE stimulation t (old eNE vs old HC; this table), y = DREADD stimulation t (old C21 vs old PBS; held-out fits,
###   real labeling, from Buffering_DREADD_vs_NE units). No aging axis is used for the axes; NE's calls only pick the
###   coloured TFs (top DR_TOP_PER_CT buffered / amplified per cell type by buffering score).
###   "Same direction" = both stimulations change the TF's activity with the same sign.
###   Also reported (CSV only): buffered/amplified on DREADD's own aging axis = -(DREADD-cohort aging t x DREADD stim t).
### Output: Figure6_TF_Score/dreadd_recap/  scatter_<B|D>, scatter_bars_<B|D>              (stimulation t on both axes)
###         score_NEaxis_<B|D>, score_NEaxis_bars_<B|D>    (buffering scores; DREADD on eNE's aging axis)
###         score_OwnAxis_<B|D>, score_OwnAxis_bars_<B|D>  (buffering scores; DREADD on its own aging axis)
###         score_OwnAxis_buffered_<B|D>                   (own axis, eNE-buffered TFs only, DR_N_LABEL_BUF labels), bars_AllApproaches,
###         DREADD_Recap_Summary.csv, DREADD_Recap_SourceData.csv
################################################################################

if ("dreadd_recap" %in% RUN) {
  cat("\n-- dreadd_recap --\n")
  dr_root <- file.path(DR_BASE, "Buffering_DREADD_vs_NE"); if (!dir.exists(dr_root)) dr_root <- NA
  drt <- if (is.na(dr_root)) NULL else dplyr::bind_rows(lapply(DR_APPROACHES, function(ap) {
    fs <- list.files(file.path(dr_root, paste0("Approach_", ap), "units"), pattern="\\.RDS$", full.names=TRUE)
    dplyr::bind_rows(lapply(fs, function(f) { u <- readRDS(f)
      if (is.null(u$tf) || !nrow(u$tf) || !all(c("aging_act", "c21_act", "real") %in% names(u$tf))) return(NULL)
      u$tf %>% dplyr::filter(real) %>% dplyr::transmute(approach=ap, celltype=u$celltype, source,
                                                        DR_aging_t=aging_act, DR_stim_t=c21_act) })) }))
  if (is.null(drt) || !nrow(drt)) {
    cat("  DREADD TF results not found under", file.path(DR_BASE, "Buffering_DREADD_vs_NE"), "- skipped\n")
  } else {
    cat("  DREADD TF results:", dr_root, "\n")
    ne_f <- d %>% dplyr::filter(sex == "female", age == "old") %>%
      dplyr::select(celltype, ct, source, NE_aging_t=aging_act, NE_stim_t=ne_act, S, dir, score)
    rc <- drt %>% dplyr::inner_join(ne_f, by=c("celltype", "source")) %>%
      dplyr::group_by(approach, celltype) %>%
      dplyr::mutate(top=dplyr::case_when(dir == "Buffered"  & rank(-ifelse(dir %in% "Buffered",  score, -Inf)) <= DR_TOP_PER_CT ~ "Buffered",
                                         dir == "Amplified" & rank(-ifelse(dir %in% "Amplified", score, -Inf)) <= DR_TOP_PER_CT ~ "Amplified"),
                    top_buf=dir %in% "Buffered" & rank(-ifelse(dir %in% "Buffered", score, -Inf)) <= DR_TOP_PER_CT_BUF) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(cells=factor(ifelse(celltype %in% DR_NEURONAL, "Neurons", "Glia"), levels=c("Neurons", "Glia")),
                    cls=factor(dplyr::case_when(top == "Buffered" ~ "Top eNE-buffered", top == "Amplified" ~ "Top eNE-amplified",
                                                is.na(dir) ~ "Other TFs", TRUE ~ "Other called"),
                               levels=c("Top eNE-buffered", "Top eNE-amplified", "Other TFs", "Other called")),
                    same_dir=sign(NE_stim_t) == sign(DR_stim_t),
                    S_DR_own=-(DR_aging_t * DR_stim_t),             # DREADD scored on its own aging axis (old PBS vs young PBS)
                    S_DR_NEax=-(NE_aging_t * DR_stim_t),            # DREADD scored on eNE's aging axis (old HC vs young HC)
                    same_call_own=dplyr::case_when(top == "Buffered" ~ S_DR_own > 0, top == "Amplified" ~ S_DR_own < 0),
                    aging_dirs_agree=sign(NE_aging_t) == sign(DR_aging_t))
    summ <- rc %>% dplyr::filter(cls != "Other called") %>% dplyr::group_by(approach, cells, cls) %>%
      dplyr::summarise(n=dplyr::n(), pct_same_direction=round(100 * mean(same_dir), 1),
                       pct_same_call_own_axis=round(100 * mean(same_call_own), 1),
                       pct_aging_dirs_agree=round(100 * mean(aging_dirs_agree), 1),
                       binom_p=if (dplyr::n() > 0) signif(stats::binom.test(sum(same_dir), dplyr::n(), 0.5)$p.value, 2) else NA_real_,
                       .groups="drop")
    rdir <- file.path(OUT, "dreadd_recap"); dir.create(rdir, recursive=TRUE, showWarnings=FALSE)
    write.csv(summ, file.path(rdir, "DREADD_Recap_Summary.csv"), row.names=FALSE)
    write.csv(rc %>% dplyr::select(approach, celltype, source, cells, NE_call=dir, highlighted=top, buffering_score=S,
                                   NE_aging_t, NE_stim_t, DR_aging_t, DR_stim_t, same_dir, S_DR_NEax, S_DR_own, same_call_own, aging_dirs_agree),
              file.path(rdir, "DREADD_Recap_SourceData.csv"), row.names=FALSE)
    print(as.data.frame(summ), row.names=FALSE)

    RC_COL <- c(Buffered=COL_BUF, Amplified=COL_AMP)
    BAR_COL <- c(`Top eNE-buffered`=COL_BUF, `Top eNE-amplified`=COL_AMP, `Other TFs`="grey70")
    has_repel <- requireNamespace("ggrepel", quietly=TRUE)
    # mode: "stim"         x = eNE stim t,          y = DREADD stim t                 (shaded = same direction)
    #       "score_NEaxis" x = eNE buffering score, y = DREADD score on eNE's aging axis  (upper right = both buffer)
    #       "score_own"    x = eNE buffering score, y = DREADD score on its own aging axis
    rc_scatter <- function(ap, base=9, mode="stim", only_buf=FALSE, n_label=DR_N_LABEL) {
      x <- rc %>% dplyr::filter(approach == ap); if (!nrow(x)) return(NULL)
      clip <- function(v) pmax(pmin(v, DR_SCORE_CLIP), -DR_SCORE_CLIP)
      x <- switch(mode,
                  stim=x %>% dplyr::mutate(xv=NE_stim_t, yv=DR_stim_t),
                  score_NEaxis=x %>% dplyr::mutate(xv=clip(S), yv=clip(S_DR_NEax)),
                  score_own=x %>% dplyr::mutate(xv=clip(S), yv=clip(S_DR_own)))
      ylab <- switch(mode, stim="TF activity change with DREADD\n(t, old C21 vs old PBS)",
                     score_NEaxis="DREADD buffering score\n(eNE aging axis)", score_own="DREADD buffering score\n(own aging axis)")
      xlab <- if (mode == "stim") "TF activity change with eNE (t, old eNE vs old HC)" else "eNE buffering score"
      q_fill <- if (mode == "stim") c("grey60", "grey60") else c(COL_BUF, COL_AMP)
      hl <- if (only_buf) x %>% dplyr::filter(top_buf) %>% dplyr::mutate(top="Buffered")   # only_buf: top DR_TOP_PER_CT_BUF eNE-buffered
            else x %>% dplyr::filter(!is.na(top))
      lab <- hl %>% dplyr::group_by(cells, top) %>% dplyr::slice_max(abs(xv) + abs(yv), n=n_label, with_ties=FALSE) %>%
        dplyr::ungroup() %>% dplyr::mutate(l=paste0(source, " (", ct, ")"))
      p <- ggplot(x, aes(xv, yv)) +
        annotate("rect", xmin=0, xmax=Inf, ymin=0, ymax=Inf, fill=q_fill[1], alpha=0.08) +
        annotate("rect", xmin=-Inf, xmax=0, ymin=-Inf, ymax=0, fill=q_fill[2], alpha=0.08) +
        geom_hline(yintercept=0, colour="grey75", linewidth=0.3) + geom_vline(xintercept=0, colour="grey75", linewidth=0.3) +
        geom_point(colour="grey82", size=0.35) +
        geom_point(data=hl, aes(fill=top), shape=21, size=1.9, stroke=0.25) +
        scale_fill_manual(values=RC_COL, name=if (only_buf) paste0("Top ", DR_TOP_PER_CT_BUF, " eNE-buffered TFs per cell type") else paste0("Top ", DR_TOP_PER_CT, " eNE TFs per cell type"),
                          breaks=if (only_buf) "Buffered" else c("Buffered", "Amplified")) +
        scale_colour_manual(values=RC_COL, guide="none") +
        facet_wrap(~ cells) + th_cl(base) +
        theme(strip.background=element_rect(fill="grey93", colour=NA), strip.text=element_text(face="bold"), legend.position="bottom") +
        labs(x=xlab, y=ylab)
      if (has_repel) p <- p + ggrepel::geom_text_repel(data=lab, aes(label=l, colour=top), size=base / 4.4, fontface="italic",
                                                       max.overlaps=Inf, segment.size=0.15, min.segment.length=0, seed=2)
      p
    }
    rc_bars <- function(aps, base=9) {
      # two measures per NE class: same stimulation direction as eNE (solid), same call on each paradigm's own aging axis (light)
      MET <- c(pct_same_direction="Direction conserved in DREADD", pct_same_call_own_axis="Category conserved in DREADD")
      b <- summ %>% dplyr::filter(approach %in% aps) %>%
        tidyr::pivot_longer(c(pct_same_direction, pct_same_call_own_axis), names_to="metric", values_to="pct") %>%
        dplyr::filter(!is.na(pct)) %>%
        dplyr::mutate(cls=factor(sub("^Top ", "", cls), levels=sub("^Top ", "", names(BAR_COL))),
                      metric=factor(MET[metric], levels=MET))
      if (!nrow(b)) return(NULL)
      bcol <- setNames(BAR_COL, sub("^Top ", "", names(BAR_COL)))
      b <- b %>% dplyr::group_by(approach, cells, cls) %>%           # bar positions: pairs side by side, single bars centred
        dplyr::mutate(xpos=as.integer(cls) + if (dplyr::n() == 2) c(-0.2, 0.2)[as.integer(metric)] else 0) %>% dplyr::ungroup()
      p <- ggplot(b, aes(xpos, pct, fill=cls, alpha=metric)) +
        geom_col(width=0.38, colour="grey30", linewidth=0.2) +
        geom_hline(yintercept=50, linetype="dashed", colour="grey40", linewidth=0.3) +
        geom_text(aes(label=paste0(round(pct), "%")), vjust=-0.4, size=base / 4.4, alpha=1) +
        scale_x_continuous(breaks=seq_along(levels(b$cls)), labels=levels(b$cls)) +
        scale_fill_manual(values=bcol, guide="none") +
        scale_alpha_manual(values=c(1, 0.35), name=NULL, drop=FALSE) +
        scale_y_continuous(limits=c(0, 112), breaks=c(0, 25, 50, 75, 100), expand=c(0, 0)) +
        th_cl(base) + theme(axis.text.x=element_text(angle=35, hjust=1), strip.background=element_rect(fill="grey93", colour=NA),
                            strip.text=element_text(face="bold"), legend.position="bottom", legend.direction="vertical") +
        guides(alpha=guide_legend(override.aes=list(fill="grey40"))) +
        labs(x="Top eNE TFs per cell type", y="% of TF x cell types")
      if (length(aps) > 1) p + facet_grid(paste("Approach", approach) ~ cells) else p + facet_wrap(~ cells)
    }
    for (ap in unique(rc$approach)) {
      for (md in c("stim", "score_NEaxis", "score_own")) {
        p <- rc_scatter(ap, mode=md); if (is.null(p)) next
        tag <- c(stim="scatter", score_NEaxis="score_NEaxis", score_own="score_OwnAxis")[[md]]
        save_fig(p, "dreadd_recap", paste0(tag, "_", ap), 7, 4)
        save_fig((p | rc_bars(ap)) + plot_layout(widths=c(2.1, 1.2)), "dreadd_recap", paste0(tag, "_bars_", ap), 11, 4.4)
      }
    }
    # own-axis scatter, eNE-buffered TFs only, more labels
    for (ap in unique(rc$approach)) {
      p <- rc_scatter(ap, mode="score_own", only_buf=TRUE, n_label=DR_N_LABEL_BUF); if (is.null(p)) next
      save_fig(p, "dreadd_recap", paste0("score_OwnAxis_buffered_", ap), 7.5, 4.4)
    }
    save_fig(rc_bars(unique(rc$approach)), "dreadd_recap", "bars_AllApproaches", 4.6, 6)
  }
}

# TF-level ranking table (same ranking as the TF-level figures)
tf_tab <- dplyr::bind_rows(lapply(groups, function(g) dplyr::bind_rows(lapply(c("Buffered", "Amplified"), function(dn)
  tf_summary(g, dn) %>% dplyr::mutate(group = g, rank = dplyr::row_number())))))
write.csv(tf_tab %>% tidyr::separate(group, c("sex", "age"), sep=" ") %>%
            dplyr::select(direction=dir, sex, age, rank, source, n_celltypes=n_ct, summed_score=sum_score, top_score,
                          top_celltype=top_ct, top_aging_t, top_stim_t, top_aging_p=top_p_aging, top_stim_p=top_p_stim),
          file.path(OUT, paste0("TF_Ranking_by_", RANK_BY, ".csv")), row.names=FALSE)
cat("\nDone:", OUT, "\n")
