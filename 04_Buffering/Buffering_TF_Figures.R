################################################################################
### Buffering_TF_Figures.R
###
### NE: which TFs are buffered / amplified at middle and old age, and how each changes
###     (maintained, lost, age-emergent, amplified), in which cell types; each sex separately.
###   T1  Trajectory map: TF (rows) x cell type (columns), colour = middle -> old trajectory
###   T2  Middle vs old side by side: dot = TF score S (colour), filled when called
###   T3  Counts: number of TFs per trajectory in each cell type
###   T4  Across cell types: each TF's position at middle and old age (arrow), shortlist
###
### DREADD vs NE (old females): do TFs move the same way?
###   D1  Agreement map: TF x cell type, colour = NE direction x whether DREADD matches
###       (dark = DREADD also called it at rank <= 5 of 20; light = same sign only;
###        orange = DREADD goes the opposite way), Approach B and D
###   D2  Per cell type: share of NE-called TFs that DREADD moves the same way (stacked bar)
###
### Inputs: Buffering_Candidates_Master/ (Step5, Step6 CSVs; Buffering_Pipeline.R step 2)
###         Buffering_MiddleAge_Robustness/ (Check 1 / Check 2 CSVs; Buffering_Pipeline.R step 1)
###         Buffering_DREADD_vs_NE/Approach_{B,D}/Part4_TFActivity_PerCellType_withNE.csv
### Outputs: Paper_Figures/<date>/TF/ (PDF + PNG per panel, source-data CSVs)
################################################################################

suppressPackageStartupMessages({
  library(conflicted); library(dplyr); library(tidyr); library(ggplot2)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::mutate, .quiet=TRUE)
})

NE_BASE <- "PATH_TO_NE_OUTPUT"
DR_BASE <- "PATH_TO_DREADD_OUTPUT"
NE_MASTER <- file.path(NE_BASE, "Buffering_Candidates_Master")
DR_VS_NE  <- file.path(DR_BASE, "Buffering_DREADD_vs_NE")
FIG_DATE  <- format(Sys.Date(), "%Y-%m-%d")   # set to an earlier date to fill in that folder
OVERWRITE <- FALSE                              # FALSE = skip panels whose PDF already exists
OUT <- file.path(NE_BASE, "Paper_Figures", FIG_DATE, "TF")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

APPROACHES <- c("B","D")
ACT_RANK   <- 5               # call = rank <= 5 of 20 (as in the pipeline)
MIN_CT     <- 2               # T1/D1: TFs with a call in >= this many cell types
MAX_TFS    <- 45              # T1/D1: max rows per panel (most cell types first)
SHORTLIST  <- c("Nr3c1","Creb1","Sirt1","Atf6","Ddit3","Atf4","Foxo1","Foxo3","Nfkb1","Nfkb2","Rela","Relb",
                "Jun","Fos","Egr1","Egr3","Hsf1")
ALWAYS_SHOW <- c("Nr3c1","Creb1","Sirt1","Atf6")      # kept in T1/D1 even below MIN_CT

CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN",
              "Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell","Macrophage","T Cell")
CT_SHORT <- c("Astrocytes"="Astro","Oligodendrocytes"="Oligo","Cajal-Retzius Cell"="CR cell",
              "Endothelial"="Endo","Microglia"="Micro","Pericyte"="Peri")
short_ct <- function(x) ifelse(x %in% names(CT_SHORT), CT_SHORT[x], x)
SEX_LAB <- c(female="Females", male="Males")

TRAJ_COL <- c("Maintained"="#1B5E20", "Weakened"="#81C784",
              "Lost (flipped to amplified)"="#C62828", "Lost"="#EF9A9A",
              "Age-emergent"="#1565C0",
              "Amplified, both ages"="#4A148C", "Amplified at old age"="#AB47BC",
              "Amplified at middle only"="#D1C4E9")
AGR_COL <- c("Buffered in both (DREADD called)"="#1B5E20", "Buffered in both (same sign)"="#A5D6A7",
             "Amplified in both (DREADD called)"="#4A148C", "Amplified in both (same sign)"="#CE93D8",
             "DREADD opposite"="#F57C00")

th <- theme_classic(base_size=7) +
  theme(strip.background=element_rect(fill="grey92", color="black", linewidth=0.3),
        strip.text=element_text(face="bold", size=7), axis.line=element_line(linewidth=0.3),
        axis.ticks=element_line(linewidth=0.3), plot.title=element_text(face="bold", hjust=0.5, size=8),
        plot.subtitle=element_text(color="grey35", hjust=0.5, size=6), legend.key.size=unit(3, "mm"))
mm <- function(x) x / 25.4
sv <- function(p, name, w, h) {
  if (!OVERWRITE && file.exists(file.path(OUT, paste0(name, ".pdf")))) { cat("  exists, skipped:", name, "\n"); return(invisible()) }
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width=mm(w), height=mm(h), limitsize=FALSE)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width=mm(w), height=mm(h), dpi=600, bg="white", limitsize=FALSE)
  cat("  saved:", name, "\n")
}
wsrc <- function(x, name) write.csv(x, file.path(OUT, paste0("SourceData_", name, ".csv")), row.names=FALSE)
grid_theme <- theme(axis.text.x=element_text(angle=45, hjust=1), axis.line=element_blank(),
                    panel.border=element_rect(fill=NA, color="grey40", linewidth=0.3),
                    panel.grid.major=element_line(color="grey90", linewidth=0.2))

# pick rows: TFs called in >= MIN_CT cell types (+ ALWAYS_SHOW), most cell types first, capped
pick_tfs <- function(d, n_col="n_called") {
  r <- d %>% dplyr::group_by(source) %>% dplyr::summarise(n=sum(.data[[n_col]]), .groups="drop") %>%
    dplyr::filter(n >= MIN_CT | source %in% ALWAYS_SHOW) %>% dplyr::arrange(dplyr::desc(n), source)
  head(unique(c(intersect(ALWAYS_SHOW, r$source), r$source)), MAX_TFS)
}

################################################################################
### NE
################################################################################

traj <- read.csv(file.path(NE_MASTER, "Step6_Trajectory_PerCellType.csv"), stringsAsFactors=FALSE) %>%
  dplyr::mutate(trajectory=ifelse(trajectory == "Lost" & flipped_to_amplified %in% TRUE,
                                  "Lost (flipped to amplified)", trajectory))
act  <- read.csv(file.path(NE_MASTER, "Step5_TFActivity_PerCellType.csv"), stringsAsFactors=FALSE)
tfl  <- read.csv(file.path(NE_MASTER, "Step6_Trajectory_TFlevel.csv"), stringsAsFactors=FALSE)
wsrc(traj, "NE_TF_Trajectory_PerCellType")

for (sx in c("female","male")) {
  cat("\n== NE,", SEX_LAB[[sx]], "==\n")

  # T1: trajectory map
  d <- traj %>% dplyr::filter(sex == sx, trajectory != "Neither") %>% dplyr::mutate(n_called=1L)
  rows <- pick_tfs(d)
  d1 <- d %>% dplyr::filter(source %in% rows) %>%
    dplyr::mutate(strength=pmax(abs(S_mid), abs(S_old), na.rm=TRUE),
                  strength=pmin(strength, quantile(strength, 0.95, na.rm=TRUE)),
                  trajectory=factor(trajectory, levels=names(TRAJ_COL)),
                  celltype=factor(celltype, levels=intersect(CT_ORDER, unique(traj$celltype[traj$sex == sx]))),
                  source=factor(source, levels=rev(rows)))
  stat <- traj %>% dplyr::filter(sex == sx) %>% dplyr::distinct(celltype, ct_status)
  xlab <- setNames(sprintf("%s\n(%s)", short_ct(stat$celltype),
                           dplyr::recode(stat$ct_status, `buffers at both ages`="both", `lost at old age`="lost at old",
                                         `gained at old age`="gained at old")), stat$celltype)
  p <- ggplot(d1, aes(celltype, source)) +
    geom_point(aes(fill=trajectory, size=strength), shape=21, color="black", stroke=0.25) +
    scale_fill_manual(values=TRAJ_COL, name="Middle -> old", drop=FALSE) +
    scale_size_continuous(range=c(0.8, 3.2), name="|S|\n(stronger age)") +
    guides(fill=guide_legend(override.aes=list(size=2.5))) +
    scale_x_discrete(labels=xlab, drop=FALSE) + th + grid_theme +
    labs(title=paste0("TF buffering, middle -> old age (NE, ", SEX_LAB[[sx]], ")"),
         subtitle=sprintf("TFs called (rank <= %d of 20) in >= %d cell types; column label = cell-type buffering status",
                          ACT_RANK, MIN_CT), x=NULL, y=NULL)
  sv(p, paste0("T1_TrajectoryMap_", sx), 130, 30 + 3.2 * length(rows))

  # T2: middle vs old side by side (score + call)
  d2 <- act %>% dplyr::filter(sex == sx, source %in% rows) %>%
    dplyr::mutate(call=dplyr::case_when(buf ~ "Buffered", amp ~ "Amplified", TRUE ~ "Not called"),
                  age=factor(ifelse(age == "middle", "Middle age", "Old age"), levels=c("Middle age","Old age")),
                  celltype=factor(celltype, levels=levels(d1$celltype)), source=factor(source, levels=rev(rows)),
                  Sc=pmax(pmin(S, quantile(abs(S), 0.95, na.rm=TRUE)), -quantile(abs(S), 0.95, na.rm=TRUE)))
  lim <- max(abs(d2$Sc), na.rm=TRUE)
  p <- ggplot(d2, aes(celltype, source)) +
    geom_point(aes(fill=Sc, size=call, color=call), shape=21, stroke=0.35) +
    scale_fill_gradient2(low="#762A83", mid="white", high="#1B7837", limits=c(-lim, lim),
                         name="Buffering score S\n(> 0 buffered,\n< 0 amplified)") +
    scale_size_manual(values=c(Buffered=2.2, Amplified=2.2, `Not called`=1.0), name=NULL) +
    scale_color_manual(values=c(Buffered="black", Amplified="black", `Not called`="grey75"), name=NULL) +
    scale_x_discrete(labels=short_ct) +
    facet_grid(. ~ age) + th + grid_theme +
    labs(title=paste0("TF buffering at middle and old age (NE, ", SEX_LAB[[sx]], ")"),
         subtitle=sprintf("Large outlined dot = called (rank <= %d of 20 relabelings)", ACT_RANK), x=NULL, y=NULL)
  sv(p, paste0("T2_MiddleVsOld_", sx), 183, 30 + 3.2 * length(rows))

  # T3: counts per trajectory per cell type
  d3 <- d %>% dplyr::count(celltype, trajectory) %>%
    dplyr::mutate(trajectory=factor(trajectory, levels=rev(names(TRAJ_COL))),
                  celltype=factor(celltype, levels=rev(levels(d1$celltype))))
  p <- ggplot(d3, aes(n, celltype, fill=trajectory)) + geom_col(width=0.75) +
    scale_fill_manual(values=TRAJ_COL, name=NULL, breaks=names(TRAJ_COL)) + th +
    labs(title=paste0("TFs per trajectory (NE, ", SEX_LAB[[sx]], ")"), x="Number of TFs", y=NULL)
  sv(p, paste0("T3_TrajectoryCounts_", sx), 110, 20 + 5 * nlevels(d3$celltype))
  wsrc(d3, paste0("T3_", sx))
}

# T4: across cell types, middle -> old position (shortlist), both sexes
d4 <- tfl %>% dplyr::filter(source %in% SHORTLIST) %>%
  dplyr::mutate(sex=factor(SEX_LAB[sex], levels=SEX_LAB), source=factor(source, levels=rev(SHORTLIST)))
p <- ggplot(d4, aes(y=source)) +
  geom_vline(xintercept=0, color="grey60", linewidth=0.3) +
  geom_segment(aes(x=pos_mid, xend=pos_old, yend=source), color="grey55", linewidth=0.4,
               arrow=arrow(length=unit(1.2, "mm"), type="closed")) +
  geom_point(aes(x=pos_mid), shape=21, fill="white", size=1.6, stroke=0.4) +
  geom_point(aes(x=pos_old, fill=label), shape=21, size=2, stroke=0.3) +
  scale_fill_manual(values=c(TRAJ_COL, `Not buffered`="grey85"), name="Middle -> old") +
  scale_x_continuous(limits=c(-1, 1), breaks=c(-1, 0, 1), labels=c("most\namplified", "0", "most\nbuffered")) +
  facet_wrap(~sex) + th +
  labs(title="Shortlist TFs across cell types (NE)",
       subtitle="Position among 20 labelings: open = middle age, filled = old age", x=NULL, y=NULL)
sv(p, "T4_Shortlist_AcrossCellTypes", 150, 85)
wsrc(d4, "T4")

################################################################################
### DREADD vs NE (old females)
################################################################################

cat("\n== DREADD vs NE ==\n")
dr <- dplyr::bind_rows(lapply(APPROACHES, function(ap)
  read.csv(file.path(DR_VS_NE, paste0("Approach_", ap), "Part4_TFActivity_PerCellType_withNE.csv"), stringsAsFactors=FALSE) %>%
    dplyr::mutate(approach=ap))) %>%
  dplyr::filter(!is.na(NE_S), NE_call %in% c("buffered","amplified")) %>%
  dplyr::mutate(same=sign(S_NEax) == sign(NE_S), dr_called=DR_call_NEax == NE_call,
                agree=dplyr::case_when(
                  NE_call == "buffered"  & dr_called ~ "Buffered in both (DREADD called)",
                  NE_call == "buffered"  & same      ~ "Buffered in both (same sign)",
                  NE_call == "amplified" & dr_called ~ "Amplified in both (DREADD called)",
                  NE_call == "amplified" & same      ~ "Amplified in both (same sign)",
                  TRUE ~ "DREADD opposite"))
wsrc(dr %>% dplyr::select(approach, celltype, source, NE_call, NE_S, S_NEax, DR_call_NEax, rank_buf_NEax, rank_amp_NEax, agree),
     "DREADD_vs_NE_TFs")

# D1: agreement map, rows = TFs NE called in >= MIN_CT cell types
rows <- pick_tfs(dr %>% dplyr::filter(approach == APPROACHES[1]) %>% dplyr::mutate(n_called=1L))
d5 <- dr %>% dplyr::filter(source %in% rows) %>%
  dplyr::mutate(strength=pmin(abs(S_NEax), quantile(abs(S_NEax), 0.95, na.rm=TRUE)),
                agree=factor(agree, levels=names(AGR_COL)), source=factor(source, levels=rev(rows)),
                celltype=factor(celltype, levels=intersect(CT_ORDER, unique(dr$celltype))),
                approach=paste("DREADD Approach", approach))
pct <- d5 %>% dplyr::group_by(approach) %>% dplyr::summarise(p=round(100 * mean(agree != "DREADD opposite")), .groups="drop")
p <- ggplot(d5, aes(celltype, source)) +
  geom_point(aes(fill=agree, size=strength), shape=21, color="black", stroke=0.25) +
  scale_fill_manual(values=AGR_COL, name=NULL, drop=FALSE) +
  scale_size_continuous(range=c(0.8, 3.2), name="|DREADD score|") +
  guides(fill=guide_legend(override.aes=list(size=2.5))) +
  scale_x_discrete(labels=short_ct) + facet_grid(. ~ approach) + th + grid_theme +
  labs(title="DREADD moves NE's TFs the same way (old females)",
       subtitle=paste0("Tiles = TF x cell type called by NE; same direction: ",
                       paste(sprintf("%s %d%%", sub("DREADD ", "", pct$approach), pct$p), collapse=", ")),
       x=NULL, y=NULL)
sv(p, "D1_TFAgreementMap", 183, 35 + 3.2 * length(rows))

# D2: per cell type, share of NE-called TFs by agreement class
d6 <- dr %>% dplyr::count(approach, celltype, agree) %>% dplyr::group_by(approach, celltype) %>%
  dplyr::mutate(pct=100 * n / sum(n), n_tot=sum(n)) %>% dplyr::ungroup() %>%
  dplyr::filter(n_tot >= 20) %>%
  dplyr::mutate(agree=factor(agree, levels=rev(names(AGR_COL))),
                celltype=factor(celltype, levels=rev(intersect(CT_ORDER, unique(celltype)))),
                approach=paste("DREADD Approach", approach))
p <- ggplot(d6, aes(pct, celltype, fill=agree)) + geom_col(width=0.75) +
  geom_vline(xintercept=50, linetype="dashed", color="grey30", linewidth=0.3) +
  scale_fill_manual(values=AGR_COL, name=NULL, breaks=names(AGR_COL)) +
  facet_wrap(~approach) + th + theme(legend.position="bottom") + guides(fill=guide_legend(ncol=2)) +
  labs(title="Share of NE-called TFs that DREADD moves the same way",
       subtitle="Green/purple = same direction as NE; orange = opposite", x="% of NE-called TFs", y=NULL)
sv(p, "D2_TFAgreementShare", 150, 35 + 5 * nlevels(d6$celltype))
wsrc(d6, "D2")

################################################################################
### HEATMAP-STYLE PANELS (old Figure 6 layout)
###   H1  Buffering TFs per cell type (tiles; number = rank among 20 labelings)
###         NE old female, NE old male, DREADD B / D (NE aging axis)
###   H2  Superpooled (summed across cell types): top buffered and amplified TFs, ranked
###         NE (each sex, each age) and DREADD B / D
###   H3  NE middle age only, females and males (tiles)
###   H4  Cosine trajectory middle -> old, females and males (held-out cosine)
###   H5  Lost buffering middle -> old (NE), females and males (tiles)
###         + the same layout for maintained and age-emergent TFs
### DREADD TF calls use NE's aging direction (S_NEax): its own-axis ranks are not valid.
################################################################################

cat("\n== Heatmap-style panels ==\n")
H_ROWS     <- 20                 # TFs per heatmap (most cell types first)
COL_SEX    <- c(female="#D35A3B", male="#3A62C4")
COL_DR     <- "#8E44AD"
COL_AMP_T  <- "grey35"
latest <- function(dir, pat) {
  f <- list.files(dir, pattern=pat, recursive=TRUE, full.names=TRUE)
  if (!length(f)) stop("No file matching '", pat, "' under ", dir)
  f[order(file.mtime(f), decreasing=TRUE)][1]
}
ct_cols_of <- function(x) intersect(CT_ORDER, unique(x))

# d: celltype, source, call ("Buffered"/"Amplified"/"ns"), rank (shown in called tiles)
hm_calls <- function(d, title, col_buf, ct_all, focus="Buffered", subtitle=NULL) {
  rows <- d %>% dplyr::filter(call == focus) %>% dplyr::count(source) %>%
    dplyr::arrange(dplyr::desc(n), source) %>% head(H_ROWS) %>% dplyr::pull(source)
  if (!length(rows)) { cat("  (no", focus, "TFs for:", title, ")\n"); return(NULL) }
  dd <- d %>% dplyr::filter(source %in% rows) %>% dplyr::select(celltype, source, call, rank) %>%
    tidyr::complete(celltype=ct_all, source=rows, fill=list(call="Not analysed", rank=NA)) %>%
    dplyr::mutate(celltype=factor(celltype, levels=ct_all), source=factor(source, levels=rev(rows)),
                  call=factor(call, levels=c("Buffered","Amplified","ns","Not analysed")),
                  lab=ifelse(call %in% c("Buffered","Amplified") & !is.na(rank), rank, ""))
  ggplot(dd, aes(celltype, source, fill=call)) + geom_tile(color="white", linewidth=0.4) +
    geom_text(aes(label=lab), size=1.7, color="white") +
    scale_fill_manual(values=c(Buffered=col_buf, Amplified=COL_AMP_T, ns="grey90", `Not analysed`="white"),
                      name=NULL, drop=FALSE) +
    scale_x_discrete(labels=short_ct) + th +
    theme(axis.text.x=element_text(angle=45, hjust=1), axis.line=element_blank(), legend.position="top",
          axis.text.y=element_text(face="italic")) +
    labs(title=title, subtitle=subtitle, x=NULL, y=NULL)
}
h_size <- function(p) if (is.null(p)) NULL else 30 + 3.4 * length(unique(p$data$source))
ne_calls_ct <- function(ag, sx) act %>% dplyr::filter(age == ag, sex == sx) %>%
  dplyr::mutate(call=dplyr::case_when(buf ~ "Buffered", amp ~ "Amplified", TRUE ~ "ns"),
                rank=dplyr::case_when(buf ~ rank_buf, amp ~ rank_amp, TRUE ~ NA_integer_))
SUB_RANK <- sprintf("Called = rank <= %d of 20 relabelings; number = rank (1 = strongest); white = cell type not analysed", ACT_RANK)

# H1 + H3: NE per cell type (old and middle, each sex)
for (ag in c("old","middle")) for (sx in c("female","male")) {
  d <- ne_calls_ct(ag, sx)
  p <- hm_calls(d, sprintf("TF buffering NE, %s age, %s", ifelse(ag == "old", "old", "middle"), SEX_LAB[[sx]]),
                COL_SEX[[sx]], ct_cols_of(act$celltype[act$sex == sx]), subtitle=SUB_RANK)
  if (!is.null(p)) sv(p, sprintf("%s_NE_%s_%s", ifelse(ag == "old", "H1", "H3"), ag, sx), 120, h_size(p))
}

# H1: DREADD per cell type (NE aging axis)
for (ap in APPROACHES) {
  d <- read.csv(file.path(DR_VS_NE, paste0("Approach_", ap), "Part4_TFActivity_PerCellType_withNE.csv"), stringsAsFactors=FALSE) %>%
    dplyr::mutate(call=dplyr::case_when(DR_call_NEax == "buffered" ~ "Buffered", DR_call_NEax == "amplified" ~ "Amplified", TRUE ~ "ns"),
                  rank=dplyr::case_when(call == "Buffered" ~ rank_buf_NEax, call == "Amplified" ~ rank_amp_NEax, TRUE ~ NA_integer_))
  p <- hm_calls(d, paste0("TF buffering DREADD (Approach ", ap, "), old females"), COL_DR, ct_cols_of(d$celltype),
                subtitle=paste(SUB_RANK, "| scored against NE's aging direction"))
  if (!is.null(p)) sv(p, paste0("H1_DREADD_Approach", ap), 120, h_size(p))
}

# H2: superpooled, ranked (top buffered and top amplified)
pooled_bar <- function(d, title, col_buf) {
  b <- d %>% dplyr::filter(call == "Buffered") %>% dplyr::slice_max(S, n=15, with_ties=FALSE)
  a <- d %>% dplyr::filter(call == "Amplified") %>% dplyr::slice_min(S, n=15, with_ties=FALSE)
  dd <- dplyr::bind_rows(b, a) %>% dplyr::arrange(S) %>% dplyr::mutate(source=factor(source, levels=unique(source)))
  if (!nrow(dd)) return(NULL)
  ggplot(dd, aes(S, source, fill=call)) + geom_col(width=0.7) +
    geom_vline(xintercept=0, color="grey30", linewidth=0.3) +
    geom_text(aes(label=paste0(rank, "/20"), hjust=ifelse(S > 0, -0.15, 1.15)), size=1.7) +
    scale_fill_manual(values=c(Buffered=col_buf, Amplified=COL_AMP_T), name=NULL) +
    scale_x_continuous(expand=expansion(mult=0.18)) + th +
    theme(axis.text.y=element_text(face="italic"), legend.position="top") +
    labs(title=title, subtitle="Summed across cell types; label = rank among 20 labelings", x="Summed buffering score", y=NULL)
}
all_ne <- read.csv(file.path(NE_MASTER, "Step5_TFActivity_AcrossCellTypes.csv"), stringsAsFactors=FALSE)
for (ag in c("old","middle")) for (sx in c("female","male")) {
  d <- all_ne %>% dplyr::filter(age == ag, sex == sx) %>%
    dplyr::mutate(S=S_sum, call=dplyr::case_when(buf ~ "Buffered", amp ~ "Amplified", TRUE ~ "ns"),
                  rank=ifelse(S > 0, rank_buf, rank_amp))
  p <- pooled_bar(d, sprintf("Superpooled TF buffering NE, %s age, %s", ag, SEX_LAB[[sx]]), COL_SEX[[sx]])
  if (!is.null(p)) sv(p, sprintf("H2_Pooled_NE_%s_%s", ag, sx), 89, 30 + 3.2 * nrow(p$data))
}
for (ap in APPROACHES) {
  d <- read.csv(file.path(DR_VS_NE, paste0("Approach_", ap), "Part4_TFActivity_AcrossCellTypes.csv"), stringsAsFactors=FALSE) %>%
    dplyr::mutate(S=S_NEax_sum,
                  call=dplyr::case_when(S > 0 & rank_buf_NEax <= ACT_RANK ~ "Buffered", S < 0 & rank_amp_NEax <= ACT_RANK ~ "Amplified", TRUE ~ "ns"),
                  rank=ifelse(S > 0, rank_buf_NEax, rank_amp_NEax))
  p <- pooled_bar(d, paste0("Superpooled TF buffering DREADD (Approach ", ap, ")"), COL_DR)
  if (!is.null(p)) sv(p, paste0("H2_Pooled_DREADD_Approach", ap), 89, 30 + 3.2 * nrow(p$data))
}

# H4: cosine trajectory middle -> old (held-out), each sex
rob_dir <- file.path(NE_BASE, "Buffering_MiddleAge_Robustness")
c1 <- read.csv(latest(rob_dir, "_Check1_BiasCorrected\\.csv$"), stringsAsFactors=FALSE) %>% dplyr::mutate(age=tolower(age))
c2 <- read.csv(latest(rob_dir, "_Check2_HeldOut_PerMouse\\.csv$"), stringsAsFactors=FALSE) %>% dplyr::mutate(age=tolower(age))
cosd <- c2 %>% dplyr::group_by(age, sex, celltype) %>% dplyr::summarise(cos=mean(cosine, na.rm=TRUE), .groups="drop") %>%
  dplyr::inner_join(c1 %>% dplyr::select(age, sex, celltype, n_genes, beyond_null), by=c("age","sex","celltype")) %>%
  dplyr::filter(n_genes >= 1000) %>%
  dplyr::mutate(pass=cos < 0 & beyond_null %in% TRUE, age=factor(ifelse(age == "middle", "Middle", "Old"), levels=c("Middle","Old")),
                sex=factor(SEX_LAB[sex], levels=SEX_LAB))
wsrc(cosd, "H4_CosineTrajectory")
CT_COLS <- c("DG ExNs"="#E53935","CA1 ExNs"="#CC0000","CA2 ExNs"="#8E24AA","CA3 ExNs"="#5E35B1","Sub ExNs"="#80D8FF",
             "IN"="#00897B","Astrocytes"="#43A047","Oligodendrocytes"="#1E88E5","OPC"="#00ACC1","Microglia"="#7CB342",
             "Endothelial"="#FB8C00","Pericyte"="#6D4C41","APC"="#546E7A","Stromal"="#27408B","Macrophage"="#C62828",
             "Cajal-Retzius Cell"="#FDD835","T Cell"="#757575")
lab_end <- cosd %>% dplyr::filter(age == "Old") %>% dplyr::mutate(lab=short_ct(celltype))
p <- ggplot(cosd, aes(age, cos, group=celltype, color=celltype)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40", linewidth=0.3) +
  geom_line(linewidth=0.45, alpha=0.85) +
  geom_point(aes(fill=ifelse(pass, celltype, NA)), shape=21, size=1.8, stroke=0.5) +
  scale_color_manual(values=CT_COLS, guide="none") +
  scale_fill_manual(values=CT_COLS, na.value="white", guide="none") +
  facet_wrap(~sex) + th +
  labs(title="Cosine trajectory: middle -> old age (NE)",
       subtitle="Held-out cosine; filled = passes the buffering gate, open = does not",
       x=NULL, y="Cosine similarity (aging, stimulation)")
if (requireNamespace("ggrepel", quietly=TRUE)) {
  p <- p + ggrepel::geom_text_repel(data=lab_end, aes(label=lab), size=1.8, nudge_x=0.25, direction="y", hjust=0,
                                    segment.size=0.15, show.legend=FALSE)
} else {
  p <- p + geom_text(data=lab_end, aes(label=lab), size=1.8, hjust=-0.2, show.legend=FALSE)
}
sv(p + scale_x_discrete(expand=expansion(add=c(0.3, 0.9))), "H4_CosineTrajectory", 140, 80)

# H5: lost buffering (and maintained / age-emergent), each sex
TR_SETS <- list(
  Lost        = list(cats=c("Lost (flipped to amplified)","Lost","Weakened"), title="Lost buffering: middle -> old age"),
  Maintained  = list(cats=c("Maintained"), title="Maintained buffering: middle and old age"),
  AgeEmergent = list(cats=c("Age-emergent"), title="Age-emergent buffering: old age only"))
for (nm in names(TR_SETS)) for (sx in c("female","male")) {
  cats <- TR_SETS[[nm]]$cats
  d <- traj %>% dplyr::filter(sex == sx)
  rows <- d %>% dplyr::filter(trajectory %in% cats) %>% dplyr::count(source) %>%
    dplyr::arrange(dplyr::desc(n), source) %>% head(H_ROWS) %>% dplyr::pull(source)
  if (!length(rows)) next
  ct_all <- ct_cols_of(d$celltype)
  cols <- c(setNames(c(COL_SEX[[sx]], grDevices::adjustcolor(COL_SEX[[sx]], alpha.f=0.55),
                       grDevices::adjustcolor(COL_SEX[[sx]], alpha.f=0.25))[seq_along(cats)], cats),
            Other="grey90", `Not analysed`="white")
  if (nm == "Lost") cols[["Lost (flipped to amplified)"]] <- "#7F1D1D"
  dd <- d %>% dplyr::filter(source %in% rows) %>%
    dplyr::mutate(cat=ifelse(trajectory %in% cats, trajectory, "Other")) %>% dplyr::select(celltype, source, cat) %>%
    tidyr::complete(celltype=ct_all, source=rows, fill=list(cat="Not analysed")) %>%
    dplyr::mutate(celltype=factor(celltype, levels=ct_all), source=factor(source, levels=rev(rows)),
                  cat=factor(cat, levels=names(cols)))
  p <- ggplot(dd, aes(celltype, source, fill=cat)) + geom_tile(color="white", linewidth=0.4) +
    scale_fill_manual(values=cols, name=NULL, drop=FALSE) + scale_x_discrete(labels=short_ct) + th +
    theme(axis.text.x=element_text(angle=45, hjust=1), axis.line=element_blank(), legend.position="top",
          axis.text.y=element_text(face="italic")) +
    labs(title=paste0(TR_SETS[[nm]]$title, " (NE, ", SEX_LAB[[sx]], ")"),
         subtitle="Cell types buffering at either age; white = not analysed", x=NULL, y=NULL)
  sv(p, paste0("H5_", nm, "_", sx), 120, 30 + 3.4 * length(rows))
}

cat("\n=== DONE ===\nOutputs in:", OUT, "\n")
