################################################################################
### TF_Buffering_Significance.R
###
### TF buffering called by SIGNIFICANCE (field-standard decoupleR use), instead of the buffering score + rank.
### Same design as the master step of Buffering_Pipeline.R: same pseudobulk, same groups, same held-out fits
### (B_k = aged HC without mouse k vs young HC; E_k = aged eNE vs aged HC mouse k), same DESeq2 settings, same signed
### CollecTRI network (UniProt-ID TFs dropped, as in the master), same cell types (master Step 1 gate).
###
### Per TF, cell type, sex, age and held-out fit k:
###   decoupleR::run_ulm on the DESeq2 Wald statistics of all expressed genes -> aging activity (t, p) and
###   stimulation activity (t, p)
###   significant in that fit = ULM p < SIG_P for that contrast (SIG_MODE = "p"), or BH FDR < SIG_P across TFs
###   within that fit and contrast (SIG_MODE = "fdr"); both are written, SIG_MODE picks the one used for figures
###   buffered in fit k  = aging significant AND stimulation significant AND opposite signs
###   amplified in fit k = both significant AND same sign
### Tiers (mirroring the gene tiers):
###   TF-T1 = buffered in 3 of 3 fits      TF-T3 = buffered in >= 2 of 3 fits      (same for amplified: TF-A1 / TF-A3)
### Chance check (RUN_NULL): the same calls on all 19 alternative assignments of the six aged mice to HC / eNE.
###   Per cell type: number of buffered TFs, real vs the 19 relabelings (empirical p).
###   Per TF: how often it is called in the relabelings (a TF called in many relabelings is not specific).
### Trajectories (middle -> old, per TF x cell type): maintained / lost / age-emergent / amplified at old age.
### Comparison with the current calls (master Step5, rank-based) is written as a table.
###
### Outputs: <base_path>/Buffering_Candidates_Master/TF_Significance/
###   units/<age>__<sex>__<celltype>.RDS   cached per-fit ULM results (real + relabelings); resumable, shared by both thresholds
###   fdr/ and p/                          everything below, once per threshold (SIG_MODES)
###   TFSig_PerFit.csv                     real labeling, every TF x fit with activities, p, FDR
###   TFSig_Calls.csv                      real labeling, per TF x cell type: fits passing, tier, mean activities
###   TFSig_Null_Counts.csv / TFSig_Null_Summary.csv     buffered-TF counts, real vs relabelings
###   TFSig_Trajectories.csv               middle -> old fate per TF x cell type
###   TFSig_vs_RankCalls.csv               overlap with the current (rank-based) calls
###   TFSig_Heatmap_<sex>_<age>.pdf/png    TFs buffered in >= 2 of 3 fits, shade = 2/3 or 3/3, bar = cell types
###   TFSig_NullCheck.pdf/png              real number of buffered TFs vs relabelings, per cell type
### Existing figures are skipped; delete a figure to redraw it.
################################################################################

Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({
  library(conflicted); library(DESeq2); library(dplyr); library(tidyr); library(ggplot2); library(parallel); library(decoupleR)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::mutate, .quiet=TRUE); conflicts_prefer(base::intersect, .quiet=TRUE)
  conflicts_prefer(base::setdiff, .quiet=TRUE); conflicts_prefer(base::union, .quiet=TRUE)
  conflicts_prefer(base::unique, .quiet=TRUE); conflicts_prefer(base::Filter, .quiet=TRUE)
  conflicts_prefer(base::setequal, .quiet=TRUE); conflicts_prefer(base::unname, .quiet=TRUE)
})
if (requireNamespace("RhpcBLASctl", quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(1)

################################################################################
### SETTINGS
################################################################################

base_path <- "PATH_TO_NE_OUTPUT"
COLLECTRI_RDS <- "PATH_TO_COLLECTRI/collectri_mouse.rds"
CACHE_DIR  <- file.path(base_path, "HeldOut_FitCache")                       # optional DESeq2 fit cache (used if present)
MASTER     <- file.path(base_path, "Buffering_Candidates_Master")            # for Step1 cell types and the rank-based calls
OUT        <- file.path(MASTER, "TF_Significance")

SIG_P      <- 0.05      # significance threshold for each ULM activity
SIG_MODES  <- c("fdr", "p")   # both are run from the same ULM results, each into its own subfolder of OUT:
                         #   fdr/ = BH FDR < SIG_P (across TFs, within each fit and contrast) - recommended for the main figure
                         #   p/   = ULM p < SIG_P with replication across fits (same logic as the gene tiers) - inclusive version
MIN_REGULON <- 5        # TFs need >= 5 targets among the tested genes (as in the master)
CT_GATE    <- TRUE      # TRUE = only the cell types that pass the master's Step 1 cosine gate (as in Figure 6); FALSE = all
RUN_NULL   <- TRUE      # run the 19 relabelings (slow: ~20x the DESeq2 fits); FALSE = real labeling only
N_CORES    <- 16
HEAT_TOP   <- 15        # TFs per heatmap panel (ranked by cell types in TF-T3, then TF-T1)
RANK_TOP   <- 5         # ranked plot: top buffered and top amplified TFs shown per cell type (significant ones only)
RANK_AGES  <- c("old", "middle")   # ages drawn in the ranked plot
CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN","Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell")
CT_SHORT <- c("DG ExNs"="DG","CA3 ExNs"="CA3","CA1 ExNs"="CA1","CA2 ExNs"="CA2","Sub ExNs"="Sub","IN"="IN","Astrocytes"="Astro",
              "Oligodendrocytes"="Oligo","OPC"="OPC","Microglia"="Micro","Endothelial"="Endo","Pericyte"="Peri","APC"="APC",
              "Stromal"="Stromal","Cajal-Retzius Cell"="CR")

dir.create(file.path(OUT, "units"), recursive=TRUE, showWarnings=FALSE)
wcsv <- function(x, name) write.csv(x, file.path(OUT, name), row.names=FALSE)
save_both <- function(p, stem, w, h) {
  f <- file.path(OUT, stem)
  if (file.exists(paste0(f, ".pdf"))) { cat("  exists, skipped:", basename(f), "\n"); return(invisible()) }
  ggsave(paste0(f, ".pdf"), p, width=w, height=h, limitsize=FALSE)
  ggsave(paste0(f, ".png"), p, width=w, height=h, dpi=300, bg="white", limitsize=FALSE)
  cat("  saved:", basename(f), "\n")
}
banner <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep="")

################################################################################
### DATA, GROUPS, FITS, NETWORK (identical to Buffering_Candidates_Master.R)
################################################################################

banner("Loading")
pb <- readRDS(file.path(base_path, "Buffering_LOO_Split", "pseudobulk_counts_and_meta.RDS"))
pb_all <- pb$pb_all; pb_meta <- pb$pb_meta
age_vals <- unique(pb_meta$age_norm)
AGES <- c(middle = age_vals[!age_vals %in% c("young","old")][1], old = "old")
groups_of <- function(ag, sx, ct) {
  m <- pb_meta[pb_meta$celltype==ct & pb_meta$sex_norm==sx, ]
  list(m=m, ref=m$sample[m$age_norm=="young" & m$condition=="HC"],
       hc=m$sample[m$age_norm==AGES[[ag]] & m$condition=="HC"],
       ne=m$sample[m$age_norm==AGES[[ag]] & m$condition=="EE"])
}
fit_key <- function(num, den) paste(paste(sort(den), collapse="+"), paste(sort(num), collapse="+"), sep=" -> ")
FIT_CACHE <- list()
if (dir.exists(CACHE_DIR)) {
  for (f in list.files(CACHE_DIR, pattern="\\.RDS$", full.names=TRUE)) FIT_CACHE <- c(FIT_CACHE, readRDS(f)$fits)
  cat(sprintf("Fit cache: %d fits\n", length(FIT_CACHE)))
}
fit_one <- function(num, den) {
  hit <- FIT_CACHE[[fit_key(num, den)]]; if (!is.null(hit)) return(hit[, c("gene","lfc","stat","pvalue")])
  cols <- c(den, num)
  cnt  <- round(as.matrix(pb_all[, cols, drop=FALSE]))
  cnt  <- cnt[rowSums(cnt >= 10) >= 2, , drop=FALSE]
  if (nrow(cnt) < 10) return(NULL)
  cd <- data.frame(grp=factor(rep(c("den","num"), c(length(den), length(num))), levels=c("den","num")), row.names=cols)
  tryCatch({
    dds <- suppressMessages(DESeqDataSetFromMatrix(cnt, cd, design=~grp))
    dds <- tryCatch(DESeq(dds, quiet=TRUE), error=function(e) {
      d <- estimateSizeFactors(dds); d <- estimateDispersionsGeneEst(d)
      dispersions(d) <- mcols(d)$dispGeneEst; nbinomWaldTest(d, quiet=TRUE) })
    r <- results(dds, name="grp_num_vs_den")
    data.frame(gene=rownames(r), lfc=r$log2FoldChange, stat=r$stat, pvalue=r$pvalue, stringsAsFactors=FALSE)
  }, error=function(e) NULL)
}
wald <- function(f) if (is.null(f)) NULL else { s <- setNames(f$stat, f$gene); s[!is.na(s)] }
uniprot <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
net_signed <- readRDS(COLLECTRI_RDS) %>%
  dplyr::filter(!grepl(uniprot, source), !grepl(uniprot, target), target %in% rownames(pb_all)) %>%
  dplyr::distinct(source, target, .keep_all=TRUE)
cat(sprintf("CollecTRI: %d TFs, %d interactions\n", dplyr::n_distinct(net_signed$source), nrow(net_signed)))

# ULM on both contrasts of one held-out fit; keeps activity (t) AND p-value
ulm_fit <- function(B, E) {
  g <- intersect(names(B), names(E)); if (length(g) < 500) return(NULL)
  one <- function(v, lab) {                           # one contrast at a time; condition labelled here
    m <- matrix(v[g], ncol=1, dimnames=list(g, lab))
    r <- tryCatch(decoupleR::run_ulm(mat=m, net=net_signed, .source="source", .target="target", .mor="mor",
                                     minsize=MIN_REGULON), error=function(e) NULL)
    if (is.null(r)) return(NULL)
    r <- r[r$statistic == "ulm", c("source", "score", "p_value")]
    r$fdr <- p.adjust(r$p_value, "BH")                # BH across TFs within this fit and contrast
    setNames(r, c("source", paste0(c("a_", "p_", "fdr_"), lab)))
  }
  rb <- one(B, "B"); re <- one(E, "E")
  if (is.null(rb) || is.null(re)) return(NULL)
  dplyr::inner_join(rb, re, by="source")
}

################################################################################
### UNITS (cell types as in Figure 6)
################################################################################

all_units <- tidyr::crossing(sex=unique(pb_meta$sex_norm), celltype=intersect(CT_ORDER, unique(pb_meta$celltype)),
                             age=c("middle","old"))
if (CT_GATE) {
  f1 <- file.path(MASTER, "Step1_CellTypes.csv")
  if (!file.exists(f1)) stop("CT_GATE = TRUE but the master Step 1 table is missing: ", f1, " (set CT_GATE <- FALSE to use all cell types)")
  gate <- read.csv(f1, stringsAsFactors=FALSE)
  units <- all_units %>% dplyr::semi_join(gate %>% dplyr::select(sex, celltype), by=c("sex","celltype"))
} else units <- all_units
cat(sprintf("Units: %d (sex x cell type x age)\n", nrow(units)))

################################################################################
### ULM PER FIT, REAL + RELABELINGS (cached per unit)
################################################################################

banner("ULM per held-out fit (real labeling first, then relabelings)")
run_labeling <- function(hc, ne, ref) {
  dplyr::bind_rows(lapply(seq_along(hc), function(k) {
    a <- ulm_fit(wald(fit_one(hc[-k], ref)), wald(fit_one(ne, hc[k])))
    if (is.null(a)) NULL else dplyr::mutate(a, k=k)
  }))
}
short_id <- function(x) sub(".*--", "", x)
for (i in seq_len(nrow(units))) {
  ag <- units$age[i]; sx <- units$sex[i]; ct <- units$celltype[i]
  f_out <- file.path(OUT, "units", paste0(ag, "__", sx, "__", gsub("[^A-Za-z0-9]+", "_", ct), ".RDS"))
  have <- if (file.exists(f_out)) readRDS(f_out) else NULL
  g <- groups_of(ag, sx, ct)
  if (length(g$ref) < 2 || length(g$hc) < 3 || length(g$ne) < 3) { cat(sprintf("  skip (n): %s %s %s\n", ag, sx, ct)); next }
  pool <- c(g$hc, g$ne); mouse <- setNames(short_id(g$m$orig.ident[match(pool, g$m$sample)]), pool)
  labs <- combn(pool, length(g$ne), simplify=FALSE)
  keys <- vapply(labs, function(l) paste(sort(mouse[l]), collapse="|"), character(1))
  is_real <- vapply(labs, function(l) setequal(l, g$ne), logical(1))
  todo <- which(!keys %in% unique(have$key))                    # resume: only labelings not computed yet
  if (!RUN_NULL) todo <- intersect(todo, which(is_real))
  if (!length(todo)) { cat(sprintf("  done already: %s %s %s (%d labelings)\n", ag, sx, ct, dplyr::n_distinct(have$key))); next }
  cat(sprintf("  %s %s %s: %d labeling(s) ... ", ag, sx, ct, length(todo)))
  one_lab <- function(j) tryCatch({
    lab <- labs[[j]]; s <- run_labeling(setdiff(pool, lab), lab, g$ref)
    if (!nrow(s)) NULL else s %>% dplyr::mutate(age=ag, sex=sx, celltype=ct, real=is_real[j], key=keys[j], .before=1)
  }, error=function(e) NULL)
  res <- mclapply(todo, one_lab, mc.cores=N_CORES)
  failed <- todo[!vapply(res, is.data.frame, logical(1))]
  if (length(failed)) {                                          # parallel workers can drop jobs silently: retry one at a time
    cat(sprintf("retrying %d serially ... ", length(failed)))
    res <- c(res, lapply(failed, one_lab))
  }
  res <- Filter(is.data.frame, res)
  n_done <- dplyr::n_distinct(c(have$key, unlist(lapply(res, function(r) unique(r$key)))))
  if (n_done < ifelse(RUN_NULL, length(labs), 1)) cat(sprintf("WARNING only %d of %d labelings; rerun to fill in ... ", n_done, length(labs)))
  saveRDS(dplyr::bind_rows(have, dplyr::bind_rows(Filter(is.data.frame, res))), f_out); cat("saved\n")
}

################################################################################
### CALLS
################################################################################

per_fit_raw <- dplyr::bind_rows(lapply(list.files(file.path(OUT, "units"), pattern="\\.RDS$", full.names=TRUE), readRDS))
if (!nrow(per_fit_raw)) stop("No ULM results found in ", file.path(OUT, "units"))
OUT_BASE <- OUT
for (SIG_MODE in SIG_MODES) {                      # ---- everything below runs once per threshold, into OUT/<mode>/ ----
OUT <- file.path(OUT_BASE, SIG_MODE); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
banner(paste0("Calls: ", ifelse(SIG_MODE == "fdr", "BH FDR", "p"), " < ", SIG_P, "  ->  ", OUT))
if (exists("null_sum")) rm(null_sum)
per_fit <- per_fit_raw %>%
  dplyr::mutate(pB_use = if (SIG_MODE == "fdr") fdr_B else p_B, pE_use = if (SIG_MODE == "fdr") fdr_E else p_E,
                # intersection-union test: "aging AND stimulation significant" -> p = the larger of the two p-values
                p_joint_buf = ifelse(sign(a_B) != sign(a_E), pmax(pB_use, pE_use), 1),
                p_joint_amp = ifelse(sign(a_B) == sign(a_E), pmax(pB_use, pE_use), 1),
                sigB = if (SIG_MODE == "fdr") fdr_B < SIG_P else p_B < SIG_P,
                sigE = if (SIG_MODE == "fdr") fdr_E < SIG_P else p_E < SIG_P,
                buf_fit = sigB & sigE & sign(a_B) != sign(a_E),
                amp_fit = sigB & sigE & sign(a_B) == sign(a_E))
call_tbl <- function(d) d %>% dplyr::group_by(age, sex, celltype, real, key, source) %>%
  dplyr::summarise(n_fits=dplyr::n(), fits_buffered=sum(buf_fit, na.rm=TRUE), fits_amplified=sum(amp_fit, na.rm=TRUE),
                   aging_act=mean(a_B), stim_act=mean(a_E), max_p_aging=max(p_B), max_p_stim=max(p_E),
                   # median over the 3 fits (a missing fit counts as p = 1): median < SIG_P <=> passes in >= 2 of 3 fits
                   p_buffered=median(c(p_joint_buf, rep(1, 3 - dplyr::n()))),
                   p_amplified=median(c(p_joint_amp, rep(1, 3 - dplyr::n()))), .groups="drop") %>%
  dplyr::mutate(tier = dplyr::case_when(n_fits == 3 & fits_buffered == 3 ~ "TF-T1 buffered (3/3)",
                                        fits_buffered >= 2 ~ "TF-T3 buffered (2/3)",
                                        n_fits == 3 & fits_amplified == 3 ~ "TF-A1 amplified (3/3)",
                                        fits_amplified >= 2 ~ "TF-A3 amplified (2/3)",
                                        TRUE ~ "not called"),
                buffered = fits_buffered >= 2, amplified = fits_amplified >= 2 & !buffered)
calls_all <- call_tbl(per_fit)
calls <- calls_all %>% dplyr::filter(real) %>% dplyr::select(-real, -key)
wcsv(per_fit %>% dplyr::filter(real) %>% dplyr::select(-real, -key) %>% dplyr::arrange(sex, age, celltype, source, k), "TFSig_PerFit.csv")
wcsv(calls %>% dplyr::arrange(sex, age, celltype, dplyr::desc(fits_buffered), dplyr::desc(fits_amplified)), "TFSig_Calls.csv")
print(as.data.frame(calls %>% dplyr::count(sex, age, tier) %>% tidyr::pivot_wider(names_from=tier, values_from=n, values_fill=0)), row.names=FALSE)

# ---- chance check: number of buffered TFs, real vs relabelings ----
null_counts <- calls_all %>% dplyr::group_by(age, sex, celltype, real, key) %>%
  dplyr::summarise(n_buffered=sum(buffered), n_T1=sum(tier == "TF-T1 buffered (3/3)"), n_amplified=sum(amplified), .groups="drop") %>%
  # how many of the real eNE mice a relabeling still calls eNE: relabelings that keep 2 of 3 retain much of a real effect,
  # so the strict null uses only relabelings that keep at most 1 (well scrambled)
  dplyr::group_by(age, sex, celltype) %>%
  dplyr::mutate(kept_real_ne = vapply(strsplit(key, "|", fixed=TRUE), function(v)
                  length(intersect(v, strsplit(key[real][1], "|", fixed=TRUE)[[1]])), integer(1))) %>% dplyr::ungroup()
wcsv(null_counts, "TFSig_Null_Counts.csv")
if (any(!null_counts$real)) {
  null_sum <- null_counts %>% dplyr::group_by(age, sex, celltype) %>%
    dplyr::summarise(real_buffered=n_buffered[real][1], null_mean=mean(n_buffered[!real]), null_max=max(n_buffered[!real]),
                     n_relabelings=sum(!real), emp_p=(1 + sum(n_buffered[!real] >= n_buffered[real][1])) / (1 + sum(!real)),
                     n_scrambled=sum(!real & kept_real_ne <= 1),
                     null_mean_scrambled=mean(n_buffered[!real & kept_real_ne <= 1]),
                     emp_p_scrambled=(1 + sum(n_buffered[!real & kept_real_ne <= 1] >= n_buffered[real][1])) / (1 + sum(!real & kept_real_ne <= 1)),
                     real_amplified=n_amplified[real][1], null_amp_mean=mean(n_amplified[!real]), .groups="drop")
  wcsv(null_sum, "TFSig_Null_Summary.csv")
  cat("\n=== Buffered TFs (>= 2/3 fits), real vs relabelings ===\n"); print(as.data.frame(null_sum), row.names=FALSE)
  # per TF: how often called in relabelings (specificity)
  tf_null <- calls_all %>% dplyr::filter(!real) %>% dplyr::group_by(age, sex, celltype, source) %>%
    dplyr::summarise(relabelings_called=sum(buffered), .groups="drop")
  calls <- calls %>% dplyr::left_join(tf_null, by=c("age","sex","celltype","source")) %>%
    dplyr::mutate(relabelings_called = ifelse(is.na(relabelings_called), 0L, relabelings_called))
  wcsv(calls %>% dplyr::arrange(sex, age, celltype, dplyr::desc(fits_buffered)), "TFSig_Calls.csv")
}

# ---- trajectories middle -> old ----
traj <- calls %>% dplyr::select(age, sex, celltype, source, buffered, amplified) %>%
  tidyr::pivot_wider(names_from=age, values_from=c(buffered, amplified), values_fill=FALSE) %>%
  dplyr::mutate(fate = dplyr::case_when(buffered_middle & buffered_old ~ "Maintained",
                                        buffered_middle & amplified_old ~ "Lost (now amplified)",
                                        buffered_middle ~ "Lost",
                                        buffered_old ~ "Age-emergent",
                                        amplified_old ~ "Amplified, old",
                                        TRUE ~ NA_character_)) %>% dplyr::filter(!is.na(fate))
wcsv(traj, "TFSig_Trajectories.csv")

# ---- comparison with the current rank-based calls (master Step 5) ----
f5 <- file.path(MASTER, "Step5_TFActivity_PerCellType.csv")
if (file.exists(f5)) {
  rk <- read.csv(f5, stringsAsFactors=FALSE) %>% dplyr::transmute(age, sex, celltype, source, rank_buffered=buf %in% TRUE)
  cmp <- dplyr::full_join(calls %>% dplyr::select(age, sex, celltype, source, sig_buffered=buffered), rk,
                          by=c("age","sex","celltype","source")) %>%
    dplyr::mutate(dplyr::across(c(sig_buffered, rank_buffered), ~ .x %in% TRUE)) %>%
    dplyr::filter(sig_buffered | rank_buffered) %>%
    dplyr::mutate(agreement = dplyr::case_when(sig_buffered & rank_buffered ~ "both",
                                               sig_buffered ~ "significance only", TRUE ~ "rank only"))
  wcsv(cmp, "TFSig_vs_RankCalls.csv")
  cat("\n=== Significance-based vs current rank-based buffered calls ===\n")
  print(as.data.frame(cmp %>% dplyr::count(sex, age, agreement) %>% tidyr::pivot_wider(names_from=agreement, values_from=n, values_fill=0)),
        row.names=FALSE)
}

################################################################################
### FIGURES
################################################################################

banner("Figures")
# ---- ranked plot: per cell type, the top buffered (right) and top amplified (left) TFs by their joint p-value ----
rank_tbl <- calls %>% dplyr::filter(buffered | amplified) %>%
  dplyr::mutate(direction = ifelse(buffered, "Buffered", "Amplified"), p = ifelse(buffered, p_buffered, p_amplified),
                value = -log10(pmax(p, 1e-300)) * ifelse(buffered, 1, -1)) %>%
  dplyr::group_by(sex, age, celltype, direction) %>% dplyr::slice_min(p, n=RANK_TOP, with_ties=FALSE) %>% dplyr::ungroup()
wcsv(calls %>% dplyr::filter(buffered | amplified) %>%
       dplyr::mutate(direction=ifelse(buffered, "Buffered", "Amplified"), p=ifelse(buffered, p_buffered, p_amplified)) %>%
       dplyr::arrange(sex, age, celltype, direction, p) %>%
       dplyr::select(sex, age, celltype, direction, source, p, fits_buffered, fits_amplified, aging_act, stim_act),
     "TFSig_Ranked.csv")
for (sx in unique(rank_tbl$sex)) for (ag in intersect(RANK_AGES, unique(rank_tbl$age))) {
  d <- rank_tbl %>% dplyr::filter(sex == sx, age == ag)
  if (!nrow(d)) next
  d <- d %>% dplyr::mutate(celltype = factor(celltype, levels=intersect(CT_ORDER, unique(celltype))),
                           lab = paste(source, celltype, sep="___")) %>%
    dplyr::arrange(celltype, value) %>% dplyr::mutate(lab = factor(lab, levels=unique(lab)))
  p <- ggplot(d, aes(value, lab, fill=direction)) +
    geom_col(width=0.7) + geom_vline(xintercept=0, color="grey40", linewidth=0.3) +
    geom_text(aes(label=source, hjust=ifelse(value > 0, -0.08, 1.08)), size=2.4, fontface="italic") +
    scale_y_discrete(labels=NULL) + scale_x_continuous(expand=expansion(mult=0.35), labels=function(x) abs(x)) +
    scale_fill_manual(values=c(Buffered="#2B6CA3", Amplified="#B2402A"), name=NULL) +
    facet_wrap(~celltype, scales="free_y", ncol=4) + theme_bw(base_size=9) +
    theme(axis.ticks.y=element_blank(), panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
          strip.text=element_text(face="bold"), legend.position="bottom",
          plot.title=element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5, size=7.5, color="grey35")) +
    labs(title=sprintf("Top buffered and amplified TFs, %s %s", ag, sx),
         subtitle=sprintf(paste0("Ranked by significance: -log10 of the median over 3 held-out fits of the larger of the aging and stimulation ULM %s
",
                                 "(both significant, opposite signs = buffered; same sign = amplified). Top %d each, significant only (%s < %g)."),
                          ifelse(SIG_MODE == "fdr", "FDRs", "p-values"), RANK_TOP, ifelse(SIG_MODE == "fdr", "FDR", "p"), SIG_P),
         x=expression(-log[10]~"(joint p)"), y=NULL)
  nct <- dplyr::n_distinct(d$celltype)
  save_both(p, sprintf("TFSig_Ranked_%s_%s", sx, ag), 11, 2.3 * ceiling(nct / 4) + 1.6)
}
for (sx in unique(calls$sex)) for (ag in c("middle","old")) {
  d <- calls %>% dplyr::filter(sex == sx, age == ag, buffered)
  if (!nrow(d)) next
  rk <- d %>% dplyr::group_by(source) %>%
    dplyr::summarise(n_ct=dplyr::n(), n_t1=sum(fits_buffered == 3), .groups="drop") %>%
    dplyr::arrange(dplyr::desc(n_ct), dplyr::desc(n_t1), source) %>% head(HEAT_TOP)
  h <- d %>% dplyr::filter(source %in% rk$source) %>%
    dplyr::mutate(source=factor(source, levels=rev(rk$source)), ct=factor(CT_SHORT[celltype], levels=CT_SHORT[CT_ORDER]),
                  support=ifelse(fits_buffered == 3, "3 of 3 fits (TF-T1)", "2 of 3 fits (TF-T3)"))
  bars <- rk %>% dplyr::mutate(source=factor(source, levels=rev(rk$source)))
  bars$w <- bars$n_ct / max(bars$n_ct) * 1.6
  p <- ggplot(h, aes(ct, source)) +
    geom_tile(aes(fill=support), color="white", linewidth=0.5) +
    scale_fill_manual(values=c("3 of 3 fits (TF-T1)"="#2B4C8C", "2 of 3 fits (TF-T3)"="#9DB3DA"), name="Buffered in") +
    scale_x_discrete(drop=FALSE) +
    geom_tile(data=bars, aes(x=length(CT_ORDER) + 1.5, y=source, width=w), inherit.aes=FALSE,
              fill="grey55", height=0.6) +
    geom_text(data=bars, aes(x=length(CT_ORDER) + 2.6, y=source, label=n_ct), inherit.aes=FALSE, size=2.4, hjust=0) +
    annotate("text", x=length(CT_ORDER) + 1.5, y=nrow(bars) + 0.9, label="Cell types", size=2.4, fontface="bold") +
    coord_cartesian(clip="off", xlim=c(0.5, length(CT_ORDER) + 3)) +
    theme_classic(base_size=9) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, color="black"), axis.line=element_blank(), axis.ticks=element_blank(),
          plot.title=element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5, size=7, color="grey35"),
          legend.position="bottom") +
    labs(title=sprintf("Buffered TFs, %s %s", ag, sx),
         subtitle=sprintf("Aging and stimulation activities each significant (ULM %s < %g) with opposite signs, in >= 2 of 3 held-out fits",
                          ifelse(SIG_MODE == "fdr", "FDR", "p"), SIG_P), x=NULL, y=NULL)
  save_both(p, sprintf("TFSig_Heatmap_%s_%s", sx, ag), 6.5, max(3.5, 0.24 * nrow(bars) + 2.2))
}
if (exists("null_sum")) {
  nd <- null_counts %>% dplyr::mutate(ct=factor(CT_SHORT[celltype], levels=rev(CT_SHORT[CT_ORDER])),
                                      panel=paste(tools::toTitleCase(sex), age))
  p <- ggplot(nd, aes(n_buffered, ct)) +
    geom_point(data=nd %>% dplyr::filter(!real), aes(color=ifelse(kept_real_ne <= 1, "scrambled (keeps <= 1 real eNE mouse)",
                                                                  "partly real (keeps 2 real eNE mice)")),
               size=1.2, position=position_jitter(height=0.15, seed=1)) +
    scale_color_manual(values=c("scrambled (keeps <= 1 real eNE mouse)"="grey45", "partly real (keeps 2 real eNE mice)"="grey80"), name=NULL) +
    geom_point(data=nd %>% dplyr::filter(real), color="#B2402A", size=2.6) +
    facet_wrap(~panel, nrow=1) + theme_bw(base_size=9) +
    theme(strip.text=element_text(face="bold"), plot.title=element_text(face="bold", hjust=0.5), legend.position="bottom",
          plot.subtitle=element_text(hjust=0.5, size=7.5, color="grey35")) +
    labs(title="Number of buffered TFs: real labeling vs relabeled animals",
         subtitle="Red = real; grey = each of the 19 alternative assignments of the six aged animals to HC / eNE (dark = well scrambled)",
         x="TFs buffered (>= 2 of 3 fits)", y=NULL)
  save_both(p, "TFSig_NullCheck", 11, 4.5)
}
cat("\nDone:", OUT, "\n")
}                                                  # ---- end of per-threshold loop ----
OUT <- OUT_BASE
