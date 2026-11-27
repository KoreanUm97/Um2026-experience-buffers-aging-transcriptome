################################################################################
### Buffering_Pipeline.R
###
### NE buffering analysis and the Figure 5 panels, in one file. Run top to bottom:
###     Rscript Buffering_Pipeline.R        (or source it / Run All in RStudio)
###
###   Step 1  robustness   cell-type gate: middle/old held-out cosines        -> Buffering_MiddleAge_Robustness/<date>/
###   Step 2  master       gene tiers + TF activity                           -> Buffering_Candidates_Master/
###   Step 3  methodcheck  TF calls with MLM and DoRothEA (needed by step 4)  -> .../TF_MethodCheck/
###   Step 4  atlas        TF atlas (needed by step 5)                        -> .../TF_Atlas/
###   Step 5  fold         gene fold support (Fig 5c filter)                  -> .../TF_FoldConsistency/
###   Step 6  bi           buffering index from the DESeq2 tables             -> Figures/<date>_Buffering_Index/
###   Step 7  figure5      Figure 5 panels + ED10b, ED10d, ED10f              -> Paper_Figures/Figure5/<date>/
###   Step 8  regulon      regulon heatmaps, top candidate TFs per category   -> Paper_Figures/Figure5/Regulon_Heatmaps/
###                        (not run by default: STOP_AFTER = "figure5")
###
### Before a step runs, the pipeline checks whether its output is already there, complete, and newer than
### the step before it: if yes it is skipped; if it is out of date, the old output folder is renamed "..._superseded_<time>"
### (nothing is deleted) so stale cached pieces are not reused, and the step reruns. So after a crash, just run again.
### A log goes to <NE_BASE>/Pipeline_Logs/.
###
### Inputs
###   <NE_BASE>/Buffering_LOO_Split/pseudobulk_counts_and_meta.RDS   NE pseudobulk counts (mouse x cell type) + metadata
###   <NE_BASE>/DESeq2_Results/<sex>_ALL_<comparison>.csv            NE DESeq2 results
###   COLLECTRI_RDS                                                  CollecTRI mouse network (decoupleR::get_collectri)
###   <DR_BASE>/DREADD_Buffering_HeldOut/, DESeq2_mScarletPos_Approach{B,D}_MergedIN/   DREADD inputs (Fig 5b, 5e)
### Step 7 reuses the DREADD metrics written by Buffering_Paper_Figures.R (SourceData_Fig5_DREADD_CellType_Metrics.csv)
### when present, and otherwise computes them.
################################################################################

suppressPackageStartupMessages(library(conflicted))
conflicts_prefer(base::setdiff, .quiet=TRUE)

NE_BASE       <- "PATH_TO_NE_OUTPUT"
DR_BASE       <- "PATH_TO_DREADD_OUTPUT"
COLLECTRI_RDS <- "PATH_TO_COLLECTRI/collectri_mouse.rds"
DOROTHEA_RDS  <- "PATH_TO_COLLECTRI/dorothea_mouse_ABC.rds"   # written on the first run if missing (step 3)
FORCE   <- character(0)     # steps to rerun even if up to date, e.g. c("bi","figure5")
STOP_AFTER <- "figure5"     # stop after the Figure 5 panels (regulon step not run)
TODAY   <- format(Sys.Date(), "%Y-%m-%d")

MASTER_DIR <- file.path(NE_BASE, "Buffering_Candidates_Master")
ROB_ROOT   <- file.path(NE_BASE, "Buffering_MiddleAge_Robustness")
DESEQ_DIR  <- file.path(NE_BASE, "DESeq2_Results")

LOG_DIR <- file.path(NE_BASE, "Pipeline_Logs"); dir.create(LOG_DIR, recursive=TRUE, showWarnings=FALSE)
LOG <- file.path(LOG_DIR, paste0("Buffering_Pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
say <- function(...) { m <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", ...); cat("\n####", m, "\n"); cat(m, "\n", file=LOG, append=TRUE) }
newest <- function(root, pattern) { f <- list.files(root, pattern=pattern, recursive=TRUE, full.names=TRUE)
  if (!length(f)) NA_character_ else f[which.max(file.mtime(f))] }
mt <- function(f) if (is.na(f) || !file.exists(f)) -Inf else as.numeric(file.mtime(f))
supersede <- function(d) if (dir.exists(d)) {
  to <- paste0(d, "_superseded_", format(Sys.time(), "%Y%m%d_%H%M%S")); file.rename(d, to); say("  moved old output aside: ", to) }

for (f in c("male_ALL_MEEMHC.csv", "female_ALL_MEEMHC.csv", "male_ALL_OEEYHC.csv", "male_ALL_OHCYHC.csv", "male_ALL_MEEYHC.csv", "male_ALL_MHCYHC.csv"))
  if (!file.exists(file.path(DESEQ_DIR, f))) stop("Missing ", file.path(DESEQ_DIR, f), ": run the DESeq2 pipeline first")
say("Buffering pipeline started. Log: ", LOG)

################################################################################
### Checks: is each step's output present, complete and newer than its input?
################################################################################

f_rob2 <- function() newest(ROB_ROOT, "_Check2_HeldOut_PerMouse\\.csv$")
f_mast <- function() file.path(MASTER_DIR, "Step5_TFActivity_AcrossCellTypes.csv")
f_mc   <- function() file.path(MASTER_DIR, "TF_MethodCheck", c("MethodCheck_AllTFs_PerCellType_mlm_collectri.csv",
                                                           "MethodCheck_AllTFs_PerCellType_ulm_dorothea.csv"))
f_atl  <- function() file.path(MASTER_DIR, "TF_Atlas", "TF_Atlas_Long.csv")
f_fold <- function() file.path(MASTER_DIR, "TF_FoldConsistency", "Gene_FoldSupport_T3A3_v2.csv")
f_bi   <- function() newest(file.path(NE_BASE, "Figures"), "_BI_Summary\\.csv$")

CHECKS <- list(
  robustness = function() {
    f2 <- f_rob2(); f1 <- newest(ROB_ROOT, "_Check1_BiasCorrected\\.csv$")
    if (is.na(f2) || is.na(f1)) return("no Check 1 / Check 2 output")
    c1 <- read.csv(f1); c2 <- read.csv(f2)
    exp <- unique(paste(tolower(c1$age), c1$sex, c1$celltype)); got <- table(paste(tolower(c2$age), c2$sex, c2$celltype))
    miss <- setdiff(exp, names(got)[got >= 3])
    if (length(miss)) return(paste("Check 2 incomplete (crashed workers?):", paste(head(miss, 8), collapse="; ")))
    "ok" },
  master = function() {
    need <- file.path(MASTER_DIR, c("Step1_CellTypes.csv", "Step2_GeneCalls.RDS", "Step5_TFActivity_PerCellType.csv",
                                   "Step5_TFActivity_AcrossCellTypes.csv"))
    if (!all(file.exists(need))) return(paste("missing:", paste(basename(need[!file.exists(need)]), collapse=", ")))
    if (mt(need[1]) < mt(f_rob2())) return("out of date (older than the robustness output)")
    "ok" },
  methodcheck = function() {
    f <- f_mc(); if (!all(file.exists(f))) return("MethodCheck tables missing")
    if (min(sapply(f, mt)) < mt(f_mast())) return("out of date (older than master v2)")
    "ok" },
  atlas = function() {
    if (!file.exists(f_atl())) return("TF_Atlas_Long.csv missing")
    if (mt(f_atl()) < max(sapply(f_mc(), mt))) return("out of date (older than MethodCheck)")
    "ok" },
  fold = function() {
    if (!file.exists(f_fold())) return("Gene_FoldSupport_T3A3_v2.csv missing")
    if (mt(f_fold()) < mt(f_atl())) return("out of date (older than the TF atlas)")
    "ok" },
  bi = function() {
    f <- f_bi(); if (is.na(f)) return("no BI summary")
    if (mt(f) < mt(file.path(DESEQ_DIR, "male_ALL_MEEMHC.csv"))) return("out of date (older than DESeq2_Results)")
    "ok" })
# folder each step owns: renamed aside before an out-of-date step reruns, so its own caches are not reused
OWNS <- list(robustness = NULL,                                     # dated folder per run; nothing cached across runs
             master      = MASTER_DIR,                              # (downstream folders live inside it and rerun too)
             methodcheck = file.path(MASTER_DIR, "TF_MethodCheck"),
             atlas       = file.path(MASTER_DIR, "TF_Atlas"),
             fold        = file.path(MASTER_DIR, "TF_FoldConsistency"),
             bi          = NULL)

################################################################################
### Step 6: buffering index from the DESeq2 tables
###   BI = sum(log2FC[aged NE vs young HC]^2) / sum(log2FC[aged HC vs young HC]^2), per sex x cell type,
###   over genes with a log2 fold change in both tables:   old = OEEYHC / OHCYHC,  middle = MEEYHC / MHCYHC.
################################################################################

step_bi <- function() {
  suppressPackageStartupMessages(library(dplyr))
  BI_COMPS <- list(old = c(ne="OEEYHC", hc="OHCYHC"), middle = c(ne="MEEYHC", hc="MHCYHC"))
  OUT_BI   <- file.path(NE_BASE, "Figures", paste0(TODAY, "_Buffering_Index"))
  dir.create(OUT_BI, recursive=TRUE, showWarnings=FALSE)
  rd <- function(dir, sx, cp) { f <- file.path(dir, paste0(sx, "_ALL_", cp, ".csv")); if (!file.exists(f)) return(NULL)
    read.csv(f, stringsAsFactors=FALSE) %>% dplyr::select(celltype, gene, lfc=log2FoldChange) }
  bi_from <- function(dir) dplyr::bind_rows(lapply(c("female", "male"), function(sx) dplyr::bind_rows(lapply(names(BI_COMPS), function(ag) {
    n <- rd(dir, sx, BI_COMPS[[ag]][["ne"]]); h <- rd(dir, sx, BI_COMPS[[ag]][["hc"]]); if (is.null(n) || is.null(h)) return(NULL)
    dplyr::inner_join(n, h, by=c("celltype", "gene"), suffix=c("_ne", "_hc")) %>% dplyr::filter(!is.na(lfc_ne), !is.na(lfc_hc)) %>%
      dplyr::group_by(celltype) %>% dplyr::summarise(n_genes=dplyr::n(), ss_ne=sum(lfc_ne^2), ss_hc=sum(lfc_hc^2), .groups="drop") %>%
      dplyr::mutate(bi=ss_ne / ss_hc, label=tools::toTitleCase(ag), sex=sx, comparison=paste0(BI_COMPS[[ag]][["ne"]], " / ", BI_COMPS[[ag]][["hc"]]))
  })))) %>% dplyr::select(label, sex, celltype, bi, n_genes, ss_ne, ss_hc, comparison)
  bi <- bi_from(DESEQ_DIR)
  if (!nrow(bi)) stop("No BI computed: DESeq2 tables not found in ", DESEQ_DIR)
  for (sx in unique(bi$sex)) write.csv(bi %>% dplyr::filter(sex == sx),
                                       file.path(OUT_BI, paste0(sx, "_BI_Summary.csv")), row.names=FALSE)
  print(as.data.frame(bi %>% dplyr::mutate(bi=round(bi, 3)) %>% dplyr::select(label, sex, celltype, bi, n_genes)), row.names=FALSE)
  cat("BI written to", OUT_BI, "\n")
}

################################################################################
### The analysis steps, one function each
################################################################################

STEPS <- list()

################################################################################
### STEP: robustness  (= Buffering_MiddleAge_Robustness.R)
################################################################################

STEPS$robustness <- function() {
################################################################################
### Buffering_MiddleAge_Robustness.R
###
### The same robustness checks used for old age, run for MIDDLE age (and old age
### alongside for comparison). NE dataset, each sex separately (never pooled).
###
###   CHECK 1  Bias-corrected cosine, all mice (cos_dc)
###            cos_dc = (B.E + sum Var(mean XHC)) / (|B| |E|)
###            0 when nothing is real; reported with a simulated no-effect band
###   CHECK 2  Held-out cosine (DESeq2) per held-out XHC mouse, labelled with the
###            mouse ID and sequencing run (from orig.ident), plus
###            leave-one-NE-out: held-out cosine recomputed without each NE mouse
###   CHECK 3  Per-mouse position on the aging axis (axis built without that
###            mouse) + global relabeling test across cell types
###            (per sex: 20 labelings; joint male x female: 400, min p = 0.0025)
###
###   B = XHC vs YHC (aging)    E = XNE vs XHC (NE)    X = Middle or Old
###
### Input: Buffering_LOO_Split/pseudobulk_counts_and_meta.RDS (NE pseudobulk cache)
### Runtime: checks 1 and 3 minutes; check 2 ~15-40 min with N_CORES cores.
################################################################################

suppressPackageStartupMessages({
  library(conflicted); library(DESeq2); library(dplyr); library(tidyr)
  library(ggplot2); library(parallel)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::mutate, .quiet=TRUE); conflicts_prefer(base::intersect, .quiet=TRUE)
  conflicts_prefer(base::setdiff, .quiet=TRUE)
})

fig_date <- format(Sys.Date(),"%Y-%m-%d")
OUT_DIR  <- file.path(ROB_ROOT, fig_date)
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

N_CORES          <- 8
RUN_CHECK2       <- TRUE          # DESeq2-heavy; set FALSE to run checks 1 and 3 only
CELLTYPES_TO_RUN <- NULL          # NULL = all reportable; e.g. c("DG ExNs","CA3 ExNs")
MIN_GENES        <- 1000
N_SIM            <- 50
set.seed(20261019)

CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN",
              "Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell")
CT_COLS <- c("DG ExNs"="#E53935","CA1 ExNs"="#CC0000","CA2 ExNs"="#8E24AA",
             "CA3 ExNs"="#5E35B1","Sub ExNs"="#80D8FF","IN"="#00897B",
             "Astrocytes"="#43A047","Oligodendrocytes"="#1E88E5","OPC"="#00ACC1",
             "Microglia"="#7CB342","Endothelial"="#FB8C00","Pericyte"="#6D4C41",
             "APC"="#546E7A","Stromal"="#27408B","Cajal-Retzius Cell"="#FDD835")

ne <- readRDS(file.path(NE_BASE, "Buffering_LOO_Split", "pseudobulk_counts_and_meta.RDS"))
age_vals <- unique(ne$pb_meta$age_norm)
MIDDLE   <- age_vals[!age_vals %in% c("young","old")][1]
AGES     <- c(Middle = MIDDLE, Old = "old")
cts      <- intersect(CT_ORDER, unique(ne$pb_meta$celltype))
if (!is.null(CELLTYPES_TO_RUN)) cts <- intersect(CELLTYPES_TO_RUN, cts)
short_id <- function(x) sub(".*--", "", x)
run_id   <- function(x) regmatches(x, regexpr("SR0[0-9]+", x))

groups <- function(ct, sx, ag) {
  m <- ne$pb_meta[ne$pb_meta$celltype==ct & ne$pb_meta$sex_norm==sx, ]
  list(m = m,
       y = m$sample[m$age_norm=="young" & m$condition=="HC"],
       o = m$sample[m$age_norm==ag      & m$condition=="HC"],
       e = m$sample[m$age_norm==ag      & m$condition=="EE"])
}

################################################################################
### SHARED HELPERS
################################################################################

log_norm <- function(cnt) {
  cnt <- round(as.matrix(cnt)); cnt <- cnt[rowSums(cnt >= 10) >= 2, , drop=FALSE]
  sf  <- DESeq2::estimateSizeFactorsForMatrix(cnt)
  log2(sweep(cnt, 2, sf, "/") + 1)
}

ne_fit <- function(num, den) {                     # DESeq2 LFC, as in the main pipeline
  cols <- c(den, num)
  cnt  <- round(as.matrix(ne$pb_all[, cols, drop=FALSE]))
  cnt  <- cnt[rowSums(cnt >= 10) >= 2, , drop=FALSE]
  if (nrow(cnt) < 10) return(NULL)
  cd <- data.frame(grp=factor(rep(c("den","num"), c(length(den), length(num))),
                              levels=c("den","num")), row.names=cols)
  tryCatch({
    dds <- suppressMessages(DESeqDataSetFromMatrix(cnt, cd, design=~grp))
    dds <- tryCatch(DESeq(dds, quiet=TRUE), error=function(e) {
      d <- estimateSizeFactors(dds); d <- estimateDispersionsGeneEst(d)
      dispersions(d) <- mcols(d)$dispGeneEst; nbinomWaldTest(d, quiet=TRUE)
    })
    r <- results(dds, name="grp_num_vs_den"); setNames(r$log2FoldChange, rownames(r))
  }, error=function(e) NULL)
}
cosv <- function(A, B) {
  if (is.null(A) || is.null(B)) return(NA_real_)
  g <- intersect(names(A)[!is.na(A)], names(B)[!is.na(B)])
  if (length(g) < 20) return(NA_real_)
  sum(A[g]*B[g]) / (sqrt(sum(A[g]^2)) * sqrt(sum(B[g]^2)))
}

################################################################################
### CHECK 1: bias-corrected cosine (all mice)
################################################################################

dc_stats <- function(X, y, o, e) {
  mY <- rowMeans(X[, y, drop=FALSE]); mO <- rowMeans(X[, o, drop=FALSE]); mE <- rowMeans(X[, e, drop=FALSE])
  vO <- apply(X[, o, drop=FALSE], 1, var) / length(o)
  B <- mO - mY; E <- mE - mO
  c(cos_raw = sum(B*E) / sqrt(sum(B^2) * sum(E^2)),
    cos_dc  = (sum(B*E) + sum(vO)) / sqrt(sum(B^2) * sum(E^2)))
}
dc_null <- function(X, y, o, e) {
  s  <- sqrt(rowMeans(cbind(apply(X[, y], 1, var), apply(X[, o], 1, var), apply(X[, e], 1, var))))
  mu <- rowMeans(X)
  v  <- replicate(N_SIM, {
    Z <- matrix(rnorm(length(mu)*ncol(X), mu, s), nrow=length(mu), dimnames=dimnames(X))
    dc_stats(Z, y, o, e)[["cos_dc"]] })
  c(null_mean = mean(v), null_sd = sd(v))
}

cat("=== CHECK 1: bias-corrected cosine ===\n")
c1 <- list()
for (ag in names(AGES)) for (sx in c("male","female")) for (ct in cts) {
  g <- groups(ct, sx, AGES[[ag]])
  if (min(length(g$y), length(g$o), length(g$e)) < 2) next
  X <- log_norm(ne$pb_all[, c(g$y, g$o, g$e), drop=FALSE]); if (nrow(X) < MIN_GENES) next
  c1[[length(c1)+1]] <- data.frame(age=ag, sex=sx, celltype=ct, n_genes=nrow(X),
                                   t(dc_stats(X, g$y, g$o, g$e)), t(dc_null(X, g$y, g$o, g$e)))
}
c1 <- dplyr::bind_rows(c1) %>%
  dplyr::mutate(beyond_null = cos_dc < null_mean - 2*null_sd)
write.csv(c1, file.path(OUT_DIR, paste0(fig_date, "_Check1_BiasCorrected.csv")), row.names=FALSE)
print(as.data.frame(c1 %>% dplyr::select(age, sex, celltype, cos_raw, cos_dc, null_sd, beyond_null) %>%
                    dplyr::mutate(dplyr::across(where(is.double), ~round(.x, 3)))))

################################################################################
### CHECK 2: held-out per XHC mouse (with run) + leave-one-NE-out (DESeq2)
################################################################################

c2_fit <- NULL; c2_loo <- NULL
if (RUN_CHECK2) {
  cat("\n=== CHECK 2: held-out per mouse + leave-one-NE-out (DESeq2) ===\n")
  jobs <- expand.grid(age=names(AGES), sex=c("male","female"), celltype=cts, stringsAsFactors=FALSE)
  one_job <- function(i) {
    ag <- jobs$age[i]; sx <- jobs$sex[i]; ct <- jobs$celltype[i]
    g  <- groups(ct, sx, AGES[[ag]])
    if (length(g$y) < 2 || length(g$o) < 3 || length(g$e) < 3) return(NULL)
    Bs <- lapply(seq_along(g$o), function(k) ne_fit(g$o[-k], g$y))
    per_k <- function(ne_set) sapply(seq_along(g$o), function(k) cosv(Bs[[k]], ne_fit(ne_set, g$o[k])))
    all_k <- per_k(g$e)
    oid   <- g$m$orig.ident[match(g$o, g$m$sample)]
    fits  <- data.frame(age=ag, sex=sx, celltype=ct, held_out_HC=short_id(oid),
                        run=run_id(oid), cosine=all_k)
    eid   <- g$m$orig.ident[match(g$e, g$m$sample)]
    loo   <- data.frame(age=ag, sex=sx, celltype=ct,
                        set=c("all NE", paste("without", short_id(eid))),
                        cosine=c(mean(all_k, na.rm=TRUE),
                                 sapply(seq_along(g$e), function(j) mean(per_k(g$e[-j]), na.rm=TRUE))))
    message(sprintf("  done: %s %s %s", ag, sx, ct))
    list(fits=fits, loo=loo)
  }
  out    <- mclapply(seq_len(nrow(jobs)), one_job, mc.cores=N_CORES, mc.preschedule=FALSE)   # a crash now loses one job, not a batch
  bad    <- which(sapply(out, function(x) is.null(x) || inherits(x, "try-error")))
  if (length(bad)) { message("Rerunning ", length(bad), " jobs one at a time"); out[bad] <- lapply(bad, one_job) }
  out    <- Filter(function(x) is.list(x) && !inherits(x, "try-error"), out)
  c2_fit <- dplyr::bind_rows(lapply(out, `[[`, "fits"))
  c2_loo <- dplyr::bind_rows(lapply(out, `[[`, "loo"))
  write.csv(c2_fit, file.path(OUT_DIR, paste0(fig_date, "_Check2_HeldOut_PerMouse.csv")), row.names=FALSE)
  write.csv(c2_loo, file.path(OUT_DIR, paste0(fig_date, "_Check2_LeaveOneNEOut.csv")), row.names=FALSE)

  cat("\n-- Held-out cosine by held-out XHC mouse (with sequencing run) --\n")
  for (ag in names(AGES)) for (sx in c("male","female")) {
    tab <- c2_fit %>% dplyr::filter(age==ag, sex==sx) %>%
      dplyr::mutate(col = paste0(held_out_HC, " [", run, "]"), cosine=round(cosine, 3)) %>%
      dplyr::select(celltype, col, cosine) %>%
      tidyr::pivot_wider(names_from=col, values_from=cosine)
    cat(sprintf("\n%s %s\n", ag, sx)); print(as.data.frame(tab))
  }
  cat("\n-- Leave-one-NE-out held-out cosine --\n")
  for (ag in names(AGES)) for (sx in c("male","female")) {
    tab <- c2_loo %>% dplyr::filter(age==ag, sex==sx) %>%
      dplyr::mutate(cosine=round(cosine, 3)) %>%
      tidyr::pivot_wider(names_from=set, values_from=cosine)
    cat(sprintf("\n%s %s\n", ag, sx)); print(as.data.frame(tab))
  }
}

################################################################################
### CHECK 3: per-mouse position + global relabeling test
################################################################################

cat("\n=== CHECK 3: per-mouse position + global relabeling test ===\n")
positions <- function(X, y, old, trt) {
  ctrl <- setdiff(old, trt); ym <- rowMeans(X[, y, drop=FALSE])
  p1 <- function(i, refs) { u <- rowMeans(X[, refs, drop=FALSE]) - ym; sum((X[, i] - ym) * u) / sum(u^2) }
  sapply(old, function(i) if (i %in% ctrl) p1(i, setdiff(ctrl, i))
                          else mean(sapply(seq_along(ctrl), function(k) p1(i, ctrl[-k]))))
}
effect <- function(p, trt) mean(p[trt]) - mean(p[setdiff(names(p), trt)])

pm <- list(); glob <- list()
for (ag in names(AGES)) {
  per_sex <- list()
  for (sx in c("male","female")) {
    effs <- list(); real_key <- NULL
    for (ct in cts) {
      g <- groups(ct, sx, AGES[[ag]])
      if (length(g$y) < 2 || length(g$o) < 3 || length(g$e) < 3) next
      X <- log_norm(ne$pb_all[, c(g$y, g$o, g$e), drop=FALSE]); if (nrow(X) < MIN_GENES) next
      old  <- c(g$o, g$e)
      pos  <- positions(X, g$y, old, g$e)
      pm[[length(pm)+1]] <- data.frame(age=ag, sex=sx, celltype=ct,
        mouse=short_id(g$m$orig.ident[match(old, g$m$sample)]),
        run=run_id(g$m$orig.ident[match(old, g$m$sample)]),
        group=ifelse(old %in% g$e, "NE", "HC"), position=pos)
      labs <- combn(old, 3, simplify=FALSE)
      key  <- sapply(labs, function(t) paste(sort(g$m$orig.ident[match(t, g$m$sample)]), collapse="|"))
      effs[[ct]] <- setNames(sapply(labs, function(t) effect(positions(X, g$y, old, t), t)), key)
      real_key   <- paste(sort(g$m$orig.ident[match(g$e, g$m$sample)]), collapse="|")
    }
    common <- Reduce(intersect, lapply(effs, names))
    gl     <- rowMeans(sapply(effs, function(v) v[common]))
    per_sex[[sx]] <- list(g=gl, real=real_key)
    glob[[length(glob)+1]] <- data.frame(age=ag, sex=sx, n_celltypes=length(effs),
      global_effect=gl[real_key], rank=sum(gl <= gl[real_key]), n_labelings=length(gl),
      p_perm=mean(gl <= gl[real_key]))
  }
  tot <- outer(per_sex$male$g, per_sex$female$g, "+")
  obs <- per_sex$male$g[per_sex$male$real] + per_sex$female$g[per_sex$female$real]
  glob[[length(glob)+1]] <- data.frame(age=ag, sex="joint (M+F)", n_celltypes=NA,
    global_effect=obs, rank=sum(tot <= obs), n_labelings=length(tot), p_perm=mean(tot <= obs))
}
pm   <- dplyr::bind_rows(pm); glob <- dplyr::bind_rows(glob)
write.csv(pm,   file.path(OUT_DIR, paste0(fig_date, "_Check3_PerMouse_Positions.csv")), row.names=FALSE)
write.csv(glob, file.path(OUT_DIR, paste0(fig_date, "_Check3_GlobalTest.csv")), row.names=FALSE)

cat("\n-- Median position per mouse across cell types (0 = young, 1 = old) --\n")
print(as.data.frame(pm %>% dplyr::group_by(age, sex, group, mouse, run) %>%
  dplyr::summarise(median_position=round(median(position), 2), .groups="drop") %>%
  dplyr::arrange(age, sex, group, mouse)))
cat("\n-- Global relabeling test (effect < 0 = NE closer to young; low rank = buffering) --\n")
print(as.data.frame(glob %>% dplyr::mutate(dplyr::across(where(is.double), ~signif(.x, 3)))))

################################################################################
### FIGURE: check 1 (corrected cosine) + check 2 per-mouse fits, Middle vs Old
################################################################################

p1 <- ggplot(c1 %>% dplyr::mutate(celltype=factor(celltype, levels=rev(CT_ORDER)),
                                  sex=factor(sex, levels=c("male","female")),
                                  age=factor(age, levels=c("Middle","Old"))),
             aes(y=celltype)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey30") +
  geom_errorbar(aes(xmin=null_mean-2*null_sd, xmax=null_mean+2*null_sd),
                width=0, color="grey85", linewidth=3, orientation="y") +
  geom_point(aes(x=cos_raw), shape=5, size=2.2, color="grey45") +
  geom_point(aes(x=cos_dc, fill=celltype), shape=21, size=3.2, color="black", stroke=0.4) +
  scale_fill_manual(values=CT_COLS, guide="none") +
  facet_grid(sex ~ age) + theme_classic(base_size=10) +
  labs(title="Check 1: bias-corrected cosine (all mice)",
       subtitle="Open diamond = raw coupled | filled = corrected | grey band = no-effect null (+/- 2 SD)",
       x="Cosine (negative = opposes aging)", y=NULL)
ggsave(file.path(OUT_DIR, paste0(fig_date, "_Check1_BiasCorrected.pdf")), p1, width=10, height=9)
ggsave(file.path(OUT_DIR, paste0(fig_date, "_Check1_BiasCorrected.png")), p1, width=10, height=9, dpi=600)

if (!is.null(c2_fit)) {
  p2 <- ggplot(c2_fit %>% dplyr::mutate(celltype=factor(celltype, levels=rev(CT_ORDER)),
                                        sex=factor(sex, levels=c("male","female")),
                                        age=factor(age, levels=c("Middle","Old"))),
               aes(x=cosine, y=celltype, color=paste0(held_out_HC, " [", run, "]"))) +
    geom_vline(xintercept=0, linetype="dashed", color="grey30") +
    geom_point(size=2.4, alpha=0.9) +
    facet_grid(sex ~ age) + theme_classic(base_size=10) + theme(legend.position="right") +
    labs(title="Check 2: held-out cosine for each held-out XHC mouse",
         subtitle="Colour = held-out HC mouse [sequencing run]",
         x="Held-out cosine", y=NULL, color="Held-out HC mouse")
  ggsave(file.path(OUT_DIR, paste0(fig_date, "_Check2_HeldOut_PerMouse.pdf")), p2, width=12, height=9)
  ggsave(file.path(OUT_DIR, paste0(fig_date, "_Check2_HeldOut_PerMouse.png")), p2, width=12, height=9, dpi=600)
}
cat("\n=== DONE ===\nOutput:", OUT_DIR, "\n")

invisible(NULL)
}

################################################################################
### STEP: master  (= Buffering_Candidates_Master.R)
################################################################################

STEPS$master <- function() {
################################################################################
### Buffering_Candidates_Master.R
###
### One pipeline, run top to bottom, resumable. NE dataset; middle and old age;
### each sex separately (never pooled). Middle age: two of three HC per sex were
### sequenced in a separate run.
###
###   STEP 1  Cell types      gate per age (held-out cosine < 0 AND bias-corrected beyond null);
###                           union per sex = buffers at middle OR old age
###   STEP 2  Gene calls      3 held-out fits per unit; buffered (T1/T3) and amplified (A1/A3)
###   STEP 3  Gene sets +     buffered middle / old, lost (by reason), maintained, age-emergent,
###           CollecTRI       amplified middle / old; Fisher per cell type + Mantel-Haenszel pooled;
###                           run for each tier in TIERS_RUN (T3 = discovery, T1 = confirmation)
###   STEP 4  Timing check    per-mouse circadian genes (is the clock signal collection time?)
###   STEP 5  TF activity     decoupleR ULM on held-out Wald stats, signed CollecTRI;
###                           S = mean_k(-a_aging * a_NE); real + all relabelings of the aged mice;
###                           union cell types at both ages (slow; resumable per unit)
###   STEP 6  Trajectories    per cell type and TF level: maintained / lost (flipped) / gained /
###                           amplified, middle -> old
###   STEP 7  Candidates      targets enriched in a set (STEP 3) AND the TF's own activity moves
###                           the matching way (STEP 5)
###
### Gene tiers (3 held-out fits; per-gene DESeq2 Wald p, unadjusted):
###   T1  aging p<0.05 in 3/3, NE opposite in 3/3, NE p<0.05 in >=2/3
###   T3  aging p<0.05 in >=2/3, NE opposite in >=2/3
###   A1/A3 = same rules with NE in the SAME direction as aging (amplified)
### TF calls (candidate mode): buffered = S > 0 & rank <= ACT_RANK of its own labelings;
###   amplified = S < 0 & rank <= ACT_RANK.
### Held-out fits: B_k = XHC(-k) vs YHC,  E_k = XNE vs XHC(k).
################################################################################

Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({
  library(conflicted); library(DESeq2); library(dplyr); library(tidyr); library(ggplot2)
  library(parallel); library(decoupleR)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::mutate, .quiet=TRUE); conflicts_prefer(base::intersect, .quiet=TRUE)
  conflicts_prefer(base::setdiff, .quiet=TRUE); conflicts_prefer(base::union, .quiet=TRUE)
  conflicts_prefer(base::unique, .quiet=TRUE); conflicts_prefer(base::Filter, .quiet=TRUE)
  conflicts_prefer(base::setequal, .quiet=TRUE); conflicts_prefer(base::unname, .quiet=TRUE)
})
if (requireNamespace("RhpcBLASctl", quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(1)

################################################################################
### PARAMETERS
################################################################################

base_path <- NE_BASE
CACHE_DIR <- file.path(base_path, "HeldOut_FitCache")       # optional, used if present

STEPS    <- c(genes=TRUE, timing=TRUE, tf_activity=TRUE, trajectory=TRUE, candidates=TRUE)
N_CORES  <- 16
TIERS_RUN <- c("T3","T1")   # STEP 3 tiers (amplified tiers follow: A3, A1)
P_AGING  <- 0.05; P_NE <- 0.05
MIN_REGULON <- 5; MIN_SET <- 10; MIN_OVERLAP <- 3
ACT_RANK <- 5               # TF call: rank <= ACT_RANK of its own labelings (candidate mode)
BG_MODES <- c("aging","tested")   # STEP 3 backgrounds, each run separately. "aging": genes that COULD be in the set --
                            #   buffered/amplified/age-emergent: genes with an aging change at that age
                            #   (same aging rule as the tier); lost/maintained: genes buffered at middle age.
                            #   "tested": all genes tested (inflates large regulons whose targets change with age)
RUN_GO   <- FALSE
WATCH <- c("Creb1","Egr1","Egr3","Jun","Fos","Atf6","Xbp1","Ddit3","Atf4","Hsf1","Bmal1","Clock",
           "Foxo1","Foxo3","Sirt1","Nr3c1","Rela","Relb","Nfkb1","Nfkb2","Smad1","Smad4")
CLOCK_GENES <- c("Dbp","Nr1d1","Nr1d2","Per1","Per2","Per3","Tef","Cry1","Npas2","Arntl","Bmal1")

out_dir <- MASTER_DIR
dir.create(file.path(out_dir, "tf_units"), recursive=TRUE, showWarnings=FALSE)
CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN",
              "Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell")
SETS <- c("Buffered, middle","Buffered, old","Lost","Lost: NE no longer opposes",
          "Lost: aging direction changed","Lost: still opposed, weaker","Maintained","Age-emergent",
          "Amplified, middle","Amplified, old")
MAIN_SETS <- c("Buffered, middle","Buffered, old","Lost","Maintained","Age-emergent","Amplified, middle","Amplified, old")
SCOL <- c(`Buffered, middle`="#81C784", `Buffered, old`="#1B5E20", Lost="#C62828", Maintained="#2E7D32",
          `Age-emergent`="#1565C0", `Amplified, middle`="#CE93D8", `Amplified, old`="#6A1B9A")
th <- theme_classic(base_size=8) +
  theme(strip.background=element_rect(fill="grey92", color=NA), strip.text=element_text(face="bold"),
        plot.title=element_text(face="bold", hjust=0.5, size=9), plot.subtitle=element_text(color="grey40", hjust=0.5, size=7.5))
sv <- function(p, dir, stem, w, h) {
  ggsave(file.path(dir, paste0(stem, ".pdf")), p, width=min(w, 40), height=min(h, 40), limitsize=FALSE)
  ggsave(file.path(dir, paste0(stem, ".png")), p, width=min(w, 40), height=min(h, 40), dpi=200, bg="white", limitsize=FALSE)
}
wcsv <- function(x, name, dir=out_dir) write.csv(x, file.path(dir, name), row.names=FALSE)
banner <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep="")

################################################################################
### STEP 1  CELL TYPES
################################################################################

banner("STEP 1  Cell types")
pb <- readRDS(file.path(base_path, "Buffering_LOO_Split", "pseudobulk_counts_and_meta.RDS"))
pb_all <- pb$pb_all; pb_meta <- pb$pb_meta
age_vals <- unique(pb_meta$age_norm)
AGES <- c(middle = age_vals[!age_vals %in% c("young","old")][1], old = "old")
short_id <- function(x) sub(".*--", "", x)

rob_dir <- ROB_ROOT
latest  <- function(pat) sort(list.files(rob_dir, pattern=pat, recursive=TRUE, full.names=TRUE), decreasing=TRUE)[1]
c1  <- read.csv(latest("_Check1_BiasCorrected\\.csv$"), stringsAsFactors=FALSE)
loo <- read.csv(latest("_Check2_LeaveOneNEOut\\.csv$"), stringsAsFactors=FALSE)
ct_status <- dplyr::full_join(
  c1  %>% dplyr::transmute(age=tolower(age), sex, celltype, dc_beyond_null=beyond_null),
  loo %>% dplyr::filter(set=="all NE") %>% dplyr::transmute(age=tolower(age), sex, celltype, cos_heldout=cosine),
  by=c("age","sex","celltype")) %>%
  dplyr::filter(age %in% c("middle","old")) %>%
  dplyr::mutate(pass = !is.na(cos_heldout) & cos_heldout < 0 & dc_beyond_null %in% TRUE) %>%
  dplyr::select(age, sex, celltype, pass) %>%
  tidyr::pivot_wider(names_from=age, values_from=pass, names_prefix="pass_") %>%
  dplyr::mutate(pass_middle = pass_middle %in% TRUE, pass_old = pass_old %in% TRUE,
                ct_status = dplyr::case_when(pass_middle & pass_old ~ "buffers at both ages",
                                             pass_middle ~ "lost at old age", pass_old ~ "gained at old age",
                                             TRUE ~ NA_character_)) %>%
  dplyr::filter(!is.na(ct_status))
units <- tidyr::crossing(ct_status %>% dplyr::select(sex, celltype), age=c("middle","old"))
print(as.data.frame(ct_status %>% dplyr::arrange(sex, ct_status, factor(celltype, levels=CT_ORDER))), row.names=FALSE)
wcsv(ct_status, "Step1_CellTypes.csv")

groups_of <- function(ag, sx, ct) {
  m <- pb_meta[pb_meta$celltype==ct & pb_meta$sex_norm==sx, ]
  list(m=m, ref=m$sample[m$age_norm=="young" & m$condition=="HC"],
       hc=m$sample[m$age_norm==AGES[[ag]] & m$condition=="HC"],
       ne=m$sample[m$age_norm==AGES[[ag]] & m$condition=="EE"])
}

################################################################################
### SHARED: held-out fits (cache first) and CollecTRI
################################################################################

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

uniprot <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
net_signed <- readRDS(COLLECTRI_RDS) %>%
  dplyr::filter(!grepl(uniprot, source), !grepl(uniprot, target), target %in% rownames(pb_all)) %>%
  dplyr::distinct(source, target, .keep_all=TRUE)
REG <- split(net_signed$target, net_signed$source)       # membership, for enrichment
cat(sprintf("CollecTRI: %d TFs, %d interactions\n", length(REG), nrow(net_signed)))

################################################################################
### STEP 2  GENE CALLS (both directions; cached)
################################################################################

calls_rds <- file.path(out_dir, "Step2_GeneCalls.RDS")
if (STEPS[["genes"]] || STEPS[["candidates"]]) {
  banner("STEP 2  Gene calls")
  if (file.exists(calls_rds)) { calls <- readRDS(calls_rds); cat("Loaded:", calls_rds, "\n") } else {
    call_unit <- function(ag, sx, ct) {
      g <- groups_of(ag, sx, ct)
      if (length(g$ref) < 2 || length(g$hc) < 3 || length(g$ne) < 2) return(NULL)
      pf <- dplyr::bind_rows(lapply(seq_along(g$hc), function(k) {
        B <- fit_one(g$hc[-k], g$ref); E <- fit_one(g$ne, g$hc[k]); if (is.null(B) || is.null(E)) return(NULL)
        dplyr::inner_join(B %>% dplyr::select(gene, lfc_B=lfc, p_B=pvalue),
                          E %>% dplyr::select(gene, lfc_E=lfc, p_E=pvalue), by="gene") %>%
          dplyr::filter(!is.na(lfc_B), !is.na(lfc_E), lfc_B != 0)
      }))
      if (nrow(pf) == 0) return(NULL)
      message(sprintf("  calls: %s %s %s", ag, sx, ct))
      pf %>% dplyr::mutate(ag_sig = !is.na(p_B) & p_B < P_AGING, ne_sig = !is.na(p_E) & p_E < P_NE,
                           opp = ag_sig & sign(lfc_E) != sign(lfc_B), same = ag_sig & sign(lfc_E) == sign(lfc_B)) %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(n_fits=dplyr::n(), n_aging=sum(ag_sig), n_opp=sum(opp), n_same=sum(same),
                         n_opp_nesig=sum(opp & ne_sig), n_same_nesig=sum(same & ne_sig),
                         aging_lfc=mean(lfc_B), ne_lfc=mean(lfc_E), .groups="drop") %>%
        dplyr::mutate(T1 = n_aging >= 3 & n_opp >= 3 & n_opp_nesig >= 2, T3 = n_aging >= 2 & n_opp >= 2,
                      A1 = n_aging >= 3 & n_same >= 3 & n_same_nesig >= 2, A3 = n_aging >= 2 & n_same >= 2,
                      age=ag, sex=sx, celltype=ct)
    }
    res <- mclapply(seq_len(nrow(units)), function(i)
      tryCatch(call_unit(units$age[i], units$sex[i], units$celltype[i]), error=function(e) conditionMessage(e)),
      mc.cores=N_CORES)
    errs <- unlist(Filter(is.character, res)); if (length(errs)) { cat("Errors:\n"); print(head(errs)) }
    calls <- dplyr::bind_rows(Filter(is.data.frame, res)); saveRDS(calls, calls_rds)
  }
}

################################################################################
### STEP 3  GENE SETS + CollecTRI (per tier)
################################################################################

build_sets <- function(TIER, BG_MODE) {
  BUF <- TIER; AMP <- sub("^T", "A", TIER); N_AG <- if (TIER == "T1") 3 else 2
  cl  <- calls %>% dplyr::mutate(buf=.data[[BUF]], amp=.data[[AMP]], agechg = n_aging >= N_AG) %>%
    dplyr::select(age, sex, celltype, gene, aging_lfc, ne_lfc, buf, amp, agechg)
  mid <- cl %>% dplyr::filter(age=="middle") %>% dplyr::select(-age)
  old <- cl %>% dplyr::filter(age=="old") %>% dplyr::select(-age)
  both <- dplyr::inner_join(mid, old, by=c("sex","celltype","gene"), suffix=c("_mid","_old"))
  lost <- function(d) d$buf_mid & !d$buf_old
  same_age_dir <- function(d) sign(d$aging_lfc_old) == sign(d$aging_lfc_mid)
  aging_bg <- BG_MODE == "aging"
  dplyr::bind_rows(
    mid %>% dplyr::semi_join(ct_status %>% dplyr::filter(pass_middle), by=c("sex","celltype")) %>%
      dplyr::transmute(sex, celltype, gene, set="Buffered, middle", member=buf, bg=!aging_bg | agechg),
    old %>% dplyr::semi_join(ct_status %>% dplyr::filter(pass_old), by=c("sex","celltype")) %>%
      dplyr::transmute(sex, celltype, gene, set="Buffered, old", member=buf, bg=!aging_bg | agechg),
    both %>% dplyr::transmute(sex, celltype, gene, set="Lost", member=lost(both), bg=!aging_bg | buf_mid),
    both %>% dplyr::transmute(sex, celltype, gene, set="Lost: aging direction changed",
                              member=lost(both) & !same_age_dir(both), bg=!aging_bg | buf_mid),
    both %>% dplyr::transmute(sex, celltype, gene, set="Lost: still opposed, weaker",
                              member=lost(both) & same_age_dir(both) & sign(ne_lfc_old) != sign(aging_lfc_old),
                              bg=!aging_bg | buf_mid),
    both %>% dplyr::transmute(sex, celltype, gene, set="Lost: NE no longer opposes",
                              member=lost(both) & same_age_dir(both) & sign(ne_lfc_old) == sign(aging_lfc_old),
                              bg=!aging_bg | buf_mid),
    both %>% dplyr::transmute(sex, celltype, gene, set="Maintained", member=buf_mid & buf_old, bg=!aging_bg | buf_mid),
    both %>% dplyr::transmute(sex, celltype, gene, set="Age-emergent", member=!buf_mid & buf_old, bg=!aging_bg | agechg_old),
    mid %>% dplyr::transmute(sex, celltype, gene, set="Amplified, middle", member=amp, bg=!aging_bg | agechg),
    old %>% dplyr::transmute(sex, celltype, gene, set="Amplified, old", member=amp, bg=!aging_bg | agechg)) %>%
    dplyr::filter(bg) %>% dplyr::select(-bg) %>%
    dplyr::left_join(ct_status %>% dplyr::select(sex, celltype, ct_status), by=c("sex","celltype")) %>%
    dplyr::left_join(both %>% dplyr::select(sex, celltype, gene, aging_lfc_mid, ne_lfc_mid, aging_lfc_old, ne_lfc_old),
                     by=c("sex","celltype","gene"))
}
enrich_sets <- function(mem) {
  tab <- mem %>% dplyr::group_by(sex, ct_status, celltype, set) %>%
    dplyr::group_modify(function(d, keys) {
      bg <- d$gene; S <- d$gene[d$member]; N <- length(bg); n <- length(S); if (n == 0) return(tibble::tibble())
      dplyr::bind_rows(lapply(names(REG), function(tf) {
        tg <- intersect(REG[[tf]], bg); K <- length(tg); if (K < MIN_REGULON) return(NULL)
        a <- length(intersect(tg, S))
        data.frame(source=tf, a=a, b=n - a, c=K - a, d=N - n - K + a, n_set=n, n_regulon=K,
                   targets=paste(sort(intersect(tg, S)), collapse=","))
      }))
    }) %>% dplyr::ungroup()
  per_ct <- tab %>% dplyr::filter(n_set >= MIN_SET) %>%
    dplyr::mutate(p=stats::phyper(a - 1, n_regulon, b + d, n_set, lower.tail=FALSE),
                  odds_ratio=((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)), pct_targets=round(100 * a / n_regulon, 1)) %>%
    dplyr::group_by(sex, celltype, set) %>% dplyr::mutate(fdr=p.adjust(p, "BH")) %>% dplyr::ungroup() %>%
    dplyr::rename(n_overlap=a) %>% dplyr::select(-b, -c, -d)
  mh <- function(d) {
    d <- d[(d$a + d$c) > 0 & (d$b + d$d) > 0, ]
    if (nrow(d) == 0 || sum(d$a) == 0) return(c(or=NA, p=1))
    arr <- array(as.vector(t(as.matrix(d[, c("a","c","b","d")]))), dim=c(2, 2, nrow(d)))
    if (nrow(d) == 1) { ft <- stats::fisher.test(arr[, , 1], alternative="greater"); return(c(or=unname(ft$estimate), p=ft$p.value)) }
    r <- tryCatch(stats::mantelhaen.test(arr, alternative="greater"), error=function(e) NULL)
    if (is.null(r)) c(or=NA, p=NA) else c(or=unname(r$estimate), p=r$p.value)
  }
  pooled <- tab %>% dplyr::group_by(sex, set, source) %>%
    dplyr::group_modify(function(d, keys) { r <- mh(d)
      tibble::tibble(n_ct=nrow(d), n_overlap=sum(d$a), pct_targets=round(100 * sum(d$a) / sum(d$n_regulon), 1),
                     cts_with_targets=paste(d$celltype[d$a > 0], collapse=", "), or_mh=r[["or"]], p_mh=r[["p"]]) }) %>%
    dplyr::ungroup() %>% dplyr::group_by(sex, set) %>% dplyr::mutate(fdr_mh=p.adjust(p_mh, "BH")) %>% dplyr::ungroup() %>%
    dplyr::arrange(sex, factor(set, levels=SETS), p_mh)
  list(per_ct=per_ct, pooled=pooled)
}

enr <- list()
if (STEPS[["genes"]] || STEPS[["candidates"]]) {
  for (TIER in TIERS_RUN) for (BG_MODE in BG_MODES) {
    banner(sprintf("STEP 3  Gene sets + CollecTRI  (buffered %s, amplified %s; background: %s)",
                   TIER, sub("^T", "A", TIER), BG_MODE))
    tdir <- file.path(out_dir, paste0("Step3_GeneSets_", TIER, "_bg", BG_MODE)); dir.create(tdir, showWarnings=FALSE)
    mem <- build_sets(TIER, BG_MODE)
    wcsv(mem %>% dplyr::filter(member) %>% dplyr::select(-member), "GeneSet_Members.csv", tdir)
    sizes <- mem %>% dplyr::group_by(sex, ct_status, celltype, set) %>%
      dplyr::summarise(n=sum(member), n_background=dplyr::n(), .groups="drop") %>% dplyr::mutate(set=factor(set, levels=SETS))
    wcsv(sizes, "GeneSet_Sizes.csv", tdir)
    cat("\n-- Gene set sizes --\n")
    print(as.data.frame(sizes %>% dplyr::select(-n_background) %>% tidyr::pivot_wider(names_from=set, values_from=n) %>%
      dplyr::arrange(sex, ct_status, factor(celltype, levels=CT_ORDER)), check.names=FALSE), row.names=FALSE)
    e <- enrich_sets(mem); enr[[paste0(TIER, "_bg", BG_MODE)]] <- e
    wcsv(e$per_ct, "CollecTRI_PerCellType.csv", tdir); wcsv(e$pooled, "CollecTRI_Pooled.csv", tdir)
    cat("\n-- Pooled across cell types (Mantel-Haenszel): top 12 TFs per set --\n")
    print(as.data.frame(e$pooled %>% dplyr::filter(n_overlap >= MIN_OVERLAP) %>% dplyr::group_by(sex, set) %>%
      dplyr::slice_min(p_mh, n=12, with_ties=FALSE) %>%
      dplyr::summarise(TFs=paste(sprintf("%s (OR %.1f, p %.2g)", source, or_mh, p_mh), collapse=", "), .groups="drop") %>%
      dplyr::left_join(e$pooled %>% dplyr::group_by(sex, set) %>%
                         dplyr::summarise(n_fdr10=sum(fdr_mh < 0.10, na.rm=TRUE), .groups="drop"), by=c("sex","set")) %>%
      dplyr::relocate(n_fdr10, .after=set) %>%
      dplyr::arrange(sex, factor(set, levels=SETS))), row.names=FALSE)

    if (RUN_GO && !(requireNamespace("clusterProfiler", quietly=TRUE) && requireNamespace("org.Mm.eg.db", quietly=TRUE)))
      cat("\nGO skipped: clusterProfiler / org.Mm.eg.db not installed\n")
    if (RUN_GO && requireNamespace("clusterProfiler", quietly=TRUE) && requireNamespace("org.Mm.eg.db", quietly=TRUE)) {
      suppressPackageStartupMessages({ library(clusterProfiler); library(org.Mm.eg.db) })
      go <- dplyr::bind_rows(lapply(split(mem, list(mem$sex, mem$set), drop=TRUE), function(d) {
        genes <- unique(d$gene[d$member]); if (length(genes) < MIN_SET) return(NULL)
        eg <- tryCatch(clusterProfiler::enrichGO(genes, OrgDb=org.Mm.eg.db, keyType="SYMBOL", ont="BP",
                                                 universe=unique(d$gene), pvalueCutoff=1, qvalueCutoff=1), error=function(e) NULL)
        if (is.null(eg)) return(NULL); r <- as.data.frame(eg); if (!nrow(r)) return(NULL)
        r %>% dplyr::mutate(sex=d$sex[1], set=d$set[1], .before=1) }))
      wcsv(go, "GO_BP_Pooled.csv", tdir)
    }

    d1 <- sizes %>% dplyr::filter(set %in% MAIN_SETS) %>%
      dplyr::mutate(set=factor(set, levels=MAIN_SETS), celltype=factor(celltype, levels=rev(CT_ORDER)))
    sv(ggplot(d1, aes(n, celltype, fill=set)) + geom_col(position=position_dodge(width=0.85), width=0.8) +
         scale_fill_manual(values=SCOL, name=NULL) + facet_grid(ct_status ~ sex, scales="free_y", space="free_y") + th +
         theme(strip.text.y=element_text(angle=0)) +
         labs(title=sprintf("Gene sets per cell type (%s; background %s)", TIER, BG_MODE), x="Number of genes", y=NULL),
       tdir, "F_GeneSet_Sizes", 11, 9)
    d2 <- sizes %>% dplyr::filter(grepl("^Lost:", set)) %>% dplyr::mutate(celltype=factor(celltype, levels=rev(CT_ORDER)))
    sv(ggplot(d2, aes(n, celltype, fill=set)) + geom_col(width=0.75) +
         scale_fill_manual(values=c(`Lost: NE no longer opposes`="#C62828", `Lost: aging direction changed`="#1565C0",
                                    `Lost: still opposed, weaker`="#EF9A9A"), name=NULL) +
         facet_grid(ct_status ~ sex, scales="free_y", space="free_y") + th + theme(strip.text.y=element_text(angle=0)) +
         labs(title="Why buffered genes are lost from middle to old age",
              subtitle="Red = NE response changed; blue = aging program changed", x="Number of genes", y=NULL),
       tdir, "F_LostGenes_ByReason", 10, 7)
    top <- e$pooled %>% dplyr::filter(set %in% MAIN_SETS, n_overlap >= MIN_OVERLAP) %>% dplyr::group_by(sex, set) %>%
      dplyr::slice_min(p_mh, n=8, with_ties=FALSE) %>% dplyr::ungroup()
    tfs <- unique(c(top$source, intersect(WATCH, e$pooled$source)))
    d3 <- e$pooled %>% dplyr::filter(set %in% MAIN_SETS, source %in% tfs) %>%
      dplyr::mutate(set=factor(set, levels=MAIN_SETS), source=factor(source, levels=rev(sort(tfs))),
                    nlp=pmin(-log10(p_mh), 6), sig=fdr_mh < 0.10)
    sv(ggplot(d3, aes(set, source)) + geom_point(aes(size=pmin(or_mh, 6), color=nlp)) +
         geom_point(data=d3 %>% dplyr::filter(sig), shape=1, size=4.5, color="black") +
         scale_color_gradient(low="#FEE8C8", high="#B30000", name="-log10 p") +
         scale_size_continuous(range=c(0.4, 5), name="Odds ratio\n(capped 6)") +
         facet_wrap(~sex) + th + theme(axis.text.x=element_text(angle=40, hjust=1)) +
         labs(title=sprintf("CollecTRI over-representation per gene set (%s; background %s; pooled across cell types)", TIER, BG_MODE),
              subtitle="Ring = FDR < 0.10", x=NULL, y=NULL),
       tdir, "F_CollecTRI_Pooled", 11, max(6, 0.16 * length(tfs) + 2))
  }

  # TFs that hold up under BOTH backgrounds (pooled, per tier)
  for (TIER in TIERS_RUN) {
    ka <- paste0(TIER, "_bgaging"); kt <- paste0(TIER, "_bgtested")
    if (!all(c(ka, kt) %in% names(enr))) next
    rob <- dplyr::inner_join(
      enr[[ka]]$pooled %>% dplyr::select(sex, set, source, n_overlap_aging=n_overlap, or_aging=or_mh, p_aging=p_mh),
      enr[[kt]]$pooled %>% dplyr::select(sex, set, source, n_overlap_tested=n_overlap, or_tested=or_mh, p_tested=p_mh),
      by=c("sex","set","source")) %>%
      dplyr::mutate(robust = p_aging < 0.05 & p_tested < 0.05 & n_overlap_aging >= MIN_OVERLAP,
                    tested_only = p_tested < 0.05 & !(p_aging < 0.05))
    wcsv(rob, paste0("Step3_BackgroundComparison_", TIER, ".csv"))
    cat(sprintf("\n-- %s: TFs enriched under BOTH backgrounds (pooled p < 0.05; top 12 by aging-background p) --\n", TIER))
    print(as.data.frame(rob %>% dplyr::filter(robust, set %in% MAIN_SETS) %>% dplyr::group_by(sex, set) %>%
      dplyr::slice_min(p_aging, n=12, with_ties=FALSE) %>%
      dplyr::summarise(n_robust=dplyr::n(), TFs=paste(sprintf("%s (OR %.1f)", source, or_aging), collapse=", "), .groups="drop") %>%
      dplyr::arrange(sex, factor(set, levels=MAIN_SETS))), row.names=FALSE)
    cat(sprintf("   %s: TFs enriched only with all tested genes as background (aging-responsiveness): %d\n", TIER,
                sum(rob$tested_only & rob$set %in% MAIN_SETS)))
  }
}

################################################################################
### STEP 4  CIRCADIAN TIMING CHECK
################################################################################

if (STEPS[["timing"]]) {
  banner("STEP 4  Circadian timing check")
  cpm <- t(t(as.matrix(pb_all)) / colSums(as.matrix(pb_all))) * 1e6
  cg  <- intersect(CLOCK_GENES, rownames(cpm))
  d <- as.data.frame(t(log2(cpm[cg, , drop=FALSE] + 1))) %>% dplyr::mutate(sample=rownames(.)) %>%
    dplyr::left_join(pb_meta %>% dplyr::select(sample, orig.ident, celltype, sex_norm, age_norm, condition), by="sample") %>%
    tidyr::pivot_longer(dplyr::all_of(cg), names_to="gene", values_to="expr") %>%
    dplyr::group_by(gene, celltype, sex_norm) %>% dplyr::mutate(expr=expr - mean(expr)) %>% dplyr::ungroup() %>%
    dplyr::group_by(gene, mouse=short_id(orig.ident), sex_norm, age_norm, condition) %>%
    dplyr::summarise(z=mean(expr), .groups="drop") %>%
    tidyr::pivot_wider(names_from=gene, values_from=z)
  bm <- intersect(c("Arntl","Bmal1"), names(d))[1]
  d$phase <- rowSums(d[, intersect(c("Dbp","Nr1d1","Per3"), names(d)), drop=FALSE]) - if (!is.na(bm)) d[[bm]] else 0
  wcsv(d, "Step4_Timing_PerMouse.csv")
  cat("\n-- Per mouse (cell-type-centered log2 CPM, averaged over cell types) --\n")
  print(as.data.frame(d %>% dplyr::arrange(sex_norm, age_norm, condition) %>% dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 2)))), row.names=FALSE)
  tt <- d %>% dplyr::group_by(sex_norm, age_norm) %>%
    dplyr::summarise(phase_HC=round(mean(phase[condition=="HC"]), 2), phase_NE=round(mean(phase[condition=="EE"]), 2),
                     p_wilcox=tryCatch(signif(wilcox.test(phase[condition=="HC"], phase[condition=="EE"])$p.value, 2),
                                       error=function(e) NA_real_), .groups="drop")
  wcsv(tt, "Step4_Timing_Summary.csv")
  cat("\n-- Phase index (Dbp + Nr1d1 + Per3 - Bmal1/Arntl): NE vs HC per age --\n"); print(as.data.frame(tt), row.names=FALSE)
  cat("   Reading: NE consistently shifted vs HC at every age (incl. young-matched) -> clock signal likely collection time/activity.\n")
}

################################################################################
### STEP 5  TF ACTIVITY (union cell types, both ages, all labelings; resumable)
################################################################################

tf_dir <- file.path(out_dir, "tf_units")
if (STEPS[["tf_activity"]]) {
  banner("STEP 5  TF activity (union, both ages, all labelings) - slow, resumable")
  tf_activity <- function(B, E) {
    g <- intersect(names(B), names(E)); if (length(g) < 500) return(NULL)
    mat <- cbind(B=B[g], E=E[g]); rownames(mat) <- g
    r <- tryCatch(decoupleR::run_ulm(mat=mat, net=net_signed, .source="source", .target="target", .mor="mor",
                                     minsize=MIN_REGULON), error=function(e) NULL)
    if (is.null(r)) return(NULL)
    r %>% dplyr::filter(statistic=="ulm") %>% dplyr::select(source, condition, score) %>%
      tidyr::pivot_wider(names_from=condition, values_from=score, names_prefix="a_")
  }
  wald <- function(f) if (is.null(f)) NULL else { s <- setNames(f$stat, f$gene); s[!is.na(s)] }
  score_labeling <- function(hc, ne, ref) {
    pf <- dplyr::bind_rows(lapply(seq_along(hc), function(k) {
      a <- tf_activity(wald(fit_one(hc[-k], ref)), wald(fit_one(ne, hc[k]))); if (is.null(a)) NULL else dplyr::mutate(a, k=k) }))
    if (nrow(pf) == 0) return(NULL)
    pf %>% dplyr::group_by(source) %>%
      dplyr::summarise(n_fits=dplyr::n(), S=mean(-a_B * a_E), aging_act=mean(a_B), ne_act=mean(a_E), .groups="drop") %>%
      dplyr::filter(n_fits >= 2)
  }
  for (i in seq_len(nrow(units))) {
    ag <- units$age[i]; sx <- units$sex[i]; ct <- units$celltype[i]
    f_out <- file.path(tf_dir, paste0(ag, "__", sx, "__", gsub("[^A-Za-z0-9]+", "_", ct), ".RDS"))
    if (file.exists(f_out)) { cat(sprintf("  done already: %s %s %s\n", ag, sx, ct)); next }
    g <- groups_of(ag, sx, ct)
    if (length(g$ref) < 2 || length(g$hc) < 3 || length(g$ne) < 3) { cat(sprintf("  skip (n): %s %s %s\n", ag, sx, ct)); next }
    pool <- c(g$hc, g$ne); mouse <- setNames(short_id(g$m$orig.ident[match(pool, g$m$sample)]), pool)
    labs <- combn(pool, length(g$ne), simplify=FALSE)
    cat(sprintf("  scoring %s %s %s (%d labelings) ... ", ag, sx, ct, length(labs)))
    res <- mclapply(labs, function(lab) tryCatch({
      s <- score_labeling(setdiff(pool, lab), lab, g$ref); if (is.null(s)) NULL else
        s %>% dplyr::mutate(age=ag, sex=sx, celltype=ct, real=setequal(lab, g$ne),
                            key=paste(sort(mouse[lab]), collapse="|"), .before=1)
    }, error=function(e) NULL), mc.cores=N_CORES)
    saveRDS(dplyr::bind_rows(Filter(is.data.frame, res)), f_out); cat("saved\n")
  }
}

tf_scores <- NULL
f_units <- list.files(tf_dir, pattern="\\.RDS$", full.names=TRUE)
if (length(f_units)) {
  tf_scores <- dplyr::bind_rows(lapply(f_units, readRDS))
  emp_p <- function(x, null, upper=TRUE) { null <- sort(null); n <- length(null)
    cnt <- if (upper) n - findInterval(x, null, left.open=TRUE) else findInterval(x, null); (1 + cnt) / (1 + n) }
  act_ct <- tf_scores %>% dplyr::group_by(age, sex, celltype) %>%
    dplyr::group_modify(function(d, keys) {
      null <- d$S[!d$real]; r <- d[d$real, ]; if (nrow(r) == 0 || length(null) < 100) return(tibble::tibble())
      own <- d %>% dplyr::group_by(source) %>%
        dplyr::summarise(rank_buf=sum(S >= S[real][1]), rank_amp=sum(S <= S[real][1]), n_labelings=dplyr::n(), .groups="drop")
      r %>% dplyr::select(source, S, aging_act, ne_act) %>%
        dplyr::mutate(p_buf=emp_p(S, null, TRUE), p_amp=emp_p(S, null, FALSE)) %>% dplyr::left_join(own, by="source")
    }) %>% dplyr::ungroup() %>%
    dplyr::mutate(buf = S > 0 & rank_buf <= ACT_RANK, amp = S < 0 & rank_amp <= ACT_RANK) %>%
    dplyr::left_join(ct_status %>% dplyr::select(sex, celltype, ct_status, pass_middle, pass_old), by=c("sex","celltype"))
  act_all <- tf_scores %>% dplyr::group_by(age, sex, key, real, source) %>%
    dplyr::summarise(S_sum=sum(S), n_ct=dplyr::n(), .groups="drop") %>% dplyr::filter(n_ct >= 2) %>%
    dplyr::group_by(age, sex) %>%
    dplyr::group_modify(function(d, keys) {
      null <- d$S_sum[!d$real]; r <- d[d$real, ]; if (nrow(r) == 0) return(tibble::tibble())
      own <- d %>% dplyr::group_by(source) %>%
        dplyr::summarise(rank_buf=sum(S_sum >= S_sum[real][1]), rank_amp=sum(S_sum <= S_sum[real][1]), .groups="drop")
      r %>% dplyr::select(source, S_sum, n_ct) %>%
        dplyr::mutate(p_buf=emp_p(S_sum, null, TRUE), p_amp=emp_p(S_sum, null, FALSE)) %>% dplyr::left_join(own, by="source")
    }) %>% dplyr::ungroup() %>%
    dplyr::mutate(buf = S_sum > 0 & rank_buf <= ACT_RANK, amp = S_sum < 0 & rank_amp <= ACT_RANK)
  wcsv(act_ct, "Step5_TFActivity_PerCellType.csv"); wcsv(act_all, "Step5_TFActivity_AcrossCellTypes.csv")
  cat(sprintf("\nTF activity available for %d units\n", dplyr::n_distinct(paste(act_ct$age, act_ct$sex, act_ct$celltype))))
}

################################################################################
### STEP 6  TRAJECTORIES (per cell type and TF level)
################################################################################

if (STEPS[["trajectory"]] && !is.null(tf_scores)) {
  banner("STEP 6  TF trajectories, middle -> old")
  k <- c("sex","celltype","source")
  mid <- act_ct %>% dplyr::filter(age=="middle"); old <- act_ct %>% dplyr::filter(age=="old")
  traj <- mid %>% dplyr::select(dplyr::all_of(k), ct_status, S_mid=S, buf_mid=buf, amp_mid=amp) %>%
    dplyr::inner_join(old %>% dplyr::select(dplyr::all_of(k), S_old=S, buf_old=buf, amp_old=amp), by=k) %>%
    dplyr::mutate(trajectory=dplyr::case_when(
      buf_mid & buf_old ~ "Maintained", buf_mid & S_old > 0 ~ "Weakened", buf_mid ~ "Lost",
      !buf_mid & buf_old ~ "Age-emergent", amp_mid & amp_old ~ "Amplified, both ages",
      !amp_mid & amp_old ~ "Amplified at old age", amp_mid ~ "Amplified at middle only", TRUE ~ "Neither"),
      flipped_to_amplified = trajectory == "Lost" & amp_old)
  wcsv(traj, "Step6_Trajectory_PerCellType.csv")
  cat("\n-- Per cell type: TF trajectories by cell-type status --\n")
  print(as.data.frame(traj %>% dplyr::filter(trajectory != "Neither") %>% dplyr::count(sex, ct_status, trajectory) %>%
    tidyr::pivot_wider(names_from=trajectory, values_from=n, values_fill=0L)), row.names=FALSE)

  LAB <- c(BB="Maintained", B0="Lost", BA="Lost (flipped to amplified)", `0B`="Age-emergent", AB="Age-emergent",
           AA="Amplified, both ages", `0A`="Amplified at old age", A0="Amplified at middle only", `00`="Not buffered")
  st <- function(d) d %>% dplyr::transmute(sex, source, st=dplyr::case_when(buf ~ "B", amp ~ "A", TRUE ~ "0"), rank_buf, rank_amp)
  tf_level <- st(act_all %>% dplyr::filter(age=="middle")) %>%
    dplyr::inner_join(st(act_all %>% dplyr::filter(age=="old")), by=c("sex","source"), suffix=c("_mid","_old")) %>%
    dplyr::mutate(label=unname(LAB[paste0(st_mid, st_old)]),
                  pos_mid=(rank_amp_mid - rank_buf_mid) / 19, pos_old=(rank_amp_old - rank_buf_old) / 19)
  wcsv(tf_level, "Step6_Trajectory_TFlevel.csv")
  cat("\n-- TF level (across cell types): counts per label --\n")
  print(as.data.frame(tf_level %>% dplyr::count(sex, label) %>% tidyr::pivot_wider(names_from=sex, values_from=n, values_fill=0L)), row.names=FALSE)
  cat("\n-- Watch list --\n")
  print(as.data.frame(tf_level %>% dplyr::filter(source %in% WATCH) %>% dplyr::arrange(sex, label, source) %>%
    dplyr::select(sex, source, label, rank_buf_mid, rank_amp_mid, rank_buf_old, rank_amp_old)), row.names=FALSE)

  LCOL <- c(Maintained="#1B5E20", Lost="#EF9A9A", `Lost (flipped to amplified)`="#C62828", `Age-emergent`="#1565C0",
            `Amplified, both ages`="#6A1B9A", `Amplified at old age`="#BA68C8", `Amplified at middle only`="#D1C4E9",
            `Not buffered`="grey85")
  d6 <- tf_level %>% dplyr::mutate(lab=ifelse(source %in% WATCH, source, NA), label=factor(label, levels=names(LCOL)))
  sv(ggplot(d6, aes(pos_mid, pos_old, color=label)) +
       geom_hline(yintercept=0, color="grey75") + geom_vline(xintercept=0, color="grey75") +
       geom_jitter(width=0.01, height=0.01, size=1, alpha=0.6) +
       geom_text(aes(label=lab), size=2.3, vjust=-0.9, show.legend=FALSE, na.rm=TRUE) +
       scale_color_manual(values=LCOL, name=NULL, drop=FALSE) + facet_wrap(~sex) + coord_equal() + th +
       labs(title="TF-level buffering, middle vs old age (union cell types)",
            subtitle="Position among its own labelings: +1 = most buffered, -1 = most amplified", x="Middle age", y="Old age"),
     out_dir, "Step6_F_TFlevel", 11, 6)
}

################################################################################
### STEP 7  CANDIDATES (enrichment in a set AND matching TF activity)
################################################################################

if (STEPS[["candidates"]] && !is.null(tf_scores) && length(enr)) {
  banner("STEP 7  Candidates")
  a <- act_ct %>% dplyr::select(age, sex, celltype, source, buf, amp) %>%
    tidyr::pivot_wider(names_from=age, values_from=c(buf, amp))
  for (col in c("buf_middle","buf_old","amp_middle","amp_old")) if (!col %in% names(a)) a[[col]] <- NA
  match_set <- function(set, bm, bo, am, ao) dplyr::case_when(
    set == "Buffered, middle" ~ bm %in% TRUE, set == "Buffered, old" ~ bo %in% TRUE,
    grepl("^Lost", set) ~ bm %in% TRUE & !(bo %in% TRUE), set == "Maintained" ~ bm %in% TRUE & bo %in% TRUE,
    set == "Age-emergent" ~ !(bm %in% TRUE) & bo %in% TRUE,
    set == "Amplified, middle" ~ am %in% TRUE, set == "Amplified, old" ~ ao %in% TRUE, TRUE ~ FALSE)
  for (TIER in names(enr)) {
    cand <- enr[[TIER]]$per_ct %>% dplyr::left_join(a, by=c("sex","celltype","source")) %>%
      dplyr::mutate(activity_match=match_set(set, buf_middle, buf_old, amp_middle, amp_old),
                    flipped_to_amplified=grepl("^Lost", set) & amp_old %in% TRUE,
                    enriched=p < 0.05 & n_overlap >= MIN_OVERLAP, candidate=enriched & activity_match)
    wcsv(cand, paste0("Step7_Candidates_", TIER, ".csv"))
    cat(sprintf("\n-- %s: candidates (targets enriched AND activity matches; * = flipped to amplified) --\n", TIER))
    print(as.data.frame(cand %>% dplyr::filter(candidate) %>%
      dplyr::mutate(source=ifelse(flipped_to_amplified, paste0(source, "*"), source)) %>%
      dplyr::group_by(sex, set, celltype) %>% dplyr::arrange(p, .by_group=TRUE) %>%
      dplyr::summarise(TFs=paste(sprintf("%s(%d)", source, n_overlap), collapse=", "), .groups="drop") %>%
      dplyr::arrange(sex, factor(set, levels=SETS), factor(celltype, levels=CT_ORDER))), row.names=FALSE)
    rec <- cand %>% dplyr::filter(candidate, set %in% MAIN_SETS) %>% dplyr::group_by(sex, set, source) %>%
      dplyr::summarise(n_ct=dplyr::n(), cts=paste(celltype, collapse=", "), .groups="drop") %>%
      dplyr::arrange(sex, factor(set, levels=MAIN_SETS), dplyr::desc(n_ct))
    wcsv(rec, paste0("Step7_Candidates_", TIER, "_Recurrence.csv"))
    cat(sprintf("\n-- %s: candidate TFs recurring in >= 2 cell types --\n", TIER))
    print(as.data.frame(rec %>% dplyr::filter(n_ct >= 2)), row.names=FALSE)
    d7 <- cand %>% dplyr::filter(candidate, set %in% MAIN_SETS) %>%
      dplyr::mutate(set=factor(set, levels=MAIN_SETS), celltype=factor(celltype, levels=CT_ORDER))
    if (nrow(d7)) sv(ggplot(d7, aes(celltype, source)) +
         geom_point(aes(color=set, size=n_overlap, shape=flipped_to_amplified)) +
         scale_color_manual(values=SCOL, name=NULL) + scale_shape_manual(values=c(`FALSE`=16, `TRUE`=17), name="Flipped to\namplified") +
         facet_grid(. ~ sex, scales="free_x", space="free_x") + th + theme(axis.text.x=element_text(angle=45, hjust=1)) +
         labs(title=sprintf("Candidate regulators (%s): targets enriched AND TF activity matches", TIER),
              subtitle=paste0("Enrichment p < 0.05 with >= ", MIN_OVERLAP, " targets; activity rank <= ", ACT_RANK), x=NULL, y=NULL),
       out_dir, paste0("Step7_F_Candidates_", TIER), 12, max(5, 0.15 * dplyr::n_distinct(d7$source) + 2))
  }
}

banner(paste("Done. Outputs in:", out_dir))

invisible(NULL)
}

################################################################################
### STEP: methodcheck  (= Buffering_TF_MethodCheck.R)
################################################################################

STEPS$methodcheck <- function() {
################################################################################
### Buffering_TF_MethodCheck.R
###
### Robustness of the TF activity results to the scoring method and the regulon set.
### Main analysis (master pipeline, Step 5) = ULM + CollecTRI. This script reruns the
### identical scoring with:
###   mlm_collectri : multivariate linear model, all TFs fitted jointly, so targets shared
###                   between TFs (NF-kB, AP-1 families) are not credited to every member
###   ulm_dorothea  : same ULM, smaller high-confidence regulons (DoRothEA levels A-C)
### Everything else is unchanged: held-out Wald stats (aging B_k, NE E_k), buffering score
### S = mean_k(-a_B * a_E), all 20 labelings of the aged mice, call = rank <= ACT_RANK of 20,
### same union cell types (master Step1_CellTypes.csv), both ages, each sex.
###
### Outputs (Buffering_Candidates_Master/TF_MethodCheck/):
###   <config>/units/*.RDS               scores per unit (resumable)
###   MethodCheck_Agreement.csv          per age x sex: how many ULM-CollecTRI calls each
###                                      config reproduces (same direction / also called)
###   MethodCheck_Watchlist_PerCellType.csv  headline TFs, call per cell type per config
###   MethodCheck_Watchlist_Pooled.csv   headline TFs, superpooled call + rank per config
################################################################################

Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({
  library(conflicted); library(DESeq2); library(dplyr); library(tidyr); library(parallel); library(decoupleR)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::mutate, .quiet=TRUE); conflicts_prefer(base::intersect, .quiet=TRUE)
  conflicts_prefer(base::setdiff, .quiet=TRUE); conflicts_prefer(base::union, .quiet=TRUE)
  conflicts_prefer(base::setequal, .quiet=TRUE); conflicts_prefer(base::Filter, .quiet=TRUE)
})
if (requireNamespace("RhpcBLASctl", quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(1)

################################################################################
### PARAMETERS
################################################################################

base_path <- NE_BASE
# DOROTHEA_RDS (set at the top) is downloaded once and saved; without internet, download it elsewhere:
#   saveRDS(decoupleR::get_dorothea(organism="mouse", levels=c("A","B","C")), DOROTHEA_RDS)
CACHE_DIR <- file.path(base_path, "HeldOut_FitCache")
master    <- MASTER_DIR
out_dir   <- file.path(master, "TF_MethodCheck")
CONFIGS   <- list(mlm_collectri=c(method="mlm", net="collectri"),
                  ulm_dorothea =c(method="ulm", net="dorothea"))
N_CORES   <- 16
MIN_REGULON <- 5
ACT_RANK  <- 5
WATCH <- c("Nr3c1","Creb1","Sirt1","Atf6","Foxo1","Foxo3","Egr3","Rela","Relb","Nfkb1","Nfkb2",
           "Atf4","Ddit3","Hsf1","Jun","Fos","Egr1")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

################################################################################
### INPUTS
################################################################################

pb <- readRDS(file.path(base_path, "Buffering_LOO_Split", "pseudobulk_counts_and_meta.RDS"))
pb_all <- pb$pb_all; pb_meta <- pb$pb_meta
age_vals <- unique(pb_meta$age_norm)
AGES <- c(middle = age_vals[!age_vals %in% c("young","old")][1], old = "old")
short_id <- function(x) sub(".*--", "", x)

ct_status <- read.csv(file.path(master, "Step1_CellTypes.csv"), stringsAsFactors=FALSE)
units <- tidyr::crossing(ct_status %>% dplyr::select(sex, celltype), age=c("middle","old"))
cat(sprintf("Units (union cell types x age): %d\n", nrow(units)))

fit_key <- function(num, den) paste(paste(sort(den), collapse="+"), paste(sort(num), collapse="+"), sep=" -> ")
FIT_CACHE <- list()
if (dir.exists(CACHE_DIR)) for (f in list.files(CACHE_DIR, pattern="\\.RDS$", full.names=TRUE)) FIT_CACHE <- c(FIT_CACHE, readRDS(f)$fits)
cat(sprintf("Fit cache: %d fits\n", length(FIT_CACHE)))
fit_one <- function(num, den) {
  hit <- FIT_CACHE[[fit_key(num, den)]]; if (!is.null(hit)) return(hit[, c("gene","stat")])
  cols <- c(den, num); cnt <- round(as.matrix(pb_all[, cols, drop=FALSE])); cnt <- cnt[rowSums(cnt >= 10) >= 2, , drop=FALSE]
  if (nrow(cnt) < 10) return(NULL)
  cd <- data.frame(grp=factor(rep(c("den","num"), c(length(den), length(num))), levels=c("den","num")), row.names=cols)
  tryCatch({
    dds <- suppressMessages(DESeqDataSetFromMatrix(cnt, cd, design=~grp))
    dds <- tryCatch(DESeq(dds, quiet=TRUE), error=function(e) {
      d <- estimateSizeFactors(dds); d <- estimateDispersionsGeneEst(d); dispersions(d) <- mcols(d)$dispGeneEst
      nbinomWaldTest(d, quiet=TRUE) })
    r <- results(dds, name="grp_num_vs_den"); data.frame(gene=rownames(r), stat=r$stat, stringsAsFactors=FALSE)
  }, error=function(e) NULL)
}
wald <- function(f) if (is.null(f)) NULL else { s <- setNames(f$stat, f$gene); s[!is.na(s)] }

uniprot <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
clean_net <- function(n) {
  n <- n %>% dplyr::filter(!grepl(uniprot, source), !grepl(uniprot, target), target %in% rownames(pb_all)) %>%
    dplyr::distinct(source, target, .keep_all=TRUE)
  # drop regulons identical to another (same targets and signs): they are not separable in a joint model
  sig <- n %>% dplyr::arrange(target) %>% dplyr::group_by(source) %>%
    dplyr::summarise(key=paste(target, mor, collapse=";"), .groups="drop")
  dup <- sig$source[duplicated(sig$key)]
  if (length(dup)) cat(sprintf("  dropped %d regulons identical to another: %s\n", length(dup), paste(head(dup, 10), collapse=", ")))
  n %>% dplyr::filter(!source %in% dup)
}
NETS <- list()
NETS$collectri <- clean_net(readRDS(COLLECTRI_RDS))
get_dorothea_mouse <- function() {
  ok <- function(d) !is.null(d) && NROW(d) > 0
  # 1) the 'dorothea' Bioconductor package ships the mouse table (no internet needed)
  if (requireNamespace("dorothea", quietly=TRUE)) {
    e <- new.env(); utils::data("dorothea_mm", package="dorothea", envir=e)
    d <- e$dorothea_mm
    if (ok(d)) { cat("DoRothEA: mouse table from the 'dorothea' package\n")
      return(d %>% dplyr::filter(confidence %in% c("A","B","C")) %>% dplyr::transmute(source=tf, target, mor, confidence)) }
  }
  # 2) OmniPath, mouse
  d <- tryCatch(decoupleR::get_dorothea(organism="mouse", levels=c("A","B","C")), error=function(e) NULL)
  if (ok(d)) { cat("DoRothEA: mouse table from OmniPath\n"); return(d) }
  # 3) OmniPath, human, symbols converted to mouse capitalisation (1:1 orthologs only; approximate)
  d <- tryCatch(decoupleR::get_dorothea(organism="human", levels=c("A","B","C")), error=function(e) NULL)
  if (!ok(d)) stop("Could not obtain DoRothEA. Install the 'dorothea' package (BiocManager::install('dorothea')) and rerun.")
  tc <- function(x) paste0(substr(x, 1, 1), tolower(substring(x, 2)))
  cat("DoRothEA: human table from OmniPath, symbols converted to mouse capitalisation (approximate)\n")
  d %>% dplyr::mutate(source=tc(source), target=tc(target))
}
if (any(sapply(CONFIGS, `[[`, "net") == "dorothea")) {
  if (!file.exists(DOROTHEA_RDS)) saveRDS(get_dorothea_mouse(), DOROTHEA_RDS)
  NETS$dorothea <- clean_net(readRDS(DOROTHEA_RDS))
}
for (nm in names(NETS)) cat(sprintf("%s: %d TFs, %d interactions\n", nm, dplyr::n_distinct(NETS[[nm]]$source), nrow(NETS[[nm]])))

################################################################################
### SCORING (identical to master Step 5 except method / network)
################################################################################

tf_activity <- function(B, E, method, net) {
  g <- intersect(names(B), names(E)); if (length(g) < 500) return(NULL)
  mat <- cbind(B=B[g], E=E[g]); rownames(mat) <- g
  f <- if (method == "mlm") decoupleR::run_mlm else decoupleR::run_ulm
  r <- tryCatch(f(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=MIN_REGULON),
                error=function(e) { message("  decoupleR error: ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  r %>% dplyr::filter(statistic == method) %>% dplyr::select(source, condition, score) %>%
    tidyr::pivot_wider(names_from=condition, values_from=score, names_prefix="a_")
}
groups_of <- function(ag, sx, ct) {
  m <- pb_meta[pb_meta$celltype == ct & pb_meta$sex_norm == sx, ]
  list(m=m, ref=m$sample[m$age_norm == "young" & m$condition == "HC"],
       hc=m$sample[m$age_norm == AGES[[ag]] & m$condition == "HC"], ne=m$sample[m$age_norm == AGES[[ag]] & m$condition == "EE"])
}
for (cfg in names(CONFIGS)) {
  method <- CONFIGS[[cfg]][["method"]]; net <- NETS[[CONFIGS[[cfg]][["net"]]]]
  udir <- file.path(out_dir, cfg, "units"); dir.create(udir, recursive=TRUE, showWarnings=FALSE)
  cat(sprintf("\n== %s ==\n", cfg))
  for (i in seq_len(nrow(units))) {
    ag <- units$age[i]; sx <- units$sex[i]; ct <- units$celltype[i]
    f_out <- file.path(udir, paste0(ag, "__", sx, "__", gsub("[^A-Za-z0-9]+", "_", ct), ".RDS"))
    if (file.exists(f_out)) next
    g <- groups_of(ag, sx, ct)
    if (length(g$ref) < 2 || length(g$hc) < 3 || length(g$ne) < 3) next
    pool <- c(g$hc, g$ne); mouse <- setNames(short_id(g$m$orig.ident[match(pool, g$m$sample)]), pool)
    labs <- combn(pool, length(g$ne), simplify=FALSE)
    t0 <- Sys.time()
    res <- mclapply(labs, function(lab) tryCatch({
      hc <- setdiff(pool, lab)
      pf <- dplyr::bind_rows(lapply(seq_along(hc), function(k) {
        a <- tf_activity(wald(fit_one(hc[-k], g$ref)), wald(fit_one(lab, hc[k])), method, net)
        if (is.null(a)) NULL else dplyr::mutate(a, k=k) }))
      if (!nrow(pf)) return(NULL)
      pf %>% dplyr::group_by(source) %>%
        dplyr::summarise(n_fits=dplyr::n(), S=mean(-a_B * a_E), .groups="drop") %>% dplyr::filter(n_fits >= 2) %>%
        dplyr::mutate(age=ag, sex=sx, celltype=ct, real=setequal(lab, g$ne), key=paste(sort(mouse[lab]), collapse="|"), .before=1)
    }, error=function(e) NULL), mc.cores=N_CORES)
    saveRDS(dplyr::bind_rows(Filter(is.data.frame, res)), f_out)
    cat(sprintf("  %s %s %-20s %.1f min\n", ag, sx, ct, as.numeric(difftime(Sys.time(), t0, units="mins"))))
  }
}

################################################################################
### SUMMARIES vs the main analysis (ULM + CollecTRI)
################################################################################

calls_ct <- function(sc) sc %>% dplyr::group_by(age, sex, celltype, source) %>%
  dplyr::summarise(rank_buf=sum(S >= S[real][1]), rank_amp=sum(S <= S[real][1]), n_lab=dplyr::n(), S=S[real][1], .groups="drop") %>%   # ranks BEFORE S is overwritten
  dplyr::filter(!is.na(S)) %>%
  dplyr::mutate(call=dplyr::case_when(S > 0 & rank_buf <= ACT_RANK ~ "buffered", S < 0 & rank_amp <= ACT_RANK ~ "amplified", TRUE ~ "-"))
calls_pool <- function(sc) sc %>% dplyr::group_by(age, sex, key, real, source) %>%
  dplyr::summarise(S_sum=sum(S), n_ct=dplyr::n(), .groups="drop") %>% dplyr::filter(n_ct >= 2) %>%
  dplyr::group_by(age, sex, source) %>%
  dplyr::summarise(rank_buf=sum(S_sum >= S_sum[real][1]), rank_amp=sum(S_sum <= S_sum[real][1]), S_sum=S_sum[real][1], .groups="drop") %>%
  dplyr::filter(!is.na(S_sum)) %>%
  dplyr::mutate(call=dplyr::case_when(S_sum > 0 & rank_buf <= ACT_RANK ~ "buffered", S_sum < 0 & rank_amp <= ACT_RANK ~ "amplified", TRUE ~ "-"))

main_ct <- read.csv(file.path(master, "Step5_TFActivity_PerCellType.csv"), stringsAsFactors=FALSE) %>%
  dplyr::transmute(age, sex, celltype, source, S_main=S, call_main=dplyr::case_when(buf ~ "buffered", amp ~ "amplified", TRUE ~ "-"))
main_pool <- read.csv(file.path(master, "Step5_TFActivity_AcrossCellTypes.csv"), stringsAsFactors=FALSE) %>%
  dplyr::transmute(age, sex, source, call_main=dplyr::case_when(buf ~ "buffered", amp ~ "amplified", TRUE ~ "-"),
                   rank_main=ifelse(S_sum > 0, rank_buf, rank_amp))

agree <- list(); agree_pool <- list(); wl_ct <- list(); wl_pool <- list()
for (cfg in names(CONFIGS)) {
  sc <- dplyr::bind_rows(lapply(list.files(file.path(out_dir, cfg, "units"), pattern="\\.RDS$", full.names=TRUE), readRDS))
  if (!nrow(sc)) next
  cc <- calls_ct(sc); cp <- calls_pool(sc)
  j <- main_ct %>% dplyr::inner_join(cc %>% dplyr::select(age, sex, celltype, source, S, call), by=c("age","sex","celltype","source"))
  agree[[cfg]] <- j %>% dplyr::filter(call_main != "-") %>% dplyr::group_by(age, sex) %>%
    dplyr::summarise(n_main_calls=dplyr::n(), pct_same_sign=round(100 * mean(sign(S) == sign(S_main)), 1),
                     pct_also_called=round(100 * mean(call == call_main), 1),
                     rho_all_TFs=round(cor(j$S[j$age == age[1] & j$sex == sex[1]], j$S_main[j$age == age[1] & j$sex == sex[1]],
                                           method="spearman", use="complete.obs"), 2), .groups="drop") %>%
    dplyr::mutate(config=cfg, .before=1)
  write.csv(j %>% dplyr::rename(S_cfg=S, call_cfg=call), file.path(out_dir, paste0("MethodCheck_AllTFs_PerCellType_", cfg, ".csv")), row.names=FALSE)
  jp <- main_pool %>% dplyr::full_join(cp %>% dplyr::transmute(age, sex, source, call_cfg=call, rank_cfg=ifelse(S_sum > 0, rank_buf, rank_amp)),
                                       by=c("age","sex","source"))
  write.csv(jp, file.path(out_dir, paste0("MethodCheck_AllTFs_Pooled_", cfg, ".csv")), row.names=FALSE)
  agree_pool[[cfg]] <- jp %>% dplyr::filter(!is.na(call_main), call_main != "-") %>% dplyr::group_by(age, sex) %>%
    dplyr::summarise(main_pooled_calls=dplyr::n(), scored=sum(!is.na(call_cfg)),
                     reproduced=sum(call_cfg == call_main, na.rm=TRUE),
                     opposite=sum(call_cfg != "-" & call_cfg != call_main, na.rm=TRUE), .groups="drop") %>%
    dplyr::mutate(pct_reproduced_of_scored=round(100 * reproduced / pmax(scored, 1), 1), config=cfg, .before=1)
  wl_ct[[cfg]] <- j %>% dplyr::filter(source %in% WATCH, call_main != "-" | call != "-") %>% dplyr::mutate(config=cfg)
  wl_pool[[cfg]] <- main_pool %>% dplyr::filter(source %in% WATCH) %>%
    dplyr::left_join(cp %>% dplyr::transmute(age, sex, source, call, rank=ifelse(S_sum > 0, rank_buf, rank_amp)), by=c("age","sex","source")) %>%
    dplyr::mutate(config=cfg)
}
agree <- dplyr::bind_rows(agree); write.csv(agree, file.path(out_dir, "MethodCheck_Agreement.csv"), row.names=FALSE)
wl_ct <- dplyr::bind_rows(wl_ct); write.csv(wl_ct, file.path(out_dir, "MethodCheck_Watchlist_PerCellType.csv"), row.names=FALSE)
wl_pool <- dplyr::bind_rows(wl_pool); write.csv(wl_pool, file.path(out_dir, "MethodCheck_Watchlist_Pooled.csv"), row.names=FALSE)

cat("\n=== Agreement with the main analysis (ULM + CollecTRI) ===\n")
cat("pct_same_sign = % of main TF calls with the same direction; pct_also_called = % also called (rank <= 5);",
    "rho_all_TFs = Spearman of S over all TFs\n")
print(as.data.frame(agree), row.names=FALSE)

agree_pool <- dplyr::bind_rows(agree_pool); write.csv(agree_pool, file.path(out_dir, "MethodCheck_Agreement_Pooled.csv"), row.names=FALSE)
cat("\n=== ALL TFs, superpooled: main calls reproduced by each config ===\n")
cat("scored = main calls whose TF exists in that config's network; opposite = called in the other direction\n")
print(as.data.frame(agree_pool), row.names=FALSE)

cat("\n=== Headline TFs, superpooled (call [rank of 20]) ===\n")
fmt <- function(call, rank) ifelse(is.na(call), "not scored", ifelse(call == "-", sprintf("- [%s]", rank), sprintf("%s [%s]", call, rank)))
print(as.data.frame(wl_pool %>% dplyr::mutate(v=fmt(call, rank), main=fmt(call_main, rank_main)) %>%
  dplyr::select(sex, age, source, main, config, v) %>% tidyr::pivot_wider(names_from=config, values_from=v) %>%
  dplyr::arrange(sex, source, age)), row.names=FALSE)

cat("\n=== Headline TFs per cell type: called by main vs each config (counts) ===\n")
print(as.data.frame(wl_ct %>% dplyr::group_by(config, sex, age, source) %>%
  dplyr::summarise(main_calls=sum(call_main != "-"), reproduced=sum(call_main != "-" & call == call_main),
                   same_sign=sum(call_main != "-" & sign(S) == sign(S_main)), new_calls=sum(call_main == "-" & call != "-"), .groups="drop") %>%
  dplyr::arrange(config, sex, source, age)), row.names=FALSE)
cat("\nOutputs in:", out_dir, "\n")

invisible(NULL)
}

################################################################################
### STEP: atlas  (= Buffering_TF_Atlas.R)
################################################################################

STEPS$atlas <- function() {
################################################################################
### Buffering_TF_Atlas.R
###
### Unbiased TF atlas (NE, each sex): EVERY TF x cell type, in each of seven categories,
### ranked by a fixed evidence score. Nothing is removed; robustness and expression are
### shown as annotations only.
###
### Categories (from the main analysis: ULM + CollecTRI, rank <= 5 of 20 per age):
###   Buffered, middle    Buffered, old    Amplified, middle    Amplified, old
###   Maintained (buffered at both ages)   Lost (buffered at middle, not at old; flag if flipped
###   to amplified)   Age-emergent (buffered at old, not at middle)
###   A TF x cell type can belong to several (e.g. "Buffered, middle" and "Lost").
### Evidence per TF x cell type (e = -log10 empirical p vs relabelings, from Step 5):
###   Buffered / Amplified at one age : e at that age
###   Maintained   : min(e_buf middle, e_buf old)
###   Lost         : e_buf middle - e_buf old (>= 0), + e_amp old if it flipped
###   Age-emergent : e_buf old - e_buf middle (>= 0), + e_amp middle if it flipped
###   TF score per category = sum over cell types (breadth x strength).
### Annotations (NOT filters):
###   tier      : does the same category hold with MLM + CollecTRI and ULM + DoRothEA?
###               all three / current + MLM (not in DoRothEA) / current + MLM / current + DoRothEA / current only
###   expressed : TF mRNA detected (>= MIN_COUNT counts) in >= MIN_MICE mice of that cell type at that
###               age (HC or NE); trajectories: at either age. NA = TF not a gene in the data.
###
### Inputs : Buffering_Candidates_Master/ (Step5, Step6 CSVs; TF_MethodCheck/*AllTFs_PerCellType*.csv)
###          Buffering_LOO_Split/pseudobulk_counts_and_meta.RDS
### Outputs: Buffering_Candidates_Master/TF_Atlas/
###          TF_Atlas_Long.csv        one row per sex x category x cell type x TF
###          TF_Atlas_Ranking.csv     one row per sex x category x TF (score, cell types, robust and
###                                   expressed counts) - all TFs, ranked
###          TF_Atlas_<sex>_<category>.pdf/png   top N_TOP TFs per panel (full list in the CSV)
################################################################################

suppressPackageStartupMessages({
  library(conflicted); library(dplyr); library(tidyr); library(ggplot2)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::mutate, .quiet=TRUE)
})

base_path <- NE_BASE
master  <- MASTER_DIR
mc_dir  <- file.path(master, "TF_MethodCheck")
out_dir <- file.path(master, "TF_Atlas"); dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
MIN_COUNT <- 3; MIN_MICE <- 2     # expression annotation (lenient: TFs are lowly captured in snRNA-seq)
N_TOP <- 30                       # rows per figure panel (all TFs are in the CSV)

CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN","Astrocytes","Oligodendrocytes","OPC",
              "Microglia","Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell")
CT_SHORT <- c("Astrocytes"="Astro","Oligodendrocytes"="Oligo","Cajal-Retzius Cell"="CR cell","Endothelial"="Endo",
              "Microglia"="Micro","Pericyte"="Peri")
short_ct <- function(x) ifelse(x %in% names(CT_SHORT), CT_SHORT[x], x)
CATS <- c("Buffered, middle","Buffered, old","Amplified, middle","Amplified, old","Maintained","Lost","Age-emergent")
TIERS <- c("all three","current + MLM (not in DoRothEA)","current + MLM","current + DoRothEA","current only")
TIER_COL <- c(`all three`="#08306B", `current + MLM (not in DoRothEA)`="#2171B5", `current + MLM`="#6BAED6",
              `current + DoRothEA`="#9ECAE1", `current only`="#D9D9D9")
ev <- function(p) -log10(pmax(p, 1e-6))

################################################################################
### Calls per method, wide by age
################################################################################

act <- read.csv(file.path(master, "Step5_TFActivity_PerCellType.csv"), stringsAsFactors=FALSE)
main <- act %>% dplyr::transmute(age, sex, celltype, source, call=dplyr::case_when(buf ~ "B", amp ~ "A", TRUE ~ "0"),
                                 eB=ev(p_buf), eA=ev(p_amp))
cfg <- function(f) {
  p <- file.path(mc_dir, f); if (!file.exists(p)) { message("missing: ", p); return(NULL) }
  read.csv(p, stringsAsFactors=FALSE) %>%
    dplyr::transmute(age, sex, celltype, source, call=dplyr::recode(call_cfg, buffered="B", amplified="A", `-`="0"))
}
mlm  <- cfg("MethodCheck_AllTFs_PerCellType_mlm_collectri.csv")
doro <- cfg("MethodCheck_AllTFs_PerCellType_ulm_dorothea.csv")
k <- c("sex","celltype","source")
wide <- function(d, keep_e=FALSE) {
  cols <- if (keep_e) c("call","eB","eA") else "call"
  d %>% dplyr::select(dplyr::all_of(c(k, "age", cols))) %>%
    tidyr::pivot_wider(names_from=age, values_from=dplyr::all_of(cols), names_sep="_")
}
W  <- wide(main, TRUE)
# one value column -> pivot_wider names the columns just "middle" / "old"
Wm <- wide(mlm)  %>% dplyr::rename(m_mlm=middle,  o_mlm=old)
Wd <- wide(doro) %>% dplyr::rename(m_doro=middle, o_doro=old)
X <- W %>% dplyr::left_join(Wm, by=k) %>% dplyr::left_join(Wd, by=k)

member <- function(cat, m, o) {
  m <- ifelse(is.na(m), "0", m); o <- ifelse(is.na(o), "0", o)
  switch(cat,
    "Buffered, middle"  = m == "B", "Buffered, old" = o == "B",
    "Amplified, middle" = m == "A", "Amplified, old" = o == "A",
    "Maintained" = m == "B" & o == "B", "Lost" = m == "B" & o != "B", "Age-emergent" = m != "B" & o == "B")
}
evidence <- function(cat, d) with(d, switch(cat,
  "Buffered, middle" = eB_middle, "Buffered, old" = eB_old, "Amplified, middle" = eA_middle, "Amplified, old" = eA_old,
  "Maintained"   = pmin(eB_middle, eB_old),
  "Lost"         = pmax(eB_middle - eB_old, 0) + ifelse(call_old == "A", eA_old, 0),
  "Age-emergent" = pmax(eB_old - eB_middle, 0) + ifelse(call_middle == "A", eA_middle, 0)))
age_of <- c("Buffered, middle"="middle","Buffered, old"="old","Amplified, middle"="middle","Amplified, old"="old",
            "Maintained"="either","Lost"="either","Age-emergent"="either")

long <- dplyr::bind_rows(lapply(CATS, function(cat) {
  d <- X[member(cat, X$call_middle, X$call_old) %in% TRUE, ]
  if (!nrow(d)) return(NULL)
  in_mlm  <- member(cat, d$m_mlm, d$o_mlm); doro_scored <- !(is.na(d$m_doro) & is.na(d$o_doro))
  in_doro <- member(cat, d$m_doro, d$o_doro)
  d %>% dplyr::mutate(category=cat, evidence=evidence(cat, d),
    flipped = cat == "Lost" & call_old == "A",
    tier = dplyr::case_when(in_mlm & doro_scored & in_doro ~ TIERS[1], in_mlm & !doro_scored ~ TIERS[2],
                            in_mlm ~ TIERS[3], doro_scored & in_doro ~ TIERS[4], TRUE ~ TIERS[5]),
    age_for_expr = age_of[[cat]])
}))

################################################################################
### Expression annotation (pseudobulk; HC and NE mice of that age)
################################################################################

pb <- readRDS(file.path(base_path, "Buffering_LOO_Split", "pseudobulk_counts_and_meta.RDS"))
cnt <- as.matrix(pb$pb_all); meta <- pb$pb_meta
age_vals <- unique(meta$age_norm); AGES <- c(middle = age_vals[!age_vals %in% c("young","old")][1], old = "old")
tfs <- intersect(unique(long$source), rownames(cnt))
cpm <- t(t(cnt[tfs, , drop=FALSE]) / colSums(cnt)) * 1e6
expr <- dplyr::bind_rows(lapply(names(AGES), function(ag) {
  mm <- meta %>% dplyr::filter(age_norm == AGES[[ag]], condition %in% c("HC","EE"))
  dplyr::bind_rows(lapply(split(mm$sample, paste(mm$sex_norm, mm$celltype, sep="|")), function(s) {
    data.frame(source=tfs, n_det=rowSums(cnt[tfs, s, drop=FALSE] >= MIN_COUNT), n_mice=length(s),
               mean_cpm=rowMeans(cpm[, s, drop=FALSE]))
  }), .id="grp") %>% tidyr::separate(grp, c("sex","celltype"), sep="\\|") %>% dplyr::mutate(age=ag)
})) %>% dplyr::mutate(expressed = n_det >= MIN_MICE)
ex_w <- expr %>% dplyr::select(sex, celltype, source, age, expressed, mean_cpm) %>%
  tidyr::pivot_wider(names_from=age, values_from=c(expressed, mean_cpm))
long <- long %>% dplyr::left_join(ex_w, by=c("sex","celltype","source")) %>%
  dplyr::mutate(expressed = dplyr::case_when(!source %in% tfs ~ NA,
                                             age_for_expr == "middle" ~ expressed_middle, age_for_expr == "old" ~ expressed_old,
                                             TRUE ~ expressed_middle | expressed_old))
write.csv(long %>% dplyr::select(sex, category, celltype, source, evidence, tier, flipped, expressed,
                                 call_middle, call_old, mean_cpm_middle, mean_cpm_old) %>%
            dplyr::arrange(sex, factor(category, levels=CATS), source, factor(celltype, levels=CT_ORDER)),
          file.path(out_dir, "TF_Atlas_Long.csv"), row.names=FALSE)

################################################################################
### Ranking (all TFs) and figures (top N_TOP per panel)
################################################################################

rank_tab <- long %>% dplyr::group_by(sex, category, source) %>%
  dplyr::summarise(score=sum(evidence, na.rm=TRUE), n_ct=dplyr::n(),
                   n_ct_all3=sum(tier == TIERS[1]), n_ct_robust=sum(tier %in% TIERS[1:2]),
                   n_ct_expressed=sum(expressed %in% TRUE), n_flipped=sum(flipped),
                   celltypes=paste(short_ct(celltype[order(-evidence)]), collapse=", "),
                   celltypes_robust=paste(short_ct(celltype[tier %in% TIERS[1:2]]), collapse=", "), .groups="drop") %>%
  dplyr::group_by(sex, category) %>% dplyr::arrange(dplyr::desc(score), .by_group=TRUE) %>%
  dplyr::mutate(rank=dplyr::row_number()) %>% dplyr::ungroup() %>%
  dplyr::arrange(sex, factor(category, levels=CATS), rank)
write.csv(rank_tab, file.path(out_dir, "TF_Atlas_Ranking.csv"), row.names=FALSE)
cat("TFs per sex x category (all kept):\n")
print(as.data.frame(rank_tab %>% dplyr::count(sex, category) %>% tidyr::pivot_wider(names_from=sex, values_from=n)), row.names=FALSE)

th <- theme_classic(base_size=7) +
  theme(plot.title=element_text(face="bold", size=8), plot.subtitle=element_text(color="grey35", size=5.8),
        axis.line=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=5.5), legend.position="bottom",
        legend.box="vertical", legend.key.size=unit(3, "mm"), legend.text=element_text(size=5.5))
for (sx in c("female","male")) for (cat in CATS) {
  top <- rank_tab %>% dplyr::filter(sex == sx, category == cat, rank <= N_TOP)
  if (!nrow(top)) next
  cts <- intersect(CT_ORDER, unique(long$celltype[long$sex == sx])); K <- length(cts)
  d <- long %>% dplyr::filter(sex == sx, category == cat) %>% dplyr::inner_join(top %>% dplyr::select(source, rank), by="source") %>%
    dplyr::mutate(x=match(celltype, cts), y=-rank, tier=factor(tier, levels=TIERS),
                  expr_lab=factor(dplyr::case_when(expressed %in% TRUE ~ "detected", expressed %in% FALSE ~ "not detected",
                                                   TRUE ~ "not a gene in data"), levels=c("detected","not detected","not a gene in data")))
  smax <- max(top$score)
  bars <- top %>% dplyr::mutate(y=-rank, x0=K + 0.9, x1=K + 0.9 + 4.5 * score / smax)
  p <- ggplot() +
    geom_hline(data=top, aes(yintercept=-rank), color="grey94", linewidth=0.25) +
    geom_point(data=d, aes(x, y, size=evidence, fill=tier, shape=expr_lab, color=flipped), stroke=0.45) +
    geom_rect(data=bars, aes(xmin=x0, xmax=x1, ymin=y - 0.36, ymax=y + 0.36), fill="grey75") +
    geom_text(data=bars, aes(x=x1 + 0.15, y=y, label=sprintf("%.1f | %d CT, %d robust", score, n_ct, n_ct_robust)),
              hjust=0, size=1.6, color="grey30") +
    geom_text(data=bars, aes(x=0.35, y=y, label=paste0(rank, ". ", source)), hjust=1, size=1.9, fontface="italic") +
    scale_fill_manual(values=TIER_COL, name="Holds with other methods", drop=FALSE) +
    scale_shape_manual(values=c(detected=21, `not detected`=24, `not a gene in data`=22), name="TF mRNA", drop=FALSE) +
    scale_color_manual(values=c(`FALSE`="grey25", `TRUE`="#B2182B"), labels=c(`FALSE`="", `TRUE`="flipped to amplified"),
                       name=NULL, guide=if (cat == "Lost") "legend" else "none") +
    scale_size_continuous(range=c(0.6, 3), name="Evidence\n-log10 p") +
    scale_x_continuous(breaks=seq_len(K), labels=short_ct(cts), expand=expansion(add=c(0.4, 0.4))) +
    scale_y_continuous(expand=expansion(add=0.6)) + coord_cartesian(xlim=c(-2.8, K + 8), clip="off") + th +
    guides(fill=guide_legend(nrow=2, override.aes=list(size=2.5, shape=21)), shape=guide_legend(override.aes=list(size=2.2))) +
    labs(title=sprintf("%s TFs, NE %s (all TFs ranked; top %d shown)", cat, ifelse(sx == "female", "females", "males"), N_TOP),
         subtitle="Rank = summed evidence across cell types; dot = cell type in this category (size = evidence; colour = method robustness; shape = TF mRNA)",
         x=NULL, y=NULL)
  nm <- paste0("TF_Atlas_", sx, "_", gsub("[^A-Za-z]+", "_", cat))
  ggsave(file.path(out_dir, paste0(nm, ".pdf")), p, width=150 / 25.4, height=(40 + 3.3 * nrow(top)) / 25.4, limitsize=FALSE)
  ggsave(file.path(out_dir, paste0(nm, ".png")), p, width=150 / 25.4, height=(40 + 3.3 * nrow(top)) / 25.4, dpi=400, bg="white", limitsize=FALSE)
}
################################################################################
### SUPERPOOLED ATLAS (score summed across cell types, ranked vs relabelings; Step 5)
###   Same seven categories, same evidence rules, applied to the superpooled calls.
###   Robustness from the superpooled calls of MLM + CollecTRI and ULM + DoRothEA.
###   Expression annotation: number of (union) cell types where the TF mRNA is detected.
###   Outputs: TF_Atlas_Superpooled.csv (all TFs) and TF_Atlas_Superpooled_<sex>.pdf/png
################################################################################

pool <- read.csv(file.path(master, "Step5_TFActivity_AcrossCellTypes.csv"), stringsAsFactors=FALSE) %>%
  dplyr::transmute(age, sex, source, call=dplyr::case_when(buf ~ "B", amp ~ "A", TRUE ~ "0"),
                   eB=ev(p_buf), eA=ev(p_amp), rank=ifelse(S_sum > 0, rank_buf, rank_amp))
pcfg <- function(f) {
  p <- file.path(mc_dir, f); if (!file.exists(p)) { message("missing: ", p); return(NULL) }
  read.csv(p, stringsAsFactors=FALSE) %>% dplyr::filter(!is.na(call_cfg)) %>%
    dplyr::transmute(age, sex, source, call=dplyr::recode(call_cfg, buffered="B", amplified="A", `-`="0"))
}
kp <- c("sex","source")
pw <- function(d, cols) d %>% dplyr::select(dplyr::all_of(c(kp, "age", cols))) %>%
  tidyr::pivot_wider(names_from=age, values_from=dplyr::all_of(cols), names_sep="_")
PW  <- pw(pool, c("call","eB","eA","rank"))
PWm <- pw(pcfg("MethodCheck_AllTFs_Pooled_mlm_collectri.csv"), "call") %>% dplyr::rename(m_mlm=middle, o_mlm=old)
PWd <- pw(pcfg("MethodCheck_AllTFs_Pooled_ulm_dorothea.csv"), "call") %>% dplyr::rename(m_doro=middle, o_doro=old)
PX  <- PW %>% dplyr::left_join(PWm, by=kp) %>% dplyr::left_join(PWd, by=kp)
n_expr <- expr %>% dplyr::filter(expressed) %>% dplyr::distinct(sex, celltype, source) %>% dplyr::count(sex, source, name="n_ct_expressed")
n_ct_total <- expr %>% dplyr::distinct(sex, celltype) %>% dplyr::count(sex, name="n_ct_total")

sp <- dplyr::bind_rows(lapply(CATS, function(cat) {
  d <- PX[member(cat, PX$call_middle, PX$call_old) %in% TRUE, ]
  if (!nrow(d)) return(NULL)
  in_mlm <- member(cat, d$m_mlm, d$o_mlm); doro_scored <- !(is.na(d$m_doro) & is.na(d$o_doro)); in_doro <- member(cat, d$m_doro, d$o_doro)
  d %>% dplyr::mutate(category=cat, evidence=evidence(cat, d), flipped = cat == "Lost" & call_old == "A",
    tier = dplyr::case_when(in_mlm & doro_scored & in_doro ~ TIERS[1], in_mlm & !doro_scored ~ TIERS[2],
                            in_mlm ~ TIERS[3], doro_scored & in_doro ~ TIERS[4], TRUE ~ TIERS[5]))
})) %>%
  dplyr::left_join(n_expr, by=kp) %>% dplyr::left_join(n_ct_total, by="sex") %>%
  dplyr::mutate(n_ct_expressed=ifelse(source %in% tfs, dplyr::coalesce(n_ct_expressed, 0L), NA_integer_)) %>%
  dplyr::group_by(sex, category) %>% dplyr::arrange(dplyr::desc(evidence), .by_group=TRUE) %>%
  dplyr::mutate(pooled_rank_in_category=dplyr::row_number()) %>% dplyr::ungroup()
# add the per-cell-type support for the same TF and category
sp <- sp %>% dplyr::left_join(rank_tab %>% dplyr::select(sex, category, source, n_ct_called=n_ct, n_ct_robust, celltypes),
                              by=c("sex","category","source"))
write.csv(sp %>% dplyr::select(sex, category, pooled_rank_in_category, source, evidence, tier, flipped,
                               rank_middle, rank_old, call_middle, call_old, n_ct_expressed, n_ct_total,
                               n_ct_called, n_ct_robust, celltypes) %>%
            dplyr::arrange(sex, factor(category, levels=CATS), pooled_rank_in_category),
          file.path(out_dir, "TF_Atlas_Superpooled.csv"), row.names=FALSE)
cat("\nSuperpooled TFs per sex x category (all kept):\n")
print(as.data.frame(sp %>% dplyr::count(sex, category) %>% tidyr::pivot_wider(names_from=sex, values_from=n)), row.names=FALSE)

for (sx in c("female","male")) {
  d <- sp %>% dplyr::filter(sex == sx, pooled_rank_in_category <= 20) %>%
    dplyr::mutate(category=factor(category, levels=CATS), tier=factor(tier, levels=TIERS),
                  key=paste(category, source, sep="___"),
                  lab=sprintf("%s  (%s CT)", source, ifelse(is.na(n_ct_called), "0", n_ct_called)),
                  expr_lab=dplyr::case_when(is.na(n_ct_expressed) ~ "not a gene in data",
                                            n_ct_expressed == 0 ~ "not detected in any cell type", TRUE ~ "detected"))
  if (!nrow(d)) next
  d <- d %>% dplyr::arrange(category, evidence) %>% dplyr::mutate(key=factor(key, levels=unique(key)))
  p <- ggplot(d, aes(evidence, key, fill=tier)) +
    geom_col(width=0.72, aes(color=flipped), linewidth=0.35) +
    geom_point(data=d %>% dplyr::filter(expr_lab != "detected"), mapping=aes(x=-0.15, y=key, shape=expr_lab),
               size=1.3, inherit.aes=FALSE) +
    scale_y_discrete(labels=setNames(d$lab, d$key)) +
    scale_fill_manual(values=TIER_COL, name="Holds with other methods", drop=FALSE) +
    scale_color_manual(values=c(`FALSE`=NA, `TRUE`="#B2182B"), labels=c(`FALSE`="", `TRUE`="flipped to amplified"), name=NULL) +
    scale_shape_manual(values=c(`not detected in any cell type`=2, `not a gene in data`=0), name="TF mRNA") +
    facet_wrap(~category, scales="free_y", ncol=4) +
    theme_classic(base_size=7) +
    theme(axis.text.y=element_text(size=5, face="italic"), strip.background=element_rect(fill="grey92", color=NA),
          strip.text=element_text(face="bold"), legend.position="bottom", legend.key.size=unit(3, "mm")) +
    labs(title=sprintf("Superpooled TF categories, NE %s (top 20 per category; all TFs in TF_Atlas_Superpooled.csv)",
                       ifelse(sx == "female", "females", "males")),
         subtitle="Bar = evidence (-log10 p of the summed score vs relabelings); colour = method robustness; (n CT) = cell types where called in that category",
         x="Evidence (-log10 p)", y=NULL)
  ggsave(file.path(out_dir, paste0("TF_Atlas_Superpooled_", sx, ".pdf")), p, width=183 / 25.4, height=150 / 25.4)
  ggsave(file.path(out_dir, paste0("TF_Atlas_Superpooled_", sx, ".png")), p, width=183 / 25.4, height=150 / 25.4, dpi=400, bg="white")
}

################################################################################
### OVERVIEW FIGURES
###   O1  Where each category is: number of TFs per cell type x category (robust in brackets)
###   O2  How robust each category is: share of calls by method-robustness tier
###   O3  Middle vs old (superpooled): every TF's signed evidence at middle (x) and old (y) age;
###       quadrants = maintained / lost / age-emergent / amplified; top TFs labelled
###   O4  Females vs males: each TF's category score in females (x) and males (y);
###       shared vs sex-specific TFs; top TFs labelled
################################################################################

HAS_REPEL <- requireNamespace("ggrepel", quietly=TRUE)
txt_lab <- function(data, size=1.7) if (HAS_REPEL)
  ggrepel::geom_text_repel(data=data, aes(label=source), size=size, fontface="italic", segment.size=0.15,
                           max.overlaps=40, min.segment.length=0.1) else
  geom_text(data=data, aes(label=source), size=size, fontface="italic", vjust=-0.7)
thO <- theme_classic(base_size=7) +
  theme(strip.background=element_rect(fill="grey92", color=NA), strip.text=element_text(face="bold"),
        plot.title=element_text(face="bold", size=8), plot.subtitle=element_text(color="grey35", size=5.8),
        legend.key.size=unit(3, "mm"), legend.position="bottom")
svO <- function(p, nm, w, h) {
  ggsave(file.path(out_dir, paste0(nm, ".pdf")), p, width=w / 25.4, height=h / 25.4, limitsize=FALSE)
  ggsave(file.path(out_dir, paste0(nm, ".png")), p, width=w / 25.4, height=h / 25.4, dpi=400, bg="white", limitsize=FALSE)
}

# O1: TFs per cell type x category
o1 <- long %>% dplyr::mutate(robust=tier %in% TIERS[1:2]) %>% dplyr::group_by(sex, category, celltype) %>%
  dplyr::summarise(n=dplyr::n(), n_rob=sum(robust), .groups="drop") %>%
  dplyr::mutate(category=factor(category, levels=rev(CATS)), celltype=factor(celltype, levels=CT_ORDER),
                sex=ifelse(sex == "female", "Females", "Males"))
svO(ggplot(o1, aes(celltype, category, fill=n)) + geom_tile(color="white", linewidth=0.4) +
      geom_text(aes(label=sprintf("%d\n(%d)", n, n_rob)), size=1.6, lineheight=0.85) +
      scale_fill_gradient(low="#F7FBFF", high="#2171B5", name="TFs") + scale_x_discrete(labels=short_ct) +
      facet_wrap(~sex, ncol=1) + thO + theme(axis.text.x=element_text(angle=45, hjust=1), axis.line=element_blank()) +
      labs(title="Where TF buffering categories occur (NE)", subtitle="Number of TFs per cell type and category; (robust across methods)",
           x=NULL, y=NULL),
    "TF_Atlas_O1_Overview_CellType_x_Category", 183, 120)

# O2: robustness by category
o2 <- long %>% dplyr::count(sex, category, tier) %>% dplyr::group_by(sex, category) %>% dplyr::mutate(pct=100 * n / sum(n)) %>%
  dplyr::ungroup() %>% dplyr::mutate(tier=factor(tier, levels=rev(TIERS)), category=factor(category, levels=rev(CATS)),
                                     sex=ifelse(sex == "female", "Females", "Males"))
svO(ggplot(o2, aes(pct, category, fill=tier)) + geom_col(width=0.7) +
      scale_fill_manual(values=TIER_COL, name=NULL, breaks=TIERS) + facet_wrap(~sex) + thO +
      guides(fill=guide_legend(nrow=2)) +
      labs(title="How robust each category is to the TF scoring method",
           subtitle="Share of TF x cell-type calls that also hold with MLM + CollecTRI and ULM + DoRothEA", x="% of calls", y=NULL),
    "TF_Atlas_O2_Robustness", 150, 70)

# O3: middle vs old, superpooled, signed evidence
o3 <- PX %>% dplyr::mutate(
    x=ifelse(eB_middle >= eA_middle, eB_middle, -eA_middle), y=ifelse(eB_old >= eA_old, eB_old, -eA_old),
    traj=dplyr::case_when(call_middle == "B" & call_old == "B" ~ "Maintained",
                          call_middle == "B" & call_old == "A" ~ "Lost (flipped)",
                          call_middle == "B" ~ "Lost",
                          call_old == "B" ~ "Age-emergent",
                          call_old == "A" ~ "Amplified at old age",
                          call_middle == "A" ~ "Amplified at middle only", TRUE ~ "Not called"),
    robust_any=paste(sex, source) %in% paste(sp$sex[sp$tier %in% TIERS[1:2]], sp$source[sp$tier %in% TIERS[1:2]]),
    sex=ifelse(sex == "female", "Females", "Males"))
TRAJ_COL <- c(Maintained="#08519C", Lost="#969696", `Lost (flipped)`="#B2182B", `Age-emergent`="#4292C6",
              `Amplified at old age`="#EF6548", `Amplified at middle only`="#FDAE6B", `Not called`="grey88")
lab3 <- o3 %>% dplyr::filter(traj != "Not called", robust_any) %>% dplyr::group_by(sex) %>%
  dplyr::slice_max(abs(x) + abs(y), n=18, with_ties=FALSE) %>% dplyr::ungroup()
lim <- max(abs(c(o3$x, o3$y)), na.rm=TRUE)
svO(ggplot(o3, aes(x, y)) +
      geom_hline(yintercept=0, color="grey60", linewidth=0.3) + geom_vline(xintercept=0, color="grey60", linewidth=0.3) +
      geom_point(data=o3 %>% dplyr::filter(traj == "Not called"), color="grey88", size=0.4) +
      geom_point(data=o3 %>% dplyr::filter(traj != "Not called"), aes(color=traj), size=0.8, alpha=0.85) +
      txt_lab(lab3) + scale_color_manual(values=TRAJ_COL, name=NULL) +
      coord_equal(xlim=c(-lim, lim), ylim=c(-lim, lim)) + facet_wrap(~sex) + thO +
      guides(color=guide_legend(nrow=2, override.aes=list(size=2))) +
      annotate("text", x=lim, y=lim, label="buffered at both", hjust=1, vjust=1, size=1.8, color="grey40") +
      annotate("text", x=lim, y=-lim, label="lost / flipped", hjust=1, vjust=0, size=1.8, color="grey40") +
      annotate("text", x=-lim, y=lim, label="age-emergent", hjust=0, vjust=1, size=1.8, color="grey40") +
      annotate("text", x=-lim, y=-lim, label="amplified at both", hjust=0, vjust=0, size=1.8, color="grey40") +
      labs(title="Every TF, middle vs old age (superpooled, NE)",
           subtitle="Signed evidence: + = buffered, - = amplified (-log10 p vs relabelings); labels = robust TFs with the largest evidence",
           x="Middle age", y="Old age"),
    "TF_Atlas_O3_MiddleVsOld_Superpooled", 183, 105)

# O4: females vs males, per category (per-cell-type ranking scores)
# OFF by default: sexes were processed in separate batches, so values are not compared across sexes.
RUN_O4 <- FALSE
if (RUN_O4) {
o4 <-rank_tab %>% dplyr::select(sex, category, source, score, n_ct_robust) %>%
  tidyr::pivot_wider(names_from=sex, values_from=c(score, n_ct_robust), values_fill=0) %>%
  dplyr::mutate(category=factor(category, levels=CATS),
                where=dplyr::case_when(score_female > 0 & score_male > 0 ~ "both sexes",
                                       score_female > 0 ~ "females only", TRUE ~ "males only"),
                rob=(n_ct_robust_female + n_ct_robust_male) > 0)
lab4 <- o4 %>% dplyr::filter(rob) %>% dplyr::group_by(category) %>%
  dplyr::slice_max(score_female + score_male, n=6, with_ties=FALSE) %>% dplyr::ungroup()
svO(ggplot(o4, aes(score_female, score_male)) +
      geom_abline(slope=1, intercept=0, linetype="dashed", color="grey70", linewidth=0.3) +
      geom_point(aes(color=where), size=0.6, alpha=0.7) + txt_lab(lab4, 1.5) +
      scale_color_manual(values=c(`both sexes`="#6A3D9A", `females only`="#D35A3B", `males only`="#3A62C4"), name=NULL) +
      facet_wrap(~category, scales="free", ncol=4) + thO +
      guides(color=guide_legend(override.aes=list(size=2))) +
      labs(title="Females vs males: TF scores per category (NE)",
           subtitle="Summed evidence across cell types; points on an axis = one sex only; labels = top robust TFs per category",
           x="Score in females", y="Score in males"),
    "TF_Atlas_O4_FemalesVsMales", 183, 110)
}

cat("\nOutputs in:", out_dir, "\n")

invisible(NULL)
}

################################################################################
### STEP: fold  (= Buffering_TF_FoldConsistency.R)
################################################################################

STEPS$fold <- function() {
################################################################################
### Buffering_TF_FoldConsistency.R
###
### Is each TF result consistent across ALL THREE held-out control mice?
### (TF-level counterpart of the T1 gene tier.)
###
### The TF buffering score is S = mean over 3 held-out fits of -(a_aging x a_NE), one fit per held-out
### HC mouse k:  B_k = XHC(-k) vs YHC,  E_k = XNE vs XHC(k).  This script recomputes the three fits for the
### REAL labeling only and keeps each fit's score S_k, so every TF call can be annotated:
###   consistent  = buffered call with S_k > 0 in all 3 fits, or amplified call with S_k < 0 in all 3 fits
### At middle age this shows whether a TF result holds against every middle HC mouse, including the one
### sequenced with the NE mice (fold "co-sequenced"), or rests on the two separate-run mice.
### Same rule for both sexes and both ages. Nothing in the existing scores, calls or atlas is changed:
### the flag is an annotation (like the method-robustness tier and the expression check).
###
### Inputs (existing): Buffering_LOO_Split/pseudobulk_counts_and_meta.RDS, HeldOut_FitCache (if present),
###   Buffering_Candidates_Master/{Step1_CellTypes.csv, Step5_TFActivity_PerCellType.csv,
###   Step5_TFActivity_AcrossCellTypes.csv, TF_Atlas/TF_Atlas_Long.csv, TF_Atlas/TF_Atlas_Ranking.csv}
### Outputs: Buffering_Candidates_Master/TF_FoldConsistency/
################################################################################

Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({
  library(conflicted); library(DESeq2); library(dplyr); library(tidyr); library(parallel); library(decoupleR)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::mutate, .quiet=TRUE); conflicts_prefer(base::intersect, .quiet=TRUE)
  conflicts_prefer(base::setdiff, .quiet=TRUE); conflicts_prefer(base::unique, .quiet=TRUE)
  conflicts_prefer(base::Filter, .quiet=TRUE); conflicts_prefer(base::unname, .quiet=TRUE)
})
if (requireNamespace("RhpcBLASctl", quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(1)

################################################################################
### PARAMETERS (match Buffering_Candidates_Master.R)
################################################################################

base_path <- NE_BASE
CACHE_DIR <- file.path(base_path, "HeldOut_FitCache")
MASTER    <- MASTER_DIR
OUT       <- file.path(MASTER, "TF_FoldConsistency"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
N_CORES   <- 16
MIN_REGULON <- 5
P_AGING     <- 0.05             # per-fit aging threshold (as in Buffering_Candidates_Master.R)
P_NE        <- 0.05             # per-fit NE threshold (T1 rule)
RUN_GENES   <- TRUE             # gene-level check of T3/A3 calls (refits the held-out fits; cached afterwards)
SEPARATE_RUN <- "SR011142"      # run label in orig.ident of the separately sequenced middle HC mice
CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN","Astrocytes","Oligodendrocytes","OPC",
              "Microglia","Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell")
# which age(s) a TF atlas category rests on (the call that must be consistent)
CAT_AGES <- list(`Buffered, middle`="middle", `Buffered, old`="old", `Amplified, middle`="middle",
                 `Amplified, old`="old", Maintained=c("middle","old"), Lost="middle", `Age-emergent`="old")

################################################################################
### DATA AND FITS
################################################################################

pb <- readRDS(file.path(base_path, "Buffering_LOO_Split", "pseudobulk_counts_and_meta.RDS"))
pb_all <- pb$pb_all; pb_meta <- pb$pb_meta
age_vals <- unique(pb_meta$age_norm)
AGES <- c(middle=age_vals[!age_vals %in% c("young","old")][1], old="old")
short_id <- function(x) sub(".*--", "", x)
run_of   <- function(x) { r <- regmatches(x, regexpr("SR0[0-9]+", x)); if (length(r)) r else NA_character_ }

ct_status <- read.csv(file.path(MASTER, "Step1_CellTypes.csv"), stringsAsFactors=FALSE)
units <- tidyr::crossing(ct_status %>% dplyr::select(sex, celltype), age=c("middle","old"))
groups_of <- function(ag, sx, ct) {
  m <- pb_meta[pb_meta$celltype == ct & pb_meta$sex_norm == sx, ]
  list(m=m, ref=m$sample[m$age_norm == "young" & m$condition == "HC"],
       hc=m$sample[m$age_norm == AGES[[ag]] & m$condition == "HC"],
       ne=m$sample[m$age_norm == AGES[[ag]] & m$condition == "EE"])
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
  cnt  <- round(as.matrix(pb_all[, cols, drop=FALSE])); cnt <- cnt[rowSums(cnt >= 10) >= 2, , drop=FALSE]
  if (nrow(cnt) < 10) return(NULL)
  cd <- data.frame(grp=factor(rep(c("den","num"), c(length(den), length(num))), levels=c("den","num")), row.names=cols)
  tryCatch({
    dds <- suppressMessages(DESeqDataSetFromMatrix(cnt, cd, design=~grp))
    dds <- tryCatch(suppressMessages(DESeq(dds, quiet=TRUE)), error=function(e) {
      d <- estimateSizeFactors(dds); d <- estimateDispersionsGeneEst(d)
      dispersions(d) <- mcols(d)$dispGeneEst; nbinomWaldTest(d, quiet=TRUE) })
    r <- results(dds, name="grp_num_vs_den")
    data.frame(gene=rownames(r), lfc=r$log2FoldChange, stat=r$stat, pvalue=r$pvalue, stringsAsFactors=FALSE)
  }, error=function(e) NULL)
}
uniprot <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
net_signed <- readRDS(COLLECTRI_RDS) %>%
  dplyr::filter(!grepl(uniprot, source), !grepl(uniprot, target), target %in% rownames(pb_all)) %>%
  dplyr::distinct(source, target, .keep_all=TRUE)
wald <- function(f) if (is.null(f)) NULL else { s <- setNames(f$stat, f$gene); s[!is.na(s)] }
tf_activity <- function(B, E) {
  if (is.null(B) || is.null(E)) return(NULL)
  g <- intersect(names(B), names(E)); if (length(g) < 500) return(NULL)
  mat <- cbind(B=B[g], E=E[g]); rownames(mat) <- g
  r <- tryCatch(decoupleR::run_ulm(mat=mat, net=net_signed, .source="source", .target="target", .mor="mor",
                                   minsize=MIN_REGULON), error=function(e) NULL)
  if (is.null(r)) return(NULL)
  r %>% dplyr::filter(statistic == "ulm") %>% dplyr::select(source, condition, score) %>%
    tidyr::pivot_wider(names_from=condition, values_from=score, names_prefix="a_")
}

################################################################################
### PER-FOLD TF SCORES (real labeling only)
################################################################################

f_folds <- file.path(OUT, "TF_PerFold_Scores.csv")
folds <- NULL
if (file.exists(f_folds)) {
  folds <- tryCatch(read.csv(f_folds, stringsAsFactors=FALSE), error=function(e) NULL)
  if (is.null(folds) || !nrow(folds)) { cat("Existing", basename(f_folds), "is empty (earlier run failed); recomputing\n")
    unlink(f_folds); folds <- NULL } else cat("Loaded:", f_folds, "\n")
}
if (is.null(folds)) {
  one_unit <- function(i) {
    ag <- units$age[i]; sx <- units$sex[i]; ct <- units$celltype[i]; g <- groups_of(ag, sx, ct)
    if (length(g$ref) < 2 || length(g$hc) < 3 || length(g$ne) < 2) return(NULL)
    oid <- g$m$orig.ident[match(g$hc, g$m$sample)]
    dplyr::bind_rows(lapply(seq_along(g$hc), function(k) {
      a <- tf_activity(wald(fit_one(g$hc[-k], g$ref)), wald(fit_one(g$ne, g$hc[k]))); if (is.null(a)) return(NULL)
      a %>% dplyr::transmute(age=ag, sex=sx, celltype=ct, source, heldout_HC=short_id(oid[k]), run=run_of(oid[k]),
                             a_aging=a_B, a_NE=a_E, S_k=-a_B * a_E)
    }))
  }
  cat(sprintf("Scoring %d units (3 held-out fits each, real labeling) on %d cores...\n", nrow(units), N_CORES))
  res <- mclapply(seq_len(nrow(units)), function(i) tryCatch(one_unit(i), error=function(e) conditionMessage(e)), mc.cores=N_CORES)
  errs <- unlist(Filter(is.character, res)); if (length(errs)) { cat("Errors:\n"); print(head(errs)) }
  folds <- dplyr::bind_rows(Filter(is.data.frame, res))
  n_null <- sum(vapply(res, is.null, TRUE))
  cat(sprintf("Units scored: %d | skipped (too few mice / no fit): %d | errors: %d\n",
              length(Filter(is.data.frame, res)), n_null, length(errs)))
  if (!nrow(folds)) stop("No TF fold scores were produced; see the errors above (nothing written)")
  write.csv(folds, f_folds, row.names=FALSE)
}
folds <- folds %>% dplyr::mutate(fold=dplyr::case_when(age == "middle" & run == SEPARATE_RUN ~ "separate-run HC",
                                                       age == "middle" ~ "co-sequenced HC", TRUE ~ "single-run HC"))

################################################################################
### FLAG PER TF x CELL TYPE
################################################################################

act <- read.csv(file.path(MASTER, "Step5_TFActivity_PerCellType.csv"), stringsAsFactors=FALSE) %>%
  dplyr::mutate(call=dplyr::case_when(buf ~ "buffered", amp ~ "amplified", TRUE ~ "-")) %>%
  dplyr::select(age, sex, celltype, source, S, call)
cons <- folds %>% dplyr::group_by(age, sex, celltype, source) %>%
  dplyr::summarise(n_folds=dplyr::n(), n_pos=sum(S_k > 0), n_neg=sum(S_k < 0), S_mean_folds=mean(S_k),
                   S_min=min(S_k), S_max=max(S_k),
                   S_cosequenced=if (any(fold == "co-sequenced HC")) S_k[fold == "co-sequenced HC"][1] else NA_real_,
                   S_separate_run_mean=if (any(fold == "separate-run HC")) mean(S_k[fold == "separate-run HC"]) else NA_real_,
                   .groups="drop") %>%
  dplyr::left_join(act, by=c("age","sex","celltype","source")) %>%
  dplyr::mutate(call=dplyr::coalesce(call, "-"),
                consistent=dplyr::case_when(call == "buffered" ~ n_folds == 3 & n_pos == 3,
                                            call == "amplified" ~ n_folds == 3 & n_neg == 3, TRUE ~ NA),
                same_sign_all_folds=n_folds == 3 & (n_pos == 3 | n_neg == 3))
write.csv(cons %>% dplyr::arrange(sex, age, factor(celltype, levels=CT_ORDER), dplyr::desc(S)),
          file.path(OUT, "TF_FoldConsistency_PerCellType.csv"), row.names=FALSE)
chk <- suppressWarnings(cor(cons$S, cons$S_mean_folds, use="complete.obs"))
cat(sprintf("\nSanity: correlation of recomputed mean score with Step 5 score = %.3f (should be ~1)\n", chk))

cat("\n=== Calls consistent across all three held-out control mice ===\n")
print(as.data.frame(cons %>% dplyr::filter(call != "-") %>% dplyr::group_by(sex, age, call) %>%
  dplyr::summarise(n_calls=dplyr::n(), n_consistent=sum(consistent %in% TRUE),
                   pct_consistent=round(100 * mean(consistent %in% TRUE), 1), .groups="drop") %>%
  dplyr::arrange(sex, age, call)), row.names=FALSE)
cat("\n--- Middle age: calls by sign of the co-sequenced-HC fold ---\n")
print(as.data.frame(cons %>% dplyr::filter(call != "-", age == "middle") %>% dplyr::group_by(sex, call) %>%
  dplyr::summarise(n_calls=dplyr::n(),
                   cosequenced_agrees=sum((call == "buffered" & S_cosequenced > 0) | (call == "amplified" & S_cosequenced < 0), na.rm=TRUE),
                   .groups="drop")), row.names=FALSE)

################################################################################
### SUPERPOOLED (sum over cell types per held-out mouse; the same mouse is held out in every cell type)
################################################################################

pool_f <- folds %>% dplyr::group_by(age, sex, heldout_HC, fold, source) %>%
  dplyr::summarise(S_sum_k=sum(S_k), n_ct=dplyr::n(), .groups="drop") %>% dplyr::filter(n_ct >= 2)
pool <- read.csv(file.path(MASTER, "Step5_TFActivity_AcrossCellTypes.csv"), stringsAsFactors=FALSE) %>%
  dplyr::mutate(call=dplyr::case_when(buf ~ "buffered", amp ~ "amplified", TRUE ~ "-")) %>% dplyr::select(age, sex, source, S_sum, call)
pool_c <- pool_f %>% dplyr::group_by(age, sex, source) %>%
  dplyr::summarise(n_folds=dplyr::n(), n_pos=sum(S_sum_k > 0), n_neg=sum(S_sum_k < 0),
                   S_sum_cosequenced=if (any(fold == "co-sequenced HC")) S_sum_k[fold == "co-sequenced HC"][1] else NA_real_, .groups="drop") %>%
  dplyr::left_join(pool, by=c("age","sex","source")) %>%
  dplyr::mutate(call=dplyr::coalesce(call, "-"),
                consistent=dplyr::case_when(call == "buffered" ~ n_folds == 3 & n_pos == 3,
                                            call == "amplified" ~ n_folds == 3 & n_neg == 3, TRUE ~ NA))
write.csv(pool_c, file.path(OUT, "TF_FoldConsistency_Superpooled.csv"), row.names=FALSE)
cat("\n=== Superpooled calls consistent across all three held-out control mice ===\n")
print(as.data.frame(pool_c %>% dplyr::filter(call != "-") %>% dplyr::group_by(sex, age, call) %>%
  dplyr::summarise(n_calls=dplyr::n(), n_consistent=sum(consistent %in% TRUE), .groups="drop")), row.names=FALSE)

################################################################################
### ANNOTATE THE TF ATLAS AND RANKING (no rows removed, no values changed)
################################################################################

long <- read.csv(file.path(MASTER, "TF_Atlas", "TF_Atlas_Long.csv"), stringsAsFactors=FALSE)
cw <- cons %>% dplyr::select(age, sex, celltype, source, consistent) %>%
  tidyr::pivot_wider(names_from=age, values_from=consistent, names_prefix="consistent_")
long_c <- long %>% dplyr::left_join(cw, by=c("sex","celltype","source")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(fold_consistent=all(c(middle=consistent_middle %in% TRUE, old=consistent_old %in% TRUE)[CAT_AGES[[category]]])) %>%
  dplyr::ungroup()
write.csv(long_c, file.path(OUT, "TF_Atlas_Long_withFoldConsistency.csv"), row.names=FALSE)
rk <- read.csv(file.path(MASTER, "TF_Atlas", "TF_Atlas_Ranking.csv"), stringsAsFactors=FALSE)
rk_c <- rk %>% dplyr::left_join(long_c %>% dplyr::group_by(sex, category, source) %>%
  dplyr::summarise(n_ct_fold_consistent=sum(fold_consistent),
                   celltypes_fold_consistent=paste(celltype[fold_consistent], collapse=", "), .groups="drop"),
  by=c("sex","category","source"))
write.csv(rk_c, file.path(OUT, "TF_Atlas_Ranking_withFoldConsistency.csv"), row.names=FALSE)
cat("\n=== Atlas: TF x cell-type entries consistent across all three control mice, per category ===\n")
print(as.data.frame(long_c %>% dplyr::group_by(sex, category) %>%
  dplyr::summarise(n_entries=dplyr::n(), n_fold_consistent=sum(fold_consistent),
                   pct=round(100 * mean(fold_consistent), 1), .groups="drop") %>%
  dplyr::arrange(sex, factor(category, levels=names(CAT_AGES)))), row.names=FALSE)
################################################################################
### GENES: are T3 / A3 calls supported by the co-sequenced control mouse?
###   Per held-out fit k, a gene is "opposed" if aging p < P_AGING and NE runs opposite to aging (the
###   per-fit rule behind T1/T3), and "opposed (direction)" if NE simply runs opposite to aging.
###   T3 = opposed in >= 2 of 3 fits, so at middle age it can rest on the two separate-run fits alone.
###   For every T3 (and A3) gene this reports whether the CO-SEQUENCED fold also supports the call.
###   Old age (all controls from one run) gives the reference: how often a single fold supports a T3
###   call when no fold is special.
################################################################################

if (RUN_GENES) {
  f_gene <- file.path(OUT, "Gene_FoldSupport_T3A3_v2.csv")   # v2: adds NE p per fit
  gene <- if (file.exists(f_gene) && file.size(f_gene) > 100) read.csv(f_gene, stringsAsFactors=FALSE) else NULL
  if (is.null(gene)) {
    gene_unit <- function(i) {
      ag <- units$age[i]; sx <- units$sex[i]; ct <- units$celltype[i]; g <- groups_of(ag, sx, ct)
      if (length(g$ref) < 2 || length(g$hc) < 3 || length(g$ne) < 2) return(NULL)
      oid <- g$m$orig.ident[match(g$hc, g$m$sample)]
      pf <- dplyr::bind_rows(lapply(seq_along(g$hc), function(k) {
        B <- fit_one(g$hc[-k], g$ref); E <- fit_one(g$ne, g$hc[k]); if (is.null(B) || is.null(E)) return(NULL)
        dplyr::inner_join(B %>% dplyr::select(gene, lfc_B=lfc, p_B=pvalue), E %>% dplyr::select(gene, lfc_E=lfc, p_E=pvalue), by="gene") %>%
          dplyr::filter(!is.na(lfc_B), !is.na(lfc_E), lfc_B != 0) %>%
          dplyr::mutate(fold=ifelse(ag == "middle" & run_of(oid[k]) == SEPARATE_RUN, "separate-run HC",
                                    ifelse(ag == "middle", "co-sequenced HC", paste0("fold", k))),
                        ag_sig=!is.na(p_B) & p_B < P_AGING, opp_dir=sign(lfc_E) != sign(lfc_B),
                        opp=ag_sig & opp_dir, same=ag_sig & !opp_dir, opp_nesig=opp & !is.na(p_E) & p_E < P_NE)
      }))
      if (!nrow(pf)) return(NULL)
      pf %>% dplyr::group_by(gene) %>%
        dplyr::summarise(n_fits=dplyr::n(), n_aging=sum(ag_sig), n_opp=sum(opp), n_same=sum(same), n_opp_dir=sum(opp_dir),
                         coseq_opp=if (any(fold == "co-sequenced HC")) opp[fold == "co-sequenced HC"][1] else NA,
                         coseq_opp_dir=if (any(fold == "co-sequenced HC")) opp_dir[fold == "co-sequenced HC"][1] else NA,
                         coseq_same=if (any(fold == "co-sequenced HC")) same[fold == "co-sequenced HC"][1] else NA,
                         fold1_opp_dir=if (ag == "old") opp_dir[fold == "fold1"][1] else NA,
                         n_opp_nesig=sum(opp_nesig),
                         coseq_opp_nesig=if (any(fold == "co-sequenced HC")) opp_nesig[fold == "co-sequenced HC"][1] else NA,
                         seprun_opp_nesig_mean=if (any(fold == "separate-run HC")) mean(opp_nesig[fold == "separate-run HC"]) else NA_real_,
                         fold_opp_nesig_mean=mean(opp_nesig),
                         .groups="drop") %>%
        dplyr::mutate(T3=n_aging >= 2 & n_opp >= 2, A3=n_aging >= 2 & n_same >= 2) %>%
        dplyr::filter(T3 | A3) %>% dplyr::mutate(age=ag, sex=sx, celltype=ct, .before=1)
    }
    cat(sprintf("\nGenes: refitting %d units (3 held-out fits each) on %d cores...\n", nrow(units), N_CORES))
    rg <- mclapply(seq_len(nrow(units)), function(i) tryCatch(gene_unit(i), error=function(e) conditionMessage(e)), mc.cores=N_CORES)
    ge <- unlist(Filter(is.character, rg)); if (length(ge)) { cat("Errors:\n"); print(head(ge)) }
    gene <- dplyr::bind_rows(Filter(is.data.frame, rg))
    if (!nrow(gene)) stop("No gene fold results produced; see errors above")
    write.csv(gene, f_gene, row.names=FALSE)
  }
  # sanity: recomputed T3 counts vs Step 2
  s2 <- readRDS(file.path(MASTER, "Step2_GeneCalls.RDS")) %>% dplyr::group_by(age, sex, celltype) %>%
    dplyr::summarise(T3_step2=sum(T3), .groups="drop")
  sc <- gene %>% dplyr::group_by(age, sex, celltype) %>% dplyr::summarise(T3_here=sum(T3), .groups="drop") %>%
    dplyr::inner_join(s2, by=c("age","sex","celltype"))
  cat(sprintf("\nSanity: T3 counts recomputed vs Step 2: correlation %.3f, total %d vs %d\n",
              suppressWarnings(cor(sc$T3_here, sc$T3_step2)), sum(sc$T3_here), sum(sc$T3_step2)))

  cat("\n=== Middle age: T3 (buffered) and A3 (amplified) genes supported by the co-sequenced control fold ===\n")
  cat("   strict    = that fold alone meets the per-fit rule (aging p < 0.05 and NE opposite / same)\n")
  cat("   direction = in that fold NE simply runs opposite to aging (T3) or with it (A3)\n")
  mid <- gene %>% dplyr::filter(age == "middle")
  print(as.data.frame(dplyr::bind_rows(
    mid %>% dplyr::filter(T3) %>% dplyr::group_by(sex) %>% dplyr::summarise(set="T3 buffered", n_genes=dplyr::n(),
      pct_strict=round(100 * mean(coseq_opp %in% TRUE), 1), pct_direction=round(100 * mean(coseq_opp_dir %in% TRUE), 1),
      pct_all3_direction=round(100 * mean(n_opp_dir == 3), 1), .groups="drop"),
    mid %>% dplyr::filter(A3) %>% dplyr::group_by(sex) %>% dplyr::summarise(set="A3 amplified", n_genes=dplyr::n(),
      pct_strict=round(100 * mean(coseq_same %in% TRUE), 1), pct_direction=round(100 * mean(coseq_opp_dir %in% FALSE), 1),
      pct_all3_direction=round(100 * mean(n_opp_dir == 0), 1), .groups="drop"))), row.names=FALSE)

  cat("\n=== Old age reference (one run): T3 genes supported by a single fold (fold 1) / in all 3 folds, direction ===\n")
  print(as.data.frame(gene %>% dplyr::filter(age == "old", T3) %>% dplyr::group_by(sex) %>%
    dplyr::summarise(n_genes=dplyr::n(), pct_fold1_direction=round(100 * mean(fold1_opp_dir %in% TRUE), 1),
                     pct_all3_direction=round(100 * mean(n_opp_dir == 3), 1), .groups="drop")), row.names=FALSE)

  cat("\n--- Middle age per cell type: % of T3 genes whose co-sequenced fold runs opposite to aging ---\n")
  print(as.data.frame(mid %>% dplyr::filter(T3) %>% dplyr::group_by(sex, celltype) %>%
    dplyr::summarise(n_T3=dplyr::n(), pct_coseq_direction=round(100 * mean(coseq_opp_dir %in% TRUE), 1), .groups="drop") %>%
    dplyr::arrange(sex, factor(celltype, levels=CT_ORDER))), row.names=FALSE)

  cat("\n=== NE strength: buffered genes (aging p<0.05 and NE opposite in all 3 fits) whose NE effect is significant ===\n")
  cat("   T1 (pipeline)        = NE p<0.05 (opposite) in >= 2 of 3 fits -> at middle age the 2 separate-run fits can supply this\n")
  cat("   NE-sig, co-sequenced = NE p<0.05 (opposite) in the co-sequenced fold (middle age; free of the separate-run mice)\n")
  cat("   NE-sig, single fold  = same, expected for one fold (mean over folds); at old age this is the comparable reference\n")
  st <- gene %>% dplyr::mutate(all3=n_aging == 3 & n_opp == 3) %>% dplyr::group_by(sex, age) %>%
    dplyr::summarise(T3=sum(T3), opposed_all3=sum(all3), T1_pipeline=sum(all3 & n_opp_nesig >= 2),
                     NEsig_cosequenced=sum(all3 & coseq_opp_nesig %in% TRUE),
                     NEsig_separate_run_fold=round(sum(all3 * dplyr::coalesce(seprun_opp_nesig_mean, NA_real_)), 0),
                     NEsig_single_fold=round(sum(all3 * fold_opp_nesig_mean), 0), .groups="drop") %>%
    dplyr::mutate(NEsig_cosequenced=ifelse(age == "middle", NEsig_cosequenced, NA))
  print(as.data.frame(st %>% dplyr::arrange(sex, age)), row.names=FALSE)
  cat("   Compare middle 'NEsig_cosequenced' with old 'NEsig_single_fold' (both = one control fold, no separate-run mouse).\n")
  st_ct <- gene %>% dplyr::mutate(all3=n_aging == 3 & n_opp == 3) %>% dplyr::group_by(sex, celltype, age) %>%
    dplyr::summarise(clean=ifelse(age[1] == "middle", sum(all3 & coseq_opp_nesig %in% TRUE), round(sum(all3 * fold_opp_nesig_mean), 0)),
                     T1_pipeline=sum(all3 & n_opp_nesig >= 2), .groups="drop") %>%
    tidyr::pivot_wider(names_from=age, values_from=c(clean, T1_pipeline))
  cat("\n--- Per cell type: NE-significant buffered genes, one clean control fold (middle: co-sequenced; old: expected single fold) ---\n")
  print(as.data.frame(st_ct %>% dplyr::arrange(sex, factor(celltype, levels=CT_ORDER))), row.names=FALSE)
  write.csv(st_ct, file.path(OUT, "Gene_NEStrength_CleanFold_PerCellType.csv"), row.names=FALSE)
}

cat("\nOutputs in:", OUT, "\n")

invisible(NULL)
}

################################################################################
### STEP: figure5  (= Figure5_Buffering.R)
################################################################################

STEPS$figure5 <- function() {
################################################################################
### Figure5_Buffering.R
###
### Figure 5: NE buffers the aging hippocampal transcriptome (cell types and genes).
### Reads existing outputs only (no DESeq2 refits). One PDF + PNG per panel, plus
### source-data CSVs. Each sex is its own panel with its own scale (separate batches);
### sexes are shown side by side, never on a shared axis.
###
###   5a  Scheme: drawn by hand, not generated here
###   5b  Cosine arc diagrams (position = held-out cosine between the aging and NE vectors, scale -0.5 to +0.5, ARC_LIM)
###       b1  NE, males | females: one ring per cell type, dot = middle age, arrowhead = old age
###       b2  DREADD, old females (Approach B | D): one ring per cell type, filled diamond = DREADD,
###           open circle = NE old females on the same ring (reference)
###       b3  b1 + b2 side by side
###       b4  spoke arc diagrams (Cos_ArcDiagram.R style): one spoke per cell type at its held-out angle,
###           NE per age (sexes side by side and separately), DREADD B | D; faded = fails the gate
###   5c  Gene-set heatmaps per cell type (z-scored expression; columns = mice:
###       Young HC | Middle HC | Middle NE | Old HC | Old NE). Row blocks, each gene in one block:
###         Buffered, both ages (maintained) | Buffered at middle, lost at old | Buffered, old only |
###         Amplified, middle | Amplified, old      (split into up / down with age)
###       c1  with gene names: LAB_TIER genes (T1 only; blocks with < LAB_MIN T1 genes are not shown),
###           top LAB_N per block by aging change
###       c2  without gene names: all SET_TIER (T3) genes, up to UNLAB_MAX per block
###       c3  compact summary, all cell types: mean aging-oriented z per block and mouse group
###           (z multiplied by the sign of the gene's aging change, so higher = more old-like)
###       c4  fold-change strips, one file per cell type x gene set (buffered middle / old, lost, maintained,
###           amplified middle / old), labeled (LAB_TIER, top STRIP_N per direction by buffering strength) and
###           unlabeled (all SET_TIER): aging and NE log2 FC (mean of the 3 held-out fits). Single-age sets show
###           only their own age; lost and maintained show middle and old side by side.
### Re-running skips any panel whose PDF already exists in OUT (OVERWRITE <- FALSE); GO and DREADD results are reused.
###   5d  Pathway enrichment (GO BP, clusterProfiler) of each gene set (SET_TIER):
###         Buffered, middle | Buffered, old | Lost | Maintained | Amplified, middle | Amplified, old
###       per cell type, per cell class (neurons / glia / vascular) and superpooled (all cell types),
###       each sex separately. Universe = genes that could be in the set (same rule as Step 3):
###         buffered / amplified at an age: genes with an aging change at that age
###         lost / maintained: genes buffered at middle age and tested at old age
###       d1 superpooled dot plot   d2 per cell type heatmap (one per set)   d3 per cell class dot plot
###   h   ED10b: held-out cosine per cell type, middle and old age (females | males), gate solid / faded
###   5e  Quadrant: held-out cosine (x, -1..1) vs buffering index (y, raw, linear from 0), NE per sex x age;
###       DREADD (old females, B | D) alongside with NE old females faded
###       BI = sum(OEEYHC log2FC^2) / sum(OHCYHC log2FC^2) from the DESeq2 comparisons (step bi); < 1 = NE closer to young
###
### Gene sets follow Buffering_Candidates_Master.R Step 3 (buffered sets restricted to cell types
### that pass the gate at that age; lost / maintained need the gene tested at both ages).
###
### Inputs (existing)
###   Buffering_LOO_Split/pseudobulk_counts_and_meta.RDS
###   Buffering_MiddleAge_Robustness/*_Check1_BiasCorrected.csv, *_Check2_HeldOut_PerMouse.csv
###   Buffering_Candidates_Master/Step1_CellTypes.csv, Step2_GeneCalls.RDS
###   DREADD: DREADD_Buffering_HeldOut/*_DREADD_HeldOut_Summary.csv,
###           pseudobulk_mScarlet_50_Approach_{B,D}.RDS (only if the DREADD metrics source data
###           from Buffering_Paper_Figures.R is not found)
### Output: Paper_Figures/Figure5/<date>/
################################################################################

suppressPackageStartupMessages({
  library(conflicted); library(DESeq2); library(dplyr); library(tidyr); library(ggplot2); library(parallel)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::mutate, .quiet=TRUE); conflicts_prefer(dplyr::slice_max, .quiet=TRUE)
  conflicts_prefer(base::intersect, .quiet=TRUE); conflicts_prefer(base::setdiff, .quiet=TRUE)
  conflicts_prefer(base::union, .quiet=TRUE); conflicts_prefer(base::unname, .quiet=TRUE)
})
HAS_REPEL <- requireNamespace("ggrepel", quietly=TRUE)
HAS_PATCH <- requireNamespace("patchwork", quietly=TRUE)

################################################################################
### PARAMETERS
################################################################################

MASTER <- MASTER_DIR
FIG_DATE <- format(Sys.Date(), "%Y-%m-%d")
OUT    <- file.path(NE_BASE, "Paper_Figures", "Figure5", FIG_DATE)
OVERWRITE <- FALSE                            # FALSE = skip any panel whose PDF already exists (GO and DREADD results are reused too)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
# 5b arc / spoke scale (cosine -ARC_LIM .. +ARC_LIM). If the panels in OUT were drawn on another scale (or before the scale
# was recorded), the old Fig5b1-b4 files are removed once so they are redrawn; every other panel is still skipped.
ARC_LIM <- 0.5
arc_stamp <- file.path(OUT, "Fig5b_arc_scale.txt")
if (!file.exists(arc_stamp) || !identical(readLines(arc_stamp, warn=FALSE)[1], as.character(ARC_LIM))) {
  old_b <- list.files(OUT, pattern="^Fig5b[1-4]_", full.names=TRUE)
  if (length(old_b)) cat("Arc diagrams were drawn on another scale: removing", length(old_b), "old Fig5b1-b4 files to redraw them\n")
  file.remove(old_b); writeLines(as.character(ARC_LIM), arc_stamp)
}

PANELS      <- c("b","c","e","h")     # arcs, heatmaps + fold-change strips, quadrants, ED10b; add "d" for GO enrichment
MIN_GENES   <- 1000                   # cell types with fewer genes tested are left out (5b, 5e)
DR_APPROACHES <- c("B","D")
N_SIM       <- 50                     # DREADD simulated null draws (only if no saved DREADD metrics)
SET_TIER    <- "T3"                   # 5c unlabeled + summary, 5d: buffered tier (amplified = A3)
LAB_TIER    <- "T1"                   # 5c labeled heatmap: named genes from this tier ...
LAB_MIN     <- 5                      #   ... blocks / strips with fewer T1 genes are not shown (no fallback to SET_TIER)
LAB_N       <- 12                     # 5c labeled heatmap: genes per block, top by buffering strength (see RANK below)
UNLAB_MAX   <- 150                    # 5c unlabeled: max genes per block (c2) / per direction (c4)
STRIP_N     <- 15                     # 5c fold-change strips, labeled: genes per direction (up / down with age)
# RANK: named genes are the strongest buffered (or amplified) genes, scored per gene as
#   buffered:  -aging LFC x NE LFC  (large when both changes are big and opposite)
#   amplified: +aging LFC x NE LFC  (large when both are big and in the same direction)
#   maintained: the smaller of the middle and old scores (strong at both ages)
MID_COSEQ   <- TRUE                   # named middle-age buffered genes (buffered middle, lost, maintained): keep genes also
                                      #   NE-significant against the co-sequenced HC, when >= LAB_MIN such genes exist
FOLD_FILE   <- "TF_FoldConsistency/Gene_FoldSupport_T3A3_v2.csv"   # under MASTER (from Buffering_TF_FoldConsistency.R)
BI_SOURCE   <- "original"             # 5e: BI from the DESeq2 tables (step bi, *_BI_Summary.csv)
BI_LABEL_MAP <- NULL                  # if the label -> age mapping printed in 5e is wrong, e.g. c(old="old", middle="middle")
STRIP_LIM   <- 1.5                    # 5c strips: colour scale clipped at +/- this log2 FC
GO_ONT      <- "BP"
GO_MIN_SET  <- 5                      # set needs >= this many genes to be tested (lab pathway settings)
GO_FDR      <- 0.05
GO_TOP      <- 5                      # terms shown per set in d1 / d3, and per set in d2 (from any cell type)
GO_SIMPLIFY <- TRUE                   # collapse redundant GO terms (clusterProfiler::simplify, needs GOSemSim)
N_CORES     <- 8
CLASSES <- list(Neurons=c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN"),
                Glia=c("Astrocytes","Oligodendrocytes","OPC","Microglia"),
                Vascular=c("Endothelial","Pericyte","APC","Stromal"))

SEX_LAB <- c(male="Males", female="Females"); SEX_ORDER <- c("Males","Females")
AGE_LAB <- c(middle="Middle age", old="Old age")
CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN","Astrocytes","Oligodendrocytes","OPC",
              "Microglia","Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell")
CT_SHORT <- c("DG ExNs"="DG","CA3 ExNs"="CA3","CA1 ExNs"="CA1","CA2 ExNs"="CA2","Sub ExNs"="Sub","Astrocytes"="Astro",
              "Oligodendrocytes"="Oligo","Cajal-Retzius Cell"="CR","Endothelial"="Endo","Microglia"="Micro","Pericyte"="Peri")
short_ct <- function(x) ifelse(x %in% names(CT_SHORT), CT_SHORT[x], x)
CT_COLS <- c("DG ExNs"="#E53935","CA1 ExNs"="#CC0000","CA2 ExNs"="#8E24AA","CA3 ExNs"="#5E35B1","Sub ExNs"="#80D8FF",
             "IN"="#00897B","Astrocytes"="#43A047","Oligodendrocytes"="#1E88E5","OPC"="#00ACC1","Microglia"="#7CB342",
             "Endothelial"="#FB8C00","Pericyte"="#6D4C41","APC"="#546E7A","Stromal"="#27408B","Cajal-Retzius Cell"="#FDD835")
BLUE <- "#2166AC"; RED <- "#B2182B"
COL_AP <- c(B="#D6604D", D="#4393C3")
SETS  <- c("Buffered, middle","Buffered, old","Lost","Maintained","Amplified, middle","Amplified, old")
SET_AGE <- c(`Buffered, middle`="middle", `Buffered, old`="old", Lost="middle", Maintained="middle",
             `Amplified, middle`="middle", `Amplified, old`="old")   # age whose aging change sorts the rows
SET_COL <- c(`Buffered, middle`="#81C784", `Buffered, old`="#1B5E20", Lost="#C62828", Maintained="#2E7D32",
             `Amplified, middle`="#CE93D8", `Amplified, old`="#6A1B9A")
BLOCKS <- c("Buffered, both ages (maintained)","Buffered at middle, lost at old","Buffered, old only",
            "Amplified, middle","Amplified, old")
BLOCK_SHORT <- setNames(c("Maintained","Lost at old","Buffered old only","Amplified middle","Amplified old"), BLOCKS)
GROUPS <- c("Young HC","Middle HC","Middle NE","Old HC","Old NE")

th <- theme_classic(base_size=7) +
  theme(strip.background=element_rect(fill="grey93", color=NA), strip.text=element_text(face="bold", size=6.5),
        plot.title=element_text(face="bold", size=8), plot.subtitle=element_text(color="grey35", size=5.8),
        axis.line=element_line(linewidth=0.3), axis.ticks=element_line(linewidth=0.3),
        legend.key.size=unit(3, "mm"), legend.text=element_text(size=5.8), legend.title=element_text(size=6))
mm <- function(x) x / 25.4
index <- list()
done <- function(name, dir=OUT) !OVERWRITE && file.exists(file.path(dir, paste0(name, ".pdf")))
sv <- function(p, name, w, h, dir=OUT, note="") {
  if (done(name, dir)) { cat("  exists, skipped:", name, "\n"); return(invisible()) }
  h <- min(h, 1000)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width=mm(w), height=mm(h), limitsize=FALSE)
  ggsave(file.path(dir, paste0(name, ".png")), p, width=mm(w), height=mm(h), dpi=500, bg="white", limitsize=FALSE)
  index[[length(index) + 1]] <<- data.frame(file=file.path(basename(dir), name), width_mm=w, height_mm=h, note=note)
  cat("  saved:", name, "\n")
}
ann_th <- theme(plot.title=element_text(face="bold", size=8), plot.subtitle=element_text(color="grey35", size=5.8))
wsrc <- function(x, name) write.csv(x, file.path(OUT, paste0("SourceData_", name, ".csv")), row.names=FALSE)
latest <- function(dir, pat) { f <- list.files(dir, pattern=pat, recursive=TRUE, full.names=TRUE)
  if (!length(f)) stop("No file matching '", pat, "' under ", dir); f[order(file.mtime(f), decreasing=TRUE)][1] }
sexf <- function(x) factor(SEX_LAB[x], levels=SEX_ORDER)
tagify <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
side_by_side <- function(pl, guides="collect") {
  pl <- Filter(Negate(is.null), pl)
  if (length(pl) == 1 || !HAS_PATCH) return(pl[[1]])
  patchwork::wrap_plots(pl, nrow=1) + patchwork::plot_layout(guides=guides)
}

################################################################################
### SHARED INPUTS
################################################################################

cat("Loading inputs\n")
pb <- readRDS(file.path(NE_BASE, "Buffering_LOO_Split", "pseudobulk_counts_and_meta.RDS"))
pb_all <- pb$pb_all; pb_meta <- pb$pb_meta; rm(pb)
age_vals <- unique(pb_meta$age_norm); MID <- age_vals[!age_vals %in% c("young","old")][1]

rob <- ROB_ROOT
c1 <- read.csv(latest(rob, "_Check1_BiasCorrected\\.csv$"), stringsAsFactors=FALSE) %>% dplyr::mutate(age=tolower(age))
c2 <- read.csv(latest(rob, "_Check2_HeldOut_PerMouse\\.csv$"), stringsAsFactors=FALSE) %>% dplyr::mutate(age=tolower(age))
cosd <- c2 %>% dplyr::group_by(age, sex, celltype) %>% dplyr::summarise(cos=mean(cosine, na.rm=TRUE), .groups="drop") %>%
  dplyr::inner_join(c1 %>% dplyr::select(age, sex, celltype, n_genes, cos_dc, beyond_null), by=c("age","sex","celltype")) %>%
  dplyr::filter(n_genes >= MIN_GENES, celltype %in% CT_ORDER, age %in% c("middle","old")) %>%
  dplyr::mutate(pass=cos < 0 & beyond_null %in% TRUE, Sex=sexf(sex))

cts_gate <- read.csv(file.path(MASTER, "Step1_CellTypes.csv"), stringsAsFactors=FALSE) %>%
  dplyr::mutate(pass_middle=pass_middle %in% TRUE, pass_old=pass_old %in% TRUE)
calls <- readRDS(file.path(MASTER, "Step2_GeneCalls.RDS"))

log_norm <- function(cnt) { cnt <- round(as.matrix(cnt)); cnt <- cnt[rowSums(cnt >= 10) >= 2, , drop=FALSE]
  sf <- DESeq2::estimateSizeFactorsForMatrix(cnt); log2(sweep(cnt, 2, sf, "/") + 1) }
# B = aged ctrl - young; E = treated - aged ctrl; N = treated - young (group means, log scale)
stats3 <- function(X, y, o, e) {
  mY <- rowMeans(X[, y, drop=FALSE]); mO <- rowMeans(X[, o, drop=FALSE]); mE <- rowMeans(X[, e, drop=FALSE])
  vO <- apply(X[, o, drop=FALSE], 1, var) / length(o)
  B <- mO - mY; E <- mE - mO; N <- mE - mY
  c(cos_dc=(sum(B * E) + sum(vO)) / sqrt(sum(B^2) * sum(E^2)), BI_raw=sum(N^2) / sum(B^2))
}
sim_null <- function(X, y, o, e) {
  s  <- sqrt(rowMeans(cbind(apply(X[, y, drop=FALSE], 1, var), apply(X[, o, drop=FALSE], 1, var),
                            apply(X[, e, drop=FALSE], 1, var))))
  mu <- rowMeans(X)
  v  <- replicate(N_SIM, { Z <- matrix(rnorm(length(mu) * ncol(X), mu, s), nrow=length(mu), dimnames=dimnames(X))
                           stats3(Z, y, o, e)[["cos_dc"]] })
  c(null_mean=mean(v), null_sd=sd(v))
}

# DREADD per cell type: held-out cosine (own aging axis), gate, buffering index. Old females only.
dr <- NULL
if (any(c("b","e") %in% PANELS)) {
  ho_dir <- file.path(DR_BASE, "DREADD_Buffering_HeldOut")
  hsum <- read.csv(latest(ho_dir, "_DREADD_HeldOut_Summary\\.csv$"), stringsAsFactors=FALSE) %>%
    dplyr::filter(state %in% paste0("mScarlet_50_Approach_", DR_APPROACHES)) %>%
    dplyr::transmute(approach=as.character(approach), celltype, cos=cosine_heldout)
  f_own <- file.path(OUT, "SourceData_DREADD_CellType_Metrics.csv")
  f_src <- if (!OVERWRITE && file.exists(f_own)) f_own else
    tryCatch(latest(file.path(NE_BASE, "Paper_Figures"), "^SourceData_Fig5_DREADD_CellType_Metrics\\.csv$"), error=function(e) NA_character_)
  if (!is.na(f_src)) {
    cat("  DREADD metrics from", f_src, "\n")
    m <- read.csv(f_src, stringsAsFactors=FALSE) %>%
      dplyr::transmute(approach=as.character(approach), celltype, n_genes, cos_dc, BI_raw, null_mean, null_sd)
  } else {
    cat("  DREADD metrics: computing (bias-corrected cosine, simulated null, buffering index)\n")
    set.seed(20261022)                       # same seed and draw order as Buffering_Paper_Figures.R
    rows <- list()
    for (ap in DR_APPROACHES) {
      PB <- readRDS(file.path(ho_dir, paste0("pseudobulk_mScarlet_50_Approach_", ap, ".RDS")))
      for (ct in intersect(CT_ORDER, unique(PB$meta$celltype))) {
        d <- PB$meta[PB$meta$celltype == ct, ]
        y <- d$sample[d$age == "young" & d$condition == "PBS"]; o <- d$sample[d$age == "old" & d$condition == "PBS"]
        e <- d$sample[d$age == "old" & d$condition == "C21"]
        if (min(length(y), length(o), length(e)) < 2) next
        X <- log_norm(PB$counts[, c(y, o, e), drop=FALSE]); if (nrow(X) < MIN_GENES) next
        rows[[length(rows) + 1]] <- data.frame(approach=ap, celltype=ct, n_genes=nrow(X), t(stats3(X, y, o, e)), t(sim_null(X, y, o, e)))
      }
    }
    m <- dplyr::bind_rows(rows)
  }
  dr <- m %>% dplyr::inner_join(hsum, by=c("approach","celltype")) %>%
    dplyr::filter(n_genes >= MIN_GENES, celltype %in% CT_ORDER) %>%
    dplyr::mutate(pass=!is.na(cos) & cos < 0 & cos_dc < null_mean - 2 * null_sd)
  # DREADD buffering index as for NE, from the DREADD DESeq2 tables,
  #   BI = sum(log2FC OC21 vs YPBS ^2) / sum(log2FC OPBS vs YPBS ^2) over genes with a fold change in both
  if (BI_SOURCE == "original") {
    DR_DEG_DIRS <- c(B="DESeq2_mScarletPos_ApproachB_MergedIN", D="DESeq2_mScarletPos_ApproachD_MergedIN")
    rd <- function(ap, comp) { f <- list.files(file.path(DR_BASE, DR_DEG_DIRS[[ap]]), pattern=paste0("_", comp, "_ALL\\.csv$"),
                                                recursive=TRUE, full.names=TRUE)
      if (!length(f)) return(NULL); f <- sort(f, decreasing=TRUE)[1]; cat("  DREADD BI:", f, "\n")
      read.csv(f, stringsAsFactors=FALSE) %>% dplyr::select(celltype, gene, lfc=log2FoldChange) }
    bi_dr <- dplyr::bind_rows(lapply(DR_APPROACHES, function(ap) {
      n <- rd(ap, "OC21YPBS"); b <- rd(ap, "OPBSYPBS"); if (is.null(n) || is.null(b)) return(NULL)
      dplyr::inner_join(n, b, by=c("celltype","gene"), suffix=c("_N","_B")) %>% dplyr::filter(!is.na(lfc_N), !is.na(lfc_B)) %>%
        dplyr::group_by(celltype) %>% dplyr::summarise(BI_deseq=sum(lfc_N^2) / sum(lfc_B^2), .groups="drop") %>% dplyr::mutate(approach=ap) }))
    if (nrow(bi_dr)) dr <- dr %>% dplyr::left_join(bi_dr, by=c("approach","celltype")) %>%
      dplyr::mutate(BI_raw=dplyr::coalesce(BI_deseq, NA_real_)) %>% dplyr::select(-BI_deseq) else
      cat("  DREADD DESeq2 tables (OC21YPBS / OPBSYPBS) not found: DREADD BI recomputed from expression (differs from NE's)\n")
  }
  wsrc(dr, "DREADD_CellType_Metrics")
}

################################################################################
### 5b  COSINE ARC DIAGRAMS
################################################################################

arc_pts <- function(r, a0, a1, n=60) { a <- seq(a0, a1, length.out=n); data.frame(x=r * cos(a), y=r * sin(a)) }
# Arc / spoke scale: the half-circle spans cosine +ARC_LIM (right, along aging) to -ARC_LIM (left, opposes aging),
# linear in cosine (axis ticks are cosine values). Cell types beyond +/-ARC_LIM are drawn at the edge and marked.
ARC_TICKS <- seq(-ARC_LIM, ARC_LIM, length.out=5)       # ARC_LIM is set with the other settings at the top
cos2ang <- function(cs) pi / 2 * (1 - pmax(pmin(cs, ARC_LIM), -ARC_LIM) / ARC_LIM)   # cos = +L -> 0, 0 -> 90 deg, -L -> 180
beyond  <- function(cs) !is.na(cs) & abs(cs) > ARC_LIM
half <- dplyr::bind_rows(
  dplyr::bind_rows(data.frame(x=0, y=0), arc_pts(1.08, pi / 2, pi)) %>% dplyr::mutate(zone="Buffering (opposes aging)"),
  dplyr::bind_rows(data.frame(x=0, y=0), arc_pts(1.08, 0, pi / 2)) %>% dplyr::mutate(zone="Along aging"))
rep_half <- function(keys) dplyr::bind_rows(lapply(seq_len(nrow(keys)), function(k)
  dplyr::bind_cols(keys[rep(k, nrow(half)), , drop=FALSE], half)))
ZONE_COL <- c(`Buffering (opposes aging)`="#DCE9F5", `Along aging`="#F7E1E1")
ring_r <- function(ct) { i <- match(ct, rev(intersect(CT_ORDER, ct))); 0.35 + 0.65 * (i - 1) / pmax(length(ct) - 1, 1) }
arc_frame <- function(p, sub) p +
  annotate("segment", x=0, y=0, xend=1.12, yend=0, color=RED, linewidth=0.4, arrow=arrow(length=unit(1.1, "mm"), type="closed")) +
  annotate("text", x=1.12, y=-0.06, label="Aging", color=RED, size=1.9, hjust=1) +
  annotate("segment", x=0, y=0, xend=0, yend=1.08, color="grey55", linewidth=0.25, linetype="dashed") +
  annotate("segment", x=1.08 * cos(cos2ang(ARC_TICKS)), y=1.08 * sin(cos2ang(ARC_TICKS)),
           xend=1.13 * cos(cos2ang(ARC_TICKS)), yend=1.13 * sin(cos2ang(ARC_TICKS)), color="grey45", linewidth=0.3) +
  geom_text(data=data.frame(x=1.2 * cos(cos2ang(ARC_TICKS)), y=1.2 * sin(cos2ang(ARC_TICKS)), lab=as.character(ARC_TICKS)),
            aes(x, y, label=lab), inherit.aes=FALSE, color="grey35", size=1.6) +     # data frame: repeats correctly per facet
  annotate("text", x=0, y=1.29, label="cosine (aging vs NE)", color="grey40", size=1.6) +
  coord_equal(xlim=c(-1.25, 1.25), ylim=c(-0.42, 1.32), clip="off") + theme_void(base_size=7) +
  theme(legend.position="bottom", legend.box="vertical", strip.text=element_text(face="bold", size=7, margin=margin(b=3)),
        plot.title=element_text(face="bold", size=8), plot.subtitle=element_text(color="grey35", size=5.8)) +
  labs(subtitle=sub)

if ("b" %in% PANELS) {
  cat("\n== 5b arcs ==\n")
  # b1 NE: arc from the middle-age angle to the old-age angle, one ring per cell type
  ar <- cosd %>% dplyr::mutate(theta=cos2ang(cos)) %>% dplyr::select(Sex, celltype, age, theta, pass) %>%
    tidyr::pivot_wider(names_from=age, values_from=c(theta, pass)) %>% dplyr::filter(!is.na(theta_middle), !is.na(theta_old)) %>%
    dplyr::group_by(Sex) %>% dplyr::mutate(r=ring_r(celltype)) %>% dplyr::ungroup() %>%
    dplyr::mutate(status=factor(dplyr::case_when(pass_old ~ "passes at old age", theta_old < pi / 2 ~ "aging direction at old age",
                                                 TRUE ~ "no longer passes"),
                                levels=c("passes at old age","no longer passes","aging direction at old age")))
  wsrc(ar %>% dplyr::mutate(plot_deg_middle=theta_middle * 180 / pi, plot_deg_old=theta_old * 180 / pi), "Fig5b_NE_Angles")  # plotted position (scale +/-ARC_LIM)
  arcs <- dplyr::bind_rows(lapply(seq_len(nrow(ar)), function(k) arc_pts(ar$r[k], ar$theta_middle[k], ar$theta_old[k], 40) %>%
    dplyr::mutate(Sex=ar$Sex[k], celltype=ar$celltype[k], status=ar$status[k])))
  heads <- arcs %>% dplyr::group_by(Sex, celltype, status) %>% dplyr::slice_tail(n=2) %>%
    dplyr::summarise(x0=x[1], y0=y[1], x1=x[2], y1=y[2], .groups="drop")
  starts <- ar %>% dplyr::mutate(x=r * cos(theta_middle), y=r * sin(theta_middle))
  rings <- dplyr::bind_rows(lapply(seq_len(nrow(ar)), function(k) arc_pts(ar$r[k], 0, pi, 90) %>%
    dplyr::mutate(Sex=ar$Sex[k], celltype=ar$celltype[k])))
  p_ne <- ggplot() + geom_polygon(data=rep_half(ar %>% dplyr::distinct(Sex)), aes(x, y, fill=zone, group=zone)) +
    geom_path(data=rings, aes(x, y, group=celltype), color="white", linewidth=0.25) +
    geom_path(data=arcs, aes(x, y, group=celltype, color=status), linewidth=0.7, lineend="round") +
    geom_segment(data=heads, aes(x=x0, y=y0, xend=x1, yend=y1, color=status), linewidth=0.7,
                 arrow=arrow(length=unit(1.2, "mm"), type="closed")) +
    geom_point(data=starts, aes(x, y, fill=pass_middle), shape=21, size=1.1, stroke=0.3) +
    geom_text(data=ar, aes(x=-r - 0.02, y=-0.03, label=short_ct(celltype)), angle=90, hjust=1, size=1.6, color="grey25") +
    scale_fill_manual(values=c(ZONE_COL, `TRUE`="grey15", `FALSE`="white"),
                      breaks=c("Buffering (opposes aging)","Along aging","TRUE","FALSE"),
                      labels=c("Buffering (opposes aging)","Along aging","Middle age: passes","Middle age: does not pass"), name=NULL) +
    scale_color_manual(values=c(`passes at old age`=BLUE, `no longer passes`="grey55", `aging direction at old age`=RED),
                       name="Old age", drop=FALSE) +
    facet_wrap(~factor(paste("NE,", Sex), levels=paste("NE,", SEX_ORDER)), nrow=1)
  p_ne <- arc_frame(p_ne, paste0("One ring per cell type; dot = middle, arrowhead = old; left of 0 = buffering; cos -",
                                  ARC_LIM, " to +", ARC_LIM, " (beyond: at the edge)")) +
    labs(title="Angle between the aging and NE vectors, middle to old age")
  sv(p_ne, "Fig5b1_Arcs_NE", 125, 85)

  # b2 DREADD: point at the DREADD angle; open circle = NE old females on the same ring
  p_dr <- NULL
  if (!is.null(dr) && nrow(dr)) {
    ne_of <- cosd %>% dplyr::filter(sex == "female", age == "old") %>% dplyr::transmute(celltype, theta_ne=cos2ang(cos))
    da <- dr %>% dplyr::mutate(theta=cos2ang(cos), Panel=paste0("DREADD ", approach, ", old females")) %>%
      dplyr::group_by(Panel) %>% dplyr::mutate(r=ring_r(celltype)) %>% dplyr::ungroup() %>%
      dplyr::left_join(ne_of, by="celltype") %>%
      dplyr::mutate(x=r * cos(theta), y=r * sin(theta), xn=r * cos(theta_ne), yn=r * sin(theta_ne),
                    gate=factor(ifelse(pass, "passes", "does not pass"), levels=c("passes","does not pass")))
    wsrc(da %>% dplyr::mutate(plot_deg_DREADD=theta * 180 / pi, plot_deg_NE_old_female=theta_ne * 180 / pi), "Fig5b_DREADD_Angles")
    rings2 <- dplyr::bind_rows(lapply(seq_len(nrow(da)), function(k) arc_pts(da$r[k], 0, pi, 90) %>%
      dplyr::mutate(Panel=da$Panel[k], celltype=da$celltype[k])))
    p_dr <- ggplot() + geom_polygon(data=rep_half(da %>% dplyr::distinct(Panel)), aes(x, y, fill=zone, group=zone)) +
      geom_path(data=rings2, aes(x, y, group=celltype), color="white", linewidth=0.25) +
      geom_point(data=da %>% dplyr::filter(!is.na(xn)), aes(xn, yn), shape=21, fill="white", color="grey35", size=1.2, stroke=0.35) +
      geom_point(data=da, aes(x, y, color=approach, alpha=gate), shape=18, size=2.2) +
      geom_text(data=da, aes(x=-r - 0.02, y=-0.03, label=short_ct(celltype)), angle=90, hjust=1, size=1.6, color="grey25") +
      scale_fill_manual(values=ZONE_COL, name=NULL) + scale_color_manual(values=COL_AP, guide="none") +
      scale_alpha_manual(values=c(passes=1, `does not pass`=0.35), name="DREADD gate", drop=FALSE) +
      facet_wrap(~Panel, nrow=1)
    p_dr <- arc_frame(p_dr, paste0("Diamond = DREADD (own aging axis, held-out); open circle = NE, old females; scale cos -",
                                  ARC_LIM, " to +", ARC_LIM)) +
      labs(title="Angle between the aging and DREADD vectors")
    sv(p_dr, "Fig5b2_Arcs_DREADD", 125, 85)
  }
  if (!is.null(p_dr) && HAS_PATCH) sv(patchwork::wrap_plots(p_ne, p_dr, nrow=1), "Fig5b3_Arcs_NE_DREADD", 230, 85,
                                      note="combined; NE (middle -> old) | DREADD (old females)")

  # b4 spoke-style arc diagram (style of Cos_ArcDiagram.R): one spoke per cell type at the angle of its held-out
  #    cosine (same values as b1/b2); faded = cell type does not pass the buffering gate
  spoke_plot <- function(d, title, sub) {
    d <- d %>% dplyr::mutate(a=cos2ang(cos), x=cos(a), y=sin(a), clipped=beyond(cos),
                             gate=factor(ifelse(pass, "passes", "does not pass"), levels=c("passes","does not pass")))
    tk <- ARC_TICKS; ta <- cos2ang(tk)
    lab <- if (HAS_REPEL) ggrepel::geom_text_repel(aes(x=x * 1.12, y=y * 1.12, label=short_ct(celltype)), size=2.1, max.overlaps=30,
                                                   segment.size=0.25, segment.color="grey60", color="grey15") else
      geom_text(aes(x=x * 1.12, y=y * 1.12, label=short_ct(celltype)), size=2.1, color="grey15")
    ggplot(d, aes(color=celltype)) +
      annotate("path", x=cos(seq(0, pi, length.out=300)), y=sin(seq(0, pi, length.out=300)), color="grey80", linewidth=0.5) +
      annotate("segment", x=cos(ta), xend=cos(ta) * 1.07, y=sin(ta), yend=sin(ta) * 1.07, color="grey60", linewidth=0.3) +
      annotate("text", x=cos(ta[c(2, 4)]) * 0.9, y=sin(ta[c(2, 4)]) * 0.9, label=as.character(tk[c(2, 4)]), size=2.1, color="grey45") +
      annotate("segment", x=-1, xend=1, y=0, yend=0, color="grey60", linewidth=0.4) +
      annotate("segment", x=0, xend=0, y=0, yend=1, linetype="dashed", color="grey60", linewidth=0.3) +
      geom_segment(aes(x=0, xend=x, y=0, yend=y, alpha=gate), linewidth=0.8) +
      geom_point(aes(x, y, alpha=gate), size=2.6) +
      geom_point(data=d[d$clipped, ], aes(x, y), shape=21, fill=NA, color="black", size=3.8, stroke=0.4) + lab +
      annotate("text", x=-1.08, y=-0.10, label=paste0("Oppose\n(cos = -", ARC_LIM, ")"), size=2.2, color="#2166AC", fontface="bold") +
      annotate("text", x=1.08, y=-0.10, label=paste0("Align\n(cos = +", ARC_LIM, ")"), size=2.2, color="#D6604D", fontface="bold") +
      annotate("text", x=0, y=1.30, label="Orthogonal (cos = 0)", size=2.2, color="grey40") +
      scale_color_manual(values=CT_COLS, guide="none") +
      scale_alpha_manual(values=c(passes=0.95, `does not pass`=0.3), name="Buffering gate", drop=FALSE) +
      coord_fixed(xlim=c(-1.35, 1.35), ylim=c(-0.25, 1.38)) + theme_void(base_size=7) +
      theme(plot.title=element_text(face="bold", size=8, hjust=0.5), plot.subtitle=element_text(size=5.8, color="grey45", hjust=0.5),
            legend.position="bottom", plot.background=element_rect(fill="white", color=NA)) +
      labs(title=title, subtitle=sub)
  }
  SUB_NE <- paste0("Held-out cosine, aging vs NE vector\nleft of 0 = NE opposes aging; scale -", ARC_LIM, " to +", ARC_LIM,
                   ", ringed = beyond, drawn at the edge")
  bottom <- function(p) if (inherits(p, "patchwork")) p & theme(legend.position="bottom") else p
  wsrc(cosd %>% dplyr::mutate(plot_deg=cos2ang(cos) * 180 / pi, beyond_scale=beyond(cos)) %>%
         dplyr::select(sex, age, celltype, cos_heldout=cos, plot_deg, beyond_scale, pass), "Fig5b4_Spokes_NE")
  for (ag in c("middle","old")) {
    pl <- lapply(SEX_ORDER, function(s) { d <- cosd %>% dplyr::filter(Sex == s, age == ag)
      if (nrow(d)) spoke_plot(d, paste0("NE, ", s, ", ", tolower(AGE_LAB[[ag]])), SUB_NE) else NULL })
    if (any(!sapply(pl, is.null))) sv(bottom(side_by_side(pl)), paste0("Fig5b4_Spokes_NE_", ag), 200, 95, note="sexes side by side, own panels")
    for (s in names(SEX_LAB)) { d <- cosd %>% dplyr::filter(sex == s, age == ag)
      if (nrow(d)) sv(spoke_plot(d, paste0("NE, ", SEX_LAB[[s]], ", ", tolower(AGE_LAB[[ag]])), SUB_NE),
                      paste0("Fig5b4_Spokes_NE_", ag, "_", s), 100, 90) }
  }
  if (!is.null(dr) && nrow(dr)) {
    pl <- lapply(DR_APPROACHES, function(ap) { d <- dr %>% dplyr::filter(approach == ap, !is.na(cos))
      if (nrow(d)) spoke_plot(d, paste0("DREADD ", ap, ", old females"),
                              paste0("Held-out cosine, aging vs DREADD vector\nleft of 0 = DREADD opposes aging; scale -", ARC_LIM,
                                     " to +", ARC_LIM, ", ringed = beyond")) else NULL })
    if (any(!sapply(pl, is.null))) sv(bottom(side_by_side(pl)), "Fig5b4_Spokes_DREADD", 200, 95)
  }
}

################################################################################
### ED10b  HELD-OUT COSINE PER CELL TYPE (middle and old age; females | males)
###   Large symbol = mean of the 3 held-out fits; grey = each held-out HC mouse;
###   solid = passes the buffering gate (held-out cosine < 0 and bias-corrected cosine beyond the simulated null)
################################################################################

if ("h" %in% PANELS) {
  cat("\n== ED10b held-out cosine ==\n")
  th_ho <- theme_classic(base_size=7) +
    theme(strip.background=element_rect(fill="grey92", color="black", linewidth=0.3),
          strip.text=element_text(face="bold", size=7), axis.line=element_line(linewidth=0.3),
          axis.ticks=element_line(linewidth=0.3), plot.title=element_text(face="bold", hjust=0.5, size=8),
          plot.subtitle=element_text(color="grey35", hjust=0.5, size=6), legend.key.size=unit(3, "mm"))
  SEX_PANEL <- c(female="Females", male="Males")
  wsrc(cosd %>% dplyr::select(age, sex, celltype, n_genes, cos_heldout=cos, cos_dc, beyond_null, pass), "ED10b_HeldOutCosine")
  for (ag in c("middle","old")) {
    sm <- cosd %>% dplyr::filter(age == ag) %>%
      dplyr::mutate(panel=factor(SEX_PANEL[sex], levels=SEX_PANEL), shape_lab=SEX_PANEL[sex],
                    celltype=factor(celltype, levels=rev(intersect(CT_ORDER, unique(celltype)))),
                    gate=factor(ifelse(pass, "passes", "does not pass"), levels=c("passes","does not pass")))
    if (!nrow(sm)) next
    pts <- c2 %>% dplyr::filter(age == ag) %>% dplyr::semi_join(sm, by=c("sex","celltype")) %>%
      dplyr::mutate(panel=factor(SEX_PANEL[sex], levels=SEX_PANEL), celltype=factor(celltype, levels=levels(sm$celltype)))
    lim <- max(0.1, abs(c(sm$cos, pts$cosine)), na.rm=TRUE) * 1.05
    p <- ggplot(sm, aes(y=celltype)) + geom_vline(xintercept=0, linetype="dashed", color="grey40", linewidth=0.3) +
      geom_point(data=pts, aes(x=cosine, y=celltype), inherit.aes=FALSE, color="grey55", size=0.6) +
      geom_point(aes(x=cos, fill=celltype, shape=shape_lab, alpha=gate), size=2.3, stroke=0.3, color="black") +
      scale_fill_manual(values=CT_COLS, guide="none") +
      scale_shape_manual(values=c(Females=24, Males=21), guide="none") +
      scale_alpha_manual(values=c(passes=1, `does not pass`=0.3), name="Buffering gate", drop=FALSE) +
      scale_x_continuous(limits=c(-lim, lim)) +
      facet_grid(. ~ panel) + th_ho + theme(legend.position="bottom") +
      labs(title=paste("Direction of buffering:", AGE_LAB[[ag]]),
           subtitle="Large symbol = mean of held-out fits; grey = each held-out HC mouse",
           x="Held-out cosine similarity (aging, stimulation)", y=NULL)
    sv(p, paste0("ED10b_HeldOutCosine_", ag), 110, 25 + 5 * nlevels(sm$celltype))
  }
}

################################################################################
### GENE SETS (shared by 5c and 5d)
################################################################################

make_sets <- function(TIER) {
  AMP <- sub("^T", "A", TIER); N_AG <- if (TIER == "T1") 3 else 2
  w <- calls %>% dplyr::transmute(age, sex, celltype, gene, aging_lfc, ne_lfc, buf=.data[[TIER]] %in% TRUE,
                                  amp=.data[[AMP]] %in% TRUE, agechg=n_aging >= N_AG) %>%
    tidyr::pivot_wider(id_cols=c(sex, celltype, gene), names_from=age, values_from=c(aging_lfc, ne_lfc, buf, amp, agechg)) %>%
    dplyr::inner_join(cts_gate %>% dplyr::select(sex, celltype, pass_middle, pass_old), by=c("sex","celltype"))
  for (v in c("buf_middle","buf_old","amp_middle","amp_old","agechg_middle","agechg_old"))
    if (!v %in% names(w)) w[[v]] <- NA
  both <- !is.na(w$buf_middle) & !is.na(w$buf_old)
  T_ <- function(x) x %in% TRUE
  long <- dplyr::bind_rows(
    w %>% dplyr::mutate(set="Buffered, middle",  member=T_(buf_middle) & pass_middle, bg=T_(agechg_middle) & pass_middle),
    w %>% dplyr::mutate(set="Buffered, old",     member=T_(buf_old) & pass_old,       bg=T_(agechg_old) & pass_old),
    w %>% dplyr::mutate(set="Lost",              member=both & T_(buf_middle) & !T_(buf_old), bg=both & T_(buf_middle)),
    w %>% dplyr::mutate(set="Maintained",        member=both & T_(buf_middle) & T_(buf_old),  bg=both & T_(buf_middle)),
    w %>% dplyr::mutate(set="Amplified, middle", member=T_(amp_middle), bg=T_(agechg_middle)),
    w %>% dplyr::mutate(set="Amplified, old",    member=T_(amp_old),    bg=T_(agechg_old))) %>%
    dplyr::filter(bg) %>% dplyr::select(sex, celltype, gene, set, member)
  # one block per gene for the heatmaps (first match wins)
  blk <- w %>% dplyr::mutate(block=dplyr::case_when(
      T_(buf_middle) & T_(buf_old) ~ BLOCKS[1],
      T_(buf_middle) & pass_middle ~ BLOCKS[2],
      T_(buf_old) & pass_old ~ BLOCKS[3],
      T_(amp_middle) ~ BLOCKS[4],
      T_(amp_old) ~ BLOCKS[5], TRUE ~ NA_character_),
      ref_lfc=ifelse(block %in% BLOCKS[c(3, 5)], aging_lfc_old, dplyr::coalesce(aging_lfc_middle, aging_lfc_old)),
      score=dplyr::case_when(block == BLOCKS[1] ~ pmin(-aging_lfc_middle * ne_lfc_middle, -aging_lfc_old * ne_lfc_old),
                             block == BLOCKS[2] ~ -aging_lfc_middle * ne_lfc_middle, block == BLOCKS[3] ~ -aging_lfc_old * ne_lfc_old,
                             block == BLOCKS[4] ~ aging_lfc_middle * ne_lfc_middle, block == BLOCKS[5] ~ aging_lfc_old * ne_lfc_old)) %>%
    dplyr::filter(!is.na(block), !is.na(ref_lfc)) %>%
    dplyr::transmute(sex, celltype, gene, block, ref_lfc, score, dir=ifelse(ref_lfc > 0, "Up with age", "Down with age"))
  list(long=long, blk=blk)
}
if (any(c("c","d") %in% PANELS)) {
  S_SET <- make_sets(SET_TIER)
  S_LAB <- if (LAB_TIER == SET_TIER) S_SET else make_sets(LAB_TIER)
  sizes <- S_SET$long %>% dplyr::group_by(sex, celltype, set) %>%
    dplyr::summarise(n_genes=sum(member), n_universe=dplyr::n(), .groups="drop")
  wsrc(sizes, paste0("Fig5_GeneSet_Sizes_", SET_TIER))
  wsrc(S_SET$long %>% dplyr::filter(member) %>% dplyr::select(-member), paste0("Fig5_GeneSet_Members_", SET_TIER))
  cat("\nGene sets (", SET_TIER, "), genes per sex:\n", sep="")
  print(as.data.frame(sizes %>% dplyr::group_by(sex, set) %>% dplyr::summarise(n_genes=sum(n_genes), .groups="drop") %>%
    tidyr::pivot_wider(names_from=set, values_from=n_genes)), row.names=FALSE)
}

################################################################################
### 5c  HEATMAPS
################################################################################

group_of <- function(m) dplyr::case_when(
  m$age_norm == "young" & m$condition == "HC" ~ "Young HC",
  m$age_norm == MID & m$condition == "HC" ~ "Middle HC", m$age_norm == MID & m$condition == "EE" ~ "Middle NE",
  m$age_norm == "old" & m$condition == "HC" ~ "Old HC", m$age_norm == "old" & m$condition == "EE" ~ "Old NE")
expr_z <- function(sx, ct, genes) {
  m <- pb_meta[pb_meta$celltype == ct & pb_meta$sex_norm == sx, ]; g <- group_of(m); keep <- !is.na(g)
  X <- log_norm(pb_all[, m$sample[keep], drop=FALSE]); genes <- intersect(genes, rownames(X))
  if (length(genes) < 2) return(NULL)
  Z <- t(scale(t(X[genes, , drop=FALSE]))); Z[is.na(Z)] <- 0
  as.data.frame(as.table(Z)) %>% setNames(c("gene","sample","z")) %>%
    dplyr::mutate(gene=as.character(gene), sample=as.character(sample),
                  group=factor(setNames(g[keep], m$sample[keep])[sample], levels=GROUPS))
}
heat_plot <- function(h, title, labeled) {
  h <- h %>% dplyr::mutate(rowlab=factor(paste0(BLOCK_SHORT[as.character(block)], "\n", ifelse(dir == "Up with age", "up", "down"), " with age")),
                           rowlab=factor(rowlab, levels=unique(rowlab[order(block, dir)])))
  ggplot(h, aes(sample, gene, fill=pmax(pmin(z, 2), -2))) + geom_tile() +
    facet_grid(rowlab ~ group, scales="free", space="free",
               labeller=labeller(group=function(x) sub(" ", "\n", x))) +
    scale_fill_gradient2(low=BLUE, mid="white", high=RED, limits=c(-2, 2), name="z-score") + th +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line=element_blank(),
          axis.text.y=if (labeled) element_text(size=3.8, face="italic") else element_blank(),
          axis.ticks.y=if (labeled) element_line(linewidth=0.2) else element_blank(),
          strip.text.x=element_text(size=5), strip.text.y=element_text(size=4.6, angle=0, hjust=0),
          panel.spacing=unit(0.35, "mm")) +
    labs(title=title, x=NULL, y=NULL)
}

if ("c" %in% PANELS) {
  cat("\n== 5c heatmaps ==\n")
  hdir <- file.path(OUT, "Fig5c_Heatmaps"); dir.create(hdir, showWarnings=FALSE)
  units <- cts_gate %>% dplyr::filter(celltype %in% CT_ORDER) %>% dplyr::arrange(factor(celltype, levels=CT_ORDER))
  summ <- list(); src_lab <- list(); src_unl <- list()
  for (ct in unique(units$celltype)) {
    pl_lab <- list(); pl_unl <- list(); n_g <- 0
    for (sx in c("male","female")) {
      if (!any(units$celltype == ct & units$sex == sx)) next
      b_set <- S_SET$blk %>% dplyr::filter(sex == sx, celltype == ct)
      if (nrow(b_set) < 2) next
      # c2 unlabeled: all SET_TIER genes (capped per block)
      u <- b_set %>% dplyr::group_by(block) %>% dplyr::slice_max(abs(ref_lfc), n=UNLAB_MAX, with_ties=FALSE) %>% dplyr::ungroup()
      hz <- expr_z(sx, ct, u$gene); if (is.null(hz)) next
      ord <- u %>% dplyr::arrange(factor(block, levels=BLOCKS), dir, ref_lfc)
      hu <- hz %>% dplyr::inner_join(u %>% dplyr::select(gene, block, dir, ref_lfc), by="gene") %>%
        dplyr::mutate(gene=factor(gene, levels=ord$gene), block=factor(block, levels=BLOCKS),
                      dir=factor(dir, levels=c("Up with age","Down with age")), Sex=sexf(sx))
      src_unl[[length(src_unl) + 1]] <- hu %>% dplyr::mutate(celltype=ct)
      cnt <- u %>% dplyr::count(block)
      pl_unl[[sx]] <- heat_plot(hu, sprintf("%s, %s", ct, SEX_LAB[[sx]]), FALSE) +
        labs(subtitle=paste0(SET_TIER, " genes: ", paste(sprintf("%s %d", BLOCK_SHORT[cnt$block], cnt$n), collapse="; ")))
      # c3 summary: mean aging-oriented z per block x group (all SET_TIER genes, uncapped)
      za <- expr_z(sx, ct, b_set$gene)
      if (!is.null(za)) summ[[length(summ) + 1]] <- za %>% dplyr::inner_join(b_set %>% dplyr::select(gene, block, ref_lfc), by="gene") %>%
        dplyr::mutate(zo=z * sign(ref_lfc)) %>% dplyr::group_by(block, group) %>%
        dplyr::summarise(z_aging=mean(zo), n_genes=dplyr::n_distinct(gene), .groups="drop") %>%
        dplyr::mutate(sex=sx, celltype=ct)
      # c1 labeled: LAB_TIER genes only (no fallback)
      b_lab <- S_LAB$blk %>% dplyr::filter(sex == sx, celltype == ct)
      nlab <- b_lab %>% dplyr::count(block, name="n_lab")
      use_lab <- nlab$block[nlab$n_lab >= LAB_MIN]
      L <- b_lab %>% dplyr::filter(block %in% use_lab) %>% dplyr::mutate(tier=LAB_TIER) %>%   # no fallback: blocks with < LAB_MIN T1 genes are not shown
        dplyr::group_by(block) %>% dplyr::slice_max(score, n=LAB_N, with_ties=FALSE) %>% dplyr::ungroup()
      hz <- if (nrow(L) >= 2) expr_z(sx, ct, L$gene) else NULL
      if (!is.null(hz)) {
        ordL <- L %>% dplyr::arrange(factor(block, levels=BLOCKS), dir, ref_lfc)
        hl <- hz %>% dplyr::inner_join(L %>% dplyr::select(gene, block, dir, tier), by="gene") %>%
          dplyr::mutate(gene=factor(gene, levels=ordL$gene), block=factor(block, levels=BLOCKS),
                        dir=factor(dir, levels=c("Up with age","Down with age")), Sex=sexf(sx))
        src_lab[[length(src_lab) + 1]] <- hl %>% dplyr::mutate(celltype=ct); n_g <- max(n_g, dplyr::n_distinct(hl$gene))
        tiers <- L %>% dplyr::distinct(block, tier)
        pl_lab[[sx]] <- heat_plot(hl, sprintf("%s, %s", ct, SEX_LAB[[sx]]), TRUE) +
          labs(subtitle=sprintf("Top %d %s genes per block by buffering strength; blocks with < %d %s genes not shown", LAB_N, LAB_TIER, LAB_MIN, LAB_TIER))
      }
    }
    tag <- tagify(ct)
    if (length(pl_lab)) sv(side_by_side(pl_lab), paste0("Fig5c1_Labeled_", tag), 90 * length(pl_lab), 35 + 2.0 * n_g, dir=hdir)
    if (length(pl_unl)) sv(side_by_side(pl_unl), paste0("Fig5c2_Unlabeled_", tag), 80 * length(pl_unl), 110, dir=hdir)
  }
  if (length(src_lab)) wsrc(dplyr::bind_rows(src_lab), "Fig5c1_Labeled")
  if (length(src_unl)) wsrc(dplyr::bind_rows(src_unl) %>% dplyr::select(celltype, Sex, block, dir, gene, sample, group, z), "Fig5c2_Unlabeled")
  if (length(summ)) {
    sm <- dplyr::bind_rows(summ) %>% dplyr::mutate(Sex=sexf(sex), celltype=factor(celltype, levels=rev(CT_ORDER)),
                                                   block=factor(block, levels=BLOCKS))
    wsrc(sm, "Fig5c3_Summary")
    lim <- max(abs(sm$z_aging), na.rm=TRUE)
    pl <- lapply(SEX_ORDER, function(s) { d <- sm %>% dplyr::filter(Sex == s); if (!nrow(d)) return(NULL)
      ggplot(d, aes(group, celltype, fill=z_aging)) + geom_tile(color="white", linewidth=0.3) +
        facet_grid(. ~ block, labeller=labeller(block=label_wrap_gen(14))) +
        scale_fill_gradient2(low=BLUE, mid="white", high=RED, limits=c(-lim, lim), name="Aging-oriented z\n(+ = old-like)") +
        scale_y_discrete(labels=short_ct) + th +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=5), axis.line=element_blank(),
              strip.text=element_text(size=5), panel.spacing=unit(0.8, "mm")) +
        labs(title=s, x=NULL, y=NULL) })
    sv(side_by_side(pl) + if (HAS_PATCH) patchwork::plot_annotation(
         title="Gene sets per cell type: mean expression, oriented to aging",
         subtitle=paste0(SET_TIER, " genes; z per gene across mice x sign of its aging change (gene counts in source data)"),
         theme=ann_th) else NULL,
       "Fig5c3_Summary_AllCellTypes", 230, 50 + 6 * nlevels(droplevels(sm$celltype)), note="compact heatmap without genes")
  }

  # c4 fold-change strips, one file per cell type x gene set (labeled and unlabeled), sexes side by side.
  # Columns = aging and NE log2 FC (mean of the 3 held-out fits). Single-age sets show only their own age;
  # lost and maintained show middle and old side by side. Grey = not tested at that age.
  sdir <- file.path(OUT, "Fig5c4_FoldChangeStrips"); dir.create(sdir, showWarnings=FALSE)
  SET_COLS <- list(`Buffered, middle`=c("Aging, middle","NE, middle"), `Buffered, old`=c("Aging, old","NE, old"),
                   Lost=c("Aging, middle","NE, middle","Aging, old","NE, old"), Maintained=c("Aging, middle","NE, middle","Aging, old","NE, old"),
                   `Amplified, middle`=c("Aging, middle","NE, middle"), `Amplified, old`=c("Aging, old","NE, old"))
  MID_SETS <- c("Buffered, middle","Lost","Maintained")
  lfc_long <- calls %>% dplyr::select(sex, celltype, gene, age, aging_lfc, ne_lfc) %>%
    tidyr::pivot_longer(c(aging_lfc, ne_lfc), names_to="what", values_to="lfc") %>%
    dplyr::mutate(col=paste0(ifelse(what == "aging_lfc", "Aging, ", "NE, "), age)) %>% dplyr::select(-what, -age)
  lw <- calls %>% dplyr::select(sex, celltype, gene, age, aging_lfc, ne_lfc) %>%
    tidyr::pivot_wider(names_from=age, values_from=c(aging_lfc, ne_lfc))
  coseq <- NULL; f_fold <- file.path(MASTER, FOLD_FILE)
  if (MID_COSEQ && file.exists(f_fold)) coseq <- read.csv(f_fold, stringsAsFactors=FALSE) %>%
    dplyr::filter(age == "middle", coseq_opp_nesig %in% TRUE) %>% dplyr::distinct(sex, celltype, gene)
  if (MID_COSEQ && is.null(coseq)) cat("  MID_COSEQ: fold-support file not found, named middle-age genes not filtered:", f_fold, "\n")
  strip_genes <- function(S, sx, ct, st) {
    g <- S$long %>% dplyr::filter(sex == sx, celltype == ct, set == st, member) %>% dplyr::pull(gene)
    d <- lw %>% dplyr::filter(sex == sx, celltype == ct, gene %in% g)
    am <- d$aging_lfc_middle; ao <- d$aging_lfc_old; nm_ <- d$ne_lfc_middle; no <- d$ne_lfc_old
    d %>% dplyr::mutate(ref=if (SET_AGE[[st]] == "old") ao else am,
                        score=switch(st, `Buffered, middle`=-am * nm_, Lost=-am * nm_, `Buffered, old`=-ao * no,
                                     Maintained=pmin(-am * nm_, -ao * no), `Amplified, middle`=am * nm_, `Amplified, old`=ao * no),
                        dir=ifelse(ref > 0, "Up with age", "Down with age")) %>%
      dplyr::filter(!is.na(ref)) %>% dplyr::select(gene, ref, score, dir)
  }
  strip_plot <- function(G, sx, ct, st, labeled, sub) {
    cols <- SET_COLS[[st]]
    ord <- G %>% dplyr::arrange(dir, ref)
    d <- lfc_long %>% dplyr::filter(sex == sx, celltype == ct, gene %in% G$gene, col %in% cols) %>% dplyr::select(gene, col, lfc) %>%
      tidyr::complete(gene=G$gene, col=cols) %>% dplyr::inner_join(G %>% dplyr::select(gene, dir), by="gene") %>%
      dplyr::mutate(col=factor(col, levels=cols), gene=factor(gene, levels=ord$gene),
                    dir=factor(dir, levels=c("Up with age","Down with age")))
    p <- ggplot(d, aes(col, gene, fill=pmax(pmin(lfc, STRIP_LIM), -STRIP_LIM))) +
      geom_tile(color=if (labeled) "white" else NA, linewidth=0.15)
    if (length(cols) == 4) p <- p + geom_vline(xintercept=2.5, color="grey20", linewidth=0.35)
    p <- p + facet_grid(dir ~ ., scales="free_y", space="free_y") +
      scale_fill_gradient2(low=BLUE, mid="white", high=RED, limits=c(-STRIP_LIM, STRIP_LIM), na.value="grey85", name="log2 FC") +
      scale_x_discrete(expand=c(0, 0)) + th +
      theme(axis.text.x=element_text(angle=45, hjust=1), axis.line=element_blank(), strip.text.y=element_text(angle=0, size=5),
            axis.text.y=if (labeled) element_text(size=3.8, face="italic") else element_blank(),
            axis.ticks.y=if (labeled) element_line(linewidth=0.2) else element_blank()) +
      labs(title=SEX_LAB[[sx]], subtitle=sub, x=NULL, y=NULL)
    list(p=p, d=d %>% dplyr::mutate(sex=sx, celltype=ct), n=nrow(G))
  }
  src_strip <- list()
  for (ct in unique(units$celltype)) for (st in SETS) {
    nm <- paste0(tagify(ct), "_", tagify(st))
    f_lab <- paste0("Fig5c4_", nm, "_labeled"); f_unl <- paste0("Fig5c4_", nm, "_unlabeled")
    if (done(f_lab, sdir) && done(f_unl, sdir)) next
    lab <- list(); unl <- list(); n_lab <- 0
    for (sx in c("male","female")) {
      Gs <- strip_genes(S_SET, sx, ct, st); if (nrow(Gs) < 2) next
      Gl <- strip_genes(S_LAB, sx, ct, st); tier <- LAB_TIER
      if (nrow(Gl) < LAB_MIN) Gl <- Gl[0, ]          # no fallback: no labeled strip when < LAB_MIN T1 genes
      if (!is.null(coseq) && st %in% MID_SETS) {
        Gc <- Gl %>% dplyr::semi_join(coseq %>% dplyr::filter(sex == sx, celltype == ct), by="gene")
        if (nrow(Gc) >= LAB_MIN) { Gl <- Gc; tier <- paste(tier, "+ NE-sig vs co-sequenced HC") }
      }
      Gl <- Gl %>% dplyr::filter(score > 0) %>% dplyr::group_by(dir) %>% dplyr::slice_max(score, n=STRIP_N, with_ties=FALSE) %>% dplyr::ungroup()
      Gu <- Gs %>% dplyr::group_by(dir) %>% dplyr::slice_max(abs(ref), n=UNLAB_MAX, with_ties=FALSE) %>% dplyr::ungroup()
      if (nrow(Gl) >= 2) { a <- strip_plot(Gl, sx, ct, st, TRUE, sprintf("%s genes; top %d per direction by %s strength (n = %d)",
                                            tier, STRIP_N, if (grepl("^Amplified", st)) "amplification" else "buffering", nrow(Gl)))
                           lab[[sx]] <- a$p; n_lab <- max(n_lab, a$n) }
      b <- strip_plot(Gu, sx, ct, st, FALSE, sprintf("%s genes (n = %d)", SET_TIER, nrow(Gu)))
      unl[[sx]] <- b$p
      src_strip[[length(src_strip) + 1]] <- b$d %>% dplyr::mutate(set=st)
    }
    if (!length(unl)) next
    ttl <- if (HAS_PATCH) patchwork::plot_annotation(title=sprintf("%s: %s", ct, st),
             subtitle="Mean log2 FC, 3 held-out fits; grey = not tested",
             theme=ann_th) else NULL
    ncol_ <- length(SET_COLS[[st]])
    if (length(lab)) sv(side_by_side(lab) + ttl, f_lab, (12 * ncol_ + 40) * length(lab) + 20, 40 + 2.2 * n_lab, dir=sdir)
    sv(side_by_side(unl) + ttl, f_unl, (10 * ncol_ + 30) * length(unl) + 30, 110, dir=sdir)
  }
  if (length(src_strip)) wsrc(dplyr::bind_rows(src_strip), "Fig5c4_FoldChangeStrips_Unlabeled")
}

################################################################################
### 5d  PATHWAY ENRICHMENT (GO)
################################################################################

if ("d" %in% PANELS) {
  cat("\n== 5d GO enrichment ==\n")
  ok <- requireNamespace("clusterProfiler", quietly=TRUE) && requireNamespace("org.Mm.eg.db", quietly=TRUE)
  if (!ok) cat("  skipped: install clusterProfiler and org.Mm.eg.db (BiocManager::install(c('clusterProfiler','org.Mm.eg.db')))\n")
  if (ok) {
    suppressPackageStartupMessages({ library(clusterProfiler); library(org.Mm.eg.db) })
    L <- S_SET$long
    ct2class <- setNames(rep(names(CLASSES), lengths(CLASSES)), unlist(CLASSES))
    # jobs: per cell type, per class (union over its cell types), superpooled (union over all)
    pooled <- function(d, unit) d %>% dplyr::group_by(sex, set, gene) %>% dplyr::summarise(member=any(member), .groups="drop") %>%
      dplyr::mutate(unit=unit)
    J <- dplyr::bind_rows(
      L %>% dplyr::mutate(level="Cell type", unit=celltype) %>% dplyr::select(level, unit, sex, set, gene, member),
      dplyr::bind_rows(lapply(names(CLASSES), function(cl) pooled(L %>% dplyr::filter(celltype %in% CLASSES[[cl]]), cl))) %>%
        dplyr::mutate(level="Cell class"),
      pooled(L, "All cell types") %>% dplyr::mutate(level="Superpooled"))
    keys <- J %>% dplyr::distinct(level, unit, sex, set)
    f_go <- file.path(OUT, paste0("Fig5d_GO_", GO_ONT, "_", SET_TIER, "_AllResults.csv"))
    if (file.exists(f_go)) { go <- read.csv(f_go, stringsAsFactors=FALSE); cat("  loaded", f_go, "\n") } else {
      run_go <- function(i) {
        k <- keys[i, ]; d <- J %>% dplyr::filter(level == k$level, unit == k$unit, sex == k$sex, set == k$set)
        genes <- unique(d$gene[d$member]); if (length(genes) < GO_MIN_SET) return(NULL)
        uni <- unique(J$gene[J$level == "Cell type" & J$sex == k$sex & J$set == k$set])   # pooled across cell types, as in the lab scripts
        eg <- tryCatch(clusterProfiler::enrichGO(genes, OrgDb=org.Mm.eg.db, keyType="SYMBOL", ont=GO_ONT, universe=uni,
                                                 pvalueCutoff=GO_FDR, qvalueCutoff=0.2, minGSSize=10, maxGSSize=500, readable=TRUE), error=function(e) NULL)
        if (is.null(eg) || !nrow(as.data.frame(eg))) return(NULL)
        if (GO_SIMPLIFY && requireNamespace("GOSemSim", quietly=TRUE))
          eg <- tryCatch(clusterProfiler::simplify(eg, cutoff=0.7, by="p.adjust", select_fun=min), error=function(e) eg)
        as.data.frame(eg) %>% dplyr::mutate(level=k$level, unit=k$unit, sex=k$sex, set=k$set, n_set=length(genes),
                                            n_universe=length(uni), .before=1)
      }
      cat(sprintf("  %d gene sets on %d cores...\n", nrow(keys), N_CORES))
      rg <- mclapply(seq_len(nrow(keys)), function(i) tryCatch(run_go(i), error=function(e) conditionMessage(e)), mc.cores=N_CORES)
      er <- unlist(Filter(is.character, rg)); if (length(er)) { cat("  errors:\n"); print(head(er)) }
      go <- dplyr::bind_rows(Filter(is.data.frame, rg))
      if (nrow(go)) write.csv(go, f_go, row.names=FALSE)
    }
    if (!nrow(go)) cat("  no enrichment results\n") else {
      frac <- function(x) sapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
      go <- go %>% dplyr::mutate(gene_ratio=frac(GeneRatio), fold=gene_ratio / frac(BgRatio), nlp=-log10(p.adjust),
                                 term=ifelse(nchar(Description) > 55, paste0(substr(Description, 1, 52), "..."), Description),
                                 set=factor(set, levels=SETS))
      cat("\n  significant terms (FDR <", GO_FDR, ") per level x sex x set:\n")
      print(as.data.frame(go %>% dplyr::filter(p.adjust < GO_FDR) %>% dplyr::count(level, sex, set) %>%
        tidyr::pivot_wider(names_from=set, values_from=n, values_fill=0) %>% dplyr::arrange(level, sex)), row.names=FALSE)
      wsrc(go %>% dplyr::filter(p.adjust < GO_FDR) %>% dplyr::select(level, unit, sex, set, ID, Description, GeneRatio, BgRatio,
                                                                     fold, pvalue, p.adjust, geneID), paste0("Fig5d_GO_Significant_", SET_TIER))
      dotplot_go <- function(d, title, facet_unit=FALSE) {
        top <- d %>% dplyr::filter(p.adjust < GO_FDR) %>% dplyr::group_by(set) %>%
          dplyr::slice_min(p.adjust, n=GO_TOP, with_ties=FALSE) %>% dplyr::ungroup()
        if (!nrow(top)) return(NULL)
        dd <- d %>% dplyr::filter(ID %in% top$ID) %>%
          dplyr::mutate(term=factor(term, levels=rev(unique(top$term[order(top$set, top$p.adjust)]))),
                        sig=p.adjust < GO_FDR)
        g <- ggplot(dd, aes(set, term)) +
          geom_point(aes(size=pmin(fold, 6), fill=pmin(nlp, 10), alpha=sig), shape=21, stroke=0.2) +
          scale_fill_gradient(low="#FEE0D2", high="#A50F15", limits=c(0, 10), name="-log10 FDR") +
          scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=0.25), labels=c(`TRUE`=paste("FDR <", GO_FDR), `FALSE`="n.s."), name=NULL) +
          scale_size_continuous(range=c(0.6, 3), limits=c(0, 6), name="Fold enrichment") + scale_x_discrete(drop=FALSE) + th +
          theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(size=5), legend.position="right") +
          labs(title=title, x=NULL, y=NULL)
        if (facet_unit) g <- g + facet_grid(. ~ unit)
        g
      }
      # d1 superpooled, sexes side by side
      pl <- lapply(SEX_ORDER, function(s) dotplot_go(go %>% dplyr::filter(level == "Superpooled", SEX_LAB[sex] == s), s))
      if (any(!sapply(pl, is.null))) sv(side_by_side(pl), "Fig5d1_GO_Superpooled", 230, 35 + 3 * GO_TOP * length(SETS) / 1.6,
                                        note="all cell types pooled; universe = genes that could be in the set")
      # d3 per cell class
      for (s in names(SEX_LAB)) {
        p <- dotplot_go(go %>% dplyr::filter(level == "Cell class", sex == s) %>% dplyr::mutate(unit=factor(unit, levels=names(CLASSES))),
                        paste(SEX_LAB[[s]], "- per cell class"), facet_unit=TRUE)
        if (!is.null(p)) sv(p, paste0("Fig5d3_GO_CellClass_", s), 230, 35 + 3 * GO_TOP * length(SETS) / 1.4)
      }
      # d2 per cell type heatmap, one per set (top terms from any cell type)
      for (st in SETS) {
        pl <- lapply(SEX_ORDER, function(s) {
          d <- go %>% dplyr::filter(level == "Cell type", SEX_LAB[sex] == s, set == st)
          top <- d %>% dplyr::filter(p.adjust < GO_FDR) %>% dplyr::group_by(unit) %>% dplyr::slice_min(p.adjust, n=GO_TOP, with_ties=FALSE) %>%
            dplyr::ungroup() %>% dplyr::count(ID, term) %>% dplyr::arrange(dplyr::desc(n))
          if (!nrow(top)) return(NULL)
          top <- head(top, 25)
          dd <- d %>% dplyr::filter(ID %in% top$ID) %>%
            dplyr::mutate(term=factor(term, levels=rev(top$term)), unit=factor(unit, levels=intersect(CT_ORDER, unit)))
          ggplot(dd, aes(unit, term, fill=pmin(nlp, 10))) + geom_tile(color="white", linewidth=0.3) +
            geom_text(aes(label=ifelse(p.adjust < GO_FDR, "*", "")), size=2, vjust=0.75) +
            scale_fill_gradient(low="#FFF5F0", high="#A50F15", limits=c(0, 10), name="-log10 FDR") + scale_x_discrete(labels=short_ct) + th +
            theme(axis.text.x=element_text(angle=45, hjust=1), axis.text.y=element_text(size=5), axis.line=element_blank()) +
            labs(title=s, x=NULL, y=NULL) })
        if (any(!sapply(pl, is.null)))
          sv(side_by_side(pl) + if (HAS_PATCH) patchwork::plot_annotation(title=paste("GO", GO_ONT, "per cell type:", st),
                                                                          subtitle=paste("* FDR <", GO_FDR), theme=ann_th) else NULL,
             paste0("Fig5d2_GO_PerCellType_", tagify(st)), 230, 110)
      }
    }
  }
}

################################################################################
### 5e  QUADRANT: held-out cosine vs buffering index (raw)
################################################################################

# Same layout as the original figure: x = cosine on a fixed -1..1 scale, y = buffering index on a linear scale from 0;
# shaded zones: lower left = full buffering (opposes aging and closer to young), top band = amplification (further from young).
COL_ZONE_BUF <- "#F2F6FB"; COL_ZONE_AMP <- "#FDF3EC"; COL_TXT_BUF <- "#2F5597"; COL_TXT_AMP <- "#C55A2B"
quad_base <- function(d, ymax) ggplot(d, aes(cos, BI_raw)) +
  annotate("rect", xmin=-Inf, xmax=0, ymin=-Inf, ymax=1, fill=COL_ZONE_BUF) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=1, ymax=Inf, fill=COL_ZONE_AMP) +
  geom_hline(yintercept=1, color="grey30", linewidth=0.4) + geom_vline(xintercept=0, color="grey30", linewidth=0.4) +
  annotate("text", x=-1.05, y=0.03 * ymax, label="FULL BUFFERING\n(oppose + compress)", color=COL_TXT_BUF, size=1.9, lineheight=0.9, hjust=0, vjust=0) +
  annotate("text", x=0.98, y=0.985 * ymax, label="AMPLIFICATION\n(align + expand)", color=COL_TXT_AMP, size=2, lineheight=0.9, hjust=1, vjust=1) +
  scale_x_continuous(limits=c(-1.1, 1.1), breaks=c(-1, -0.5, 0, 0.5, 1), expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, ymax), expand=c(0, 0)) +
  th + theme(panel.border=element_rect(fill=NA, color="black", linewidth=0.4), legend.position="bottom") +
  guides(alpha=guide_legend(override.aes=list(fill="grey40", color="grey40", size=2)))
lab_layer <- function(d) if (HAS_REPEL)
  ggrepel::geom_text_repel(data=d, aes(label=short_ct(as.character(celltype)), color=celltype), size=1.8, segment.size=0.2,
                           max.overlaps=30, min.segment.length=0.2, show.legend=FALSE) else
  geom_text(data=d, aes(label=short_ct(as.character(celltype)), color=celltype), size=1.8, vjust=-0.9, show.legend=FALSE)
ymax_of <- function(x) max(1.25, max(x, na.rm=TRUE) * 1.08)

if ("e" %in% PANELS) {
  cat("\n== 5e quadrant ==\n")
  # Buffering index from step bi (*_BI_Summary.csv: bi = sum(OEEYHC LFC^2) / sum(OHCYHC LFC^2), from the DESeq2 comparisons)
  bi_orig <- NULL
  if (BI_SOURCE == "original") {
    dirs <- sort(list.dirs(file.path(NE_BASE, "Figures"), recursive=FALSE), decreasing=TRUE)
    dirs <- dirs[grepl("Buffering_Index", basename(dirs))]
    dirs <- dirs[sapply(dirs, function(d) length(list.files(d, pattern="_BI_Summary\\.csv$")) > 0)]
    dirs <- dirs[order(file.mtime(dirs), decreasing=TRUE)]   # newest folder that holds BI files
    fs <- if (length(dirs)) list.files(dirs[1], pattern="_BI_Summary\\.csv$", full.names=TRUE) else character(0)
    if (!length(fs)) stop("No *_BI_Summary.csv under ", file.path(NE_BASE, "Figures"), ": run pipeline step bi first") else {
      cat("  BI from:", paste(basename(fs), collapse=", "), "in", dirs[1], "\n")
      raw <- dplyr::bind_rows(lapply(fs, read.csv, stringsAsFactors=FALSE))
      cat("  labels in the BI files:", paste(unique(raw$label), collapse=", "), "\n")
      lab_age <- function(l) { l <- tolower(l)
        if (!is.null(BI_LABEL_MAP)) return(unname(BI_LABEL_MAP[l]))
        ifelse(grepl("mid|^m|meemhc|mhcyhc", l), "middle", ifelse(grepl("old|^o", l), "old", NA_character_)) }
      bi_orig <- raw %>% dplyr::transmute(age=lab_age(label), sex=tolower(sex), celltype, BI_raw=bi) %>%
        dplyr::filter(!is.na(age)) %>% dplyr::distinct(age, sex, celltype, .keep_all=TRUE)
      cat("  mapped:", paste(unique(paste(raw$label, "->", lab_age(raw$label))), collapse="; "), "\n")
    }
  }
  if (is.null(bi_orig)) stop("BI_SOURCE must be \"original\": the BI is read from the DESeq2 results (step bi)")
  bi <- if (!is.null(bi_orig)) bi_orig else dplyr::bind_rows(lapply(seq_len(nrow(cosd)), function(i) {
    ag <- cosd$age[i]; sx <- cosd$sex[i]; ct <- cosd$celltype[i]
    m <- pb_meta[pb_meta$celltype == ct & pb_meta$sex_norm == sx, ]
    y <- m$sample[m$age_norm == "young" & m$condition == "HC"]
    o <- m$sample[m$age_norm == (if (ag == "middle") MID else "old") & m$condition == "HC"]
    e <- m$sample[m$age_norm == (if (ag == "middle") MID else "old") & m$condition == "EE"]
    if (min(length(y), length(o), length(e)) < 2) return(NULL)
    data.frame(age=ag, sex=sx, celltype=ct, BI_raw=stats3(log_norm(pb_all[, c(y, o, e), drop=FALSE]), y, o, e)[["BI_raw"]])
  }))
  q <- cosd %>% dplyr::inner_join(bi, by=c("age","sex","celltype")) %>%
    dplyr::mutate(Age=factor(AGE_LAB[age], levels=AGE_LAB), celltype=factor(celltype, levels=CT_ORDER),
                  gate=factor(ifelse(pass, "passes", "does not pass"), levels=c("passes","does not pass")))
  wsrc(q %>% dplyr::select(sex, age, celltype, cos_heldout=cos, cos_dc, BI_raw, pass), "Fig5e_Quadrant_NE")
  pq <- function(d, title) quad_base(d, ymax_of(d$BI_raw)) +
      geom_point(aes(fill=celltype, alpha=gate), shape=21, color="white", size=2.2, stroke=0.2) + lab_layer(d) +
      scale_fill_manual(values=CT_COLS, guide="none") + scale_color_manual(values=CT_COLS, guide="none") +
      scale_alpha_manual(values=c(passes=1, `does not pass`=0.35), name="Buffering gate", drop=FALSE) +
      facet_wrap(~Age, nrow=1) +
      labs(title=title, x="Held-out cosine similarity (direction)", y="Buffering index (magnitude)")
  pl <- lapply(SEX_ORDER, function(s) { d <- q %>% dplyr::filter(Sex == s); if (nrow(d)) pq(d, paste("NE,", s)) else NULL })
  sv(side_by_side(pl), "Fig5e1_Quadrant_NE", 230, 85, note="each sex x age its own axes")
  for (s in names(SEX_LAB)) if (any(q$sex == s)) sv(pq(q %>% dplyr::filter(sex == s), paste("NE,", SEX_LAB[[s]])),
                                                    paste0("Fig5e1_Quadrant_NE_", s), 120, 75)
  if (!is.null(dr) && nrow(dr) && "BI_raw" %in% names(dr)) {
    ref <- q %>% dplyr::filter(sex == "female", age == "old")
    pl <- lapply(DR_APPROACHES, function(ap) {
      d <- dr %>% dplyr::filter(approach == ap, !is.na(BI_raw)) %>%
        dplyr::mutate(celltype=factor(celltype, levels=CT_ORDER), gate=factor(ifelse(pass, "passes", "does not pass"), levels=c("passes","does not pass")))
      if (!nrow(d)) return(NULL)
      quad_base(d, ymax_of(c(d$BI_raw, ref$BI_raw))) +
        geom_point(data=ref, aes(cos, BI_raw), shape=21, fill="grey85", color="grey60", size=1.5, stroke=0.2) +
        geom_point(aes(fill=celltype, alpha=gate), shape=23, color="white", size=2.3, stroke=0.2) + lab_layer(d) +
        scale_fill_manual(values=CT_COLS, guide="none") + scale_color_manual(values=CT_COLS, guide="none") +
        scale_alpha_manual(values=c(passes=1, `does not pass`=0.35), name="DREADD gate", drop=FALSE) +
        labs(title=paste0("DREADD ", ap, ", old females"), subtitle="Grey = NE, old females",
             x="Held-out cosine similarity (direction)", y="Buffering index (magnitude)") })
    if (any(!sapply(pl, is.null))) sv(side_by_side(pl), "Fig5e2_Quadrant_DREADD", 160, 85)
    wsrc(dr %>% dplyr::select(approach, celltype, cos_heldout=cos, cos_dc, BI_raw, pass), "Fig5e_Quadrant_DREADD")
  }
}

write.csv(dplyr::bind_rows(index), file.path(OUT, "Index.csv"), row.names=FALSE)
cat("\n=== DONE ===\nOutputs in:", OUT, "\n")

invisible(NULL)
}

################################################################################
### STEP: regulon  (TF regulon heatmaps for the top candidate TFs of every category; new, after TF_Regulon_Heatmap_EE.R)
###   Candidates: the Figure 6 rule (TF_Atlas_Ranking_withFoldConsistency.csv): per sex x category, top REG_N_CAND by
###     atlas score excluding clock TFs, then robust across methods + expressed + fold-consistent in >= 1 cell type;
###     the top REG_TOP of those are drawn. Categories: Buffered middle / old, Maintained, Lost, Age-emergent,
###     Amplified middle / old (REG_CATS).
###   Regulon: CollecTRI mouse (the network the TF scores use; protein-ID entries removed), not DoRothEA.
###   Columns (log2FC, DESeq2_Results): Aging middle (MHCYHC) | NE middle (MEEMHC) | Aging old (OHCYHC) |
###     NE old (OEEOHC). Buffering reads as NE opposite in colour to aging at that age.
###   Rows: regulon targets in each cell type where the TF is a candidate (up to REG_MAX_CT, by atlas evidence);
###     targets in the category's gene set first (filled dot), then other targets by |aging log2FC|, up to REG_MAX_GENES.
###     Gene-set membership from master Step2_GeneCalls.RDS (T3 buffered / A3 amplified calls, as in the master step).
###     +/- = CollecTRI mode of regulation (activates / represses).
###   Outputs: Paper_Figures/Figure5/Regulon_Heatmaps/
###     Regulon_<sex>_<category>_<TF>          one heatmap per TF (facets = cell types)
###     RegulonSummary_<sex>_<category>        all top TFs together: mean log2FC of the set-member targets (TF x cell type
###                                            rows), each target signed by its aging change, so NE < 0 = opposes aging
###     Regulon_TopTFs.csv, Regulon_Genes_<sex>_<category>.csv   what is drawn
###   Existing files are skipped.
################################################################################

STEPS$regulon <- function() {
  suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
  REG_CATS      <- c("Buffered, middle","Buffered, old","Maintained","Lost","Age-emergent","Amplified, middle","Amplified, old")
  REG_N_CAND    <- 25        # candidate rule: top N by atlas score before the robust / expressed / fold filters (as Figure 6)
  REG_TOP       <- 5         # TFs drawn per sex x category
  REG_MAX_CT    <- 4         # cell types per TF (where it is a candidate, by evidence)
  REG_MAX_GENES <- 50        # regulon rows per cell type
  LFC_LIM       <- 1.5       # colour scale clipped at +/- this log2FC
  CLOCK_TFS <- c("Bmal1","Bmal2","Arntl","Arntl2","Clock","Npas2","Nr1d1","Nr1d2","Dbp","Tef","Hlf","Nfil3",
                 "Bhlhe40","Bhlhe41","Rora","Rorb","Rorc","Cry1","Cry2","Per1","Per2","Per3","Arnt")
  CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN","Astrocytes","Oligodendrocytes","OPC",
                "Microglia","Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell")
  COLS <- c(MHCYHC="Aging\nmiddle", MEEMHC="NE\nmiddle", OHCYHC="Aging\nold", OEEOHC="NE\nold")
  CAT_AGE <- c(`Buffered, middle`="middle", `Buffered, old`="old", Maintained="middle", Lost="middle",
               `Age-emergent`="old", `Amplified, middle`="middle", `Amplified, old`="old")   # aging column used to sort
  OUTR <- file.path(NE_BASE, "Paper_Figures", "Figure5", "Regulon_Heatmaps"); dir.create(OUTR, recursive=TRUE, showWarnings=FALSE)
  slug <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  sv2 <- function(p, stem, w, h) for (ext in c("pdf", "png")) {
    f <- file.path(OUTR, paste0(stem, ".", ext))
    if (file.exists(f)) { cat("  exists, skipped:", basename(f), "\n"); next }
    if (ext == "pdf") ggplot2::ggsave(f, p, width=min(w, 45), height=min(h, 45), limitsize=FALSE)
    else ggplot2::ggsave(f, p, width=min(w, 45), height=min(h, 45), dpi=300, bg="white", limitsize=FALSE)
    cat("  saved:", basename(f), "\n") }
  wr <- function(x, name) { f <- file.path(OUTR, name); if (!file.exists(f)) write.csv(x, f, row.names=FALSE) }

  # ---- candidates (Figure 6 rule) ----
  FOLD <- file.path(MASTER_DIR, "TF_FoldConsistency")
  rank_tab <- read.csv(file.path(FOLD, "TF_Atlas_Ranking_withFoldConsistency.csv"), stringsAsFactors=FALSE)
  long     <- read.csv(file.path(FOLD, "TF_Atlas_Long_withFoldConsistency.csv"), stringsAsFactors=FALSE)
  cand <- rank_tab %>% dplyr::filter(category %in% REG_CATS, !source %in% CLOCK_TFS) %>% dplyr::group_by(sex, category) %>%
    dplyr::arrange(dplyr::desc(score), .by_group=TRUE) %>% dplyr::slice_head(n=REG_N_CAND) %>%
    dplyr::filter(n_ct_robust >= 1, n_ct_expressed >= 1, n_ct_fold_consistent >= 1) %>%
    dplyr::slice_head(n=REG_TOP) %>% dplyr::mutate(rank_drawn=dplyr::row_number()) %>% dplyr::ungroup()
  cts <- long %>% dplyr::semi_join(cand, by=c("sex","category","source")) %>% dplyr::filter(celltype %in% CT_ORDER) %>%
    dplyr::group_by(sex, category, source) %>% dplyr::slice_max(evidence, n=REG_MAX_CT, with_ties=FALSE) %>% dplyr::ungroup()
  wr(cand %>% dplyr::left_join(cts %>% dplyr::group_by(sex, category, source) %>%
       dplyr::summarise(celltypes_drawn=paste(celltype, collapse=", "), .groups="drop"), by=c("sex","category","source")),
     "Regulon_TopTFs.csv")
  cat("Top TFs drawn per sex x category:\n")
  print(as.data.frame(cand %>% dplyr::group_by(sex, category) %>% dplyr::summarise(TFs=paste(source, collapse=", "), .groups="drop")), row.names=FALSE)

  # ---- regulons, log2FC, set membership ----
  uniprot <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
  net <- readRDS(COLLECTRI_RDS) %>%
    dplyr::filter(!grepl(uniprot, source), !grepl(uniprot, target)) %>% dplyr::distinct(source, target, .keep_all=TRUE)
  lfc <- dplyr::bind_rows(lapply(c("male","female"), function(sx) dplyr::bind_rows(lapply(names(COLS), function(cp)
    read.csv(file.path(DESEQ_DIR, paste0(sx, "_ALL_", cp, ".csv")), stringsAsFactors=FALSE) %>%
      dplyr::transmute(sex=sx, comp=cp, celltype, gene, lfc=log2FoldChange, padj)))))
  calls <- readRDS(file.path(MASTER_DIR, "Step2_GeneCalls.RDS")) %>% dplyr::select(age, sex, celltype, gene, T3, A3)
  cw <- calls %>% tidyr::pivot_wider(names_from=age, values_from=c(T3, A3), values_fill=FALSE)
  for (cc in c("T3_middle","T3_old","A3_middle","A3_old")) if (!cc %in% names(cw)) cw[[cc]] <- FALSE
  member_of <- function(cat, d) switch(cat,
    `Buffered, middle`=d$T3_middle, `Buffered, old`=d$T3_old, Maintained=d$T3_middle & d$T3_old,
    Lost=d$T3_middle & !d$T3_old, `Age-emergent`=!d$T3_middle & d$T3_old,
    `Amplified, middle`=d$A3_middle, `Amplified, old`=d$A3_old)

  for (sx in unique(cand$sex)) for (cat in intersect(REG_CATS, unique(cand$category))) {
    tfs <- cand %>% dplyr::filter(sex == sx, category == cat) %>% dplyr::pull(source)
    if (!length(tfs)) next
    sort_comp <- if (CAT_AGE[[cat]] == "middle") "MHCYHC" else "OHCYHC"
    rows_all <- list(); summ <- list()
    for (tf in tfs) {
      ctv <- cts %>% dplyr::filter(sex == sx, category == cat, source == tf) %>% dplyr::pull(celltype)
      targ <- net %>% dplyr::filter(source == tf) %>% dplyr::select(gene=target, mor)
      if (!length(ctv) || !nrow(targ)) next
      d <- lfc %>% dplyr::filter(sex == sx, celltype %in% ctv, gene %in% targ$gene) %>%
        dplyr::select(celltype, gene, comp, lfc) %>% tidyr::pivot_wider(names_from=comp, values_from=lfc) %>%
        dplyr::left_join(targ, by="gene") %>%
        dplyr::left_join(cw %>% dplyr::filter(sex == sx) %>% dplyr::select(-sex), by=c("celltype","gene"))
      for (cc in c("T3_middle","T3_old","A3_middle","A3_old")) d[[cc]] <- d[[cc]] %in% TRUE
      for (cp in names(COLS)) if (!cp %in% names(d)) d[[cp]] <- NA_real_
      d$member <- member_of(cat, d)
      d <- d %>% dplyr::filter(!is.na(.data[[sort_comp]])) %>% dplyr::group_by(celltype) %>%
        dplyr::arrange(dplyr::desc(member), dplyr::desc(ifelse(member, .data[[sort_comp]], abs(.data[[sort_comp]]))), .by_group=TRUE) %>%
        dplyr::slice_head(n=REG_MAX_GENES) %>% dplyr::ungroup() %>% dplyr::mutate(tf=tf)
      if (!nrow(d)) next
      rows_all[[tf]] <- d
      # aging-aligned: each target's values multiplied by the sign of its aging log2FC at the category's age, so aging
      # is positive and an NE column below 0 means NE opposes aging (up and down targets no longer cancel)
      summ[[tf]] <- d %>% dplyr::filter(member) %>% dplyr::mutate(s_age=sign(.data[[sort_comp]])) %>%
        dplyr::group_by(tf, celltype) %>%
        dplyr::summarise(n_members=dplyr::n(), dplyr::across(dplyr::all_of(names(COLS)), ~ mean(.x * s_age, na.rm=TRUE)), .groups="drop")

      # one heatmap per TF: facets = cell types
      dl <- d %>% dplyr::mutate(celltype=factor(celltype, levels=intersect(CT_ORDER, ctv)),
                                row=paste(gene, celltype, sep="__")) %>%
        dplyr::group_by(celltype) %>% dplyr::mutate(ord=dplyr::row_number()) %>% dplyr::ungroup()
      hl <- dl %>% tidyr::pivot_longer(dplyr::all_of(names(COLS)), names_to="comp", values_to="val") %>%
        dplyr::mutate(x=match(comp, names(COLS)))
      side <- dl %>% dplyr::mutate(x_m=0.2, x_r=-0.35)
      n_max <- max(table(dl$celltype))
      p <- ggplot() +
        geom_tile(data=hl, aes(x=x, y=-ord, fill=pmax(pmin(val, LFC_LIM), -LFC_LIM)), color="white", linewidth=0.15) +
        geom_vline(xintercept=2.5, color="grey30", linewidth=0.4) +
        geom_point(data=side[side$member, ], aes(x=x_m, y=-ord), size=1.1, color="black") +
        geom_text(data=side, aes(x=x_r, y=-ord, label=ifelse(mor > 0, "+", "-")), size=2, color="grey35") +
        geom_text(data=dl, aes(x=4.65, y=-ord, label=gene), hjust=0, size=if (n_max > 35) 1.5 else 2, fontface="italic") +
        scale_fill_gradient2(low="#2166AC", mid="grey97", high="#D6604D", midpoint=0, limits=c(-LFC_LIM, LFC_LIM),
                             na.value="grey85", name=paste0("log2FC\n(+/-", LFC_LIM, ")")) +
        scale_x_continuous(breaks=1:4, labels=COLS, expand=expansion(add=c(0.6, 1.6))) +
        facet_wrap(~celltype, nrow=1, scales="free_y") +
        theme_classic(base_size=8) +
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank(),
              axis.text.x=element_text(color="black", size=6.5), axis.title=element_blank(),
              strip.background=element_blank(), strip.text=element_text(face="bold"),
              plot.title=element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5, size=6.5, color="grey35")) +
        labs(title=paste0(tf, " regulon: ", cat, " (", if (sx == "male") "Males" else "Females", ")"),
             subtitle=paste0("CollecTRI targets in cell types where ", tf, " is a candidate; dot = in the ", tolower(cat),
                             " gene set (listed first); +/- = activates / represses\nBuffering: NE opposite in colour to aging at that age"))
      sv2(p, paste0("Regulon_", sx, "_", slug(cat), "_", slug(tf)), 2.2 * length(ctv) + 1.6, max(3, n_max * 0.11 + 1.6))
    }
    if (!length(rows_all)) next
    wr(dplyr::bind_rows(rows_all), paste0("Regulon_Genes_", sx, "_", slug(cat), ".csv"))

    # summary: TF x cell type rows, mean log2FC of the set-member targets
    s <- dplyr::bind_rows(summ)
    if (nrow(s)) {
      s <- s %>% dplyr::mutate(row=paste0(tf, " | ", celltype, " (n=", n_members, ")"),
                               tf=factor(tf, levels=tfs), celltype=factor(celltype, levels=CT_ORDER)) %>%
        dplyr::arrange(tf, celltype) %>% dplyr::mutate(row=factor(row, levels=rev(unique(row))))
      sl <- s %>% tidyr::pivot_longer(dplyr::all_of(names(COLS)), names_to="comp", values_to="val") %>%
        dplyr::mutate(comp=factor(COLS[comp], levels=COLS))
      p2 <- ggplot(sl, aes(comp, row, fill=pmax(pmin(val, LFC_LIM), -LFC_LIM))) +
        geom_tile(color="white", linewidth=0.3) +
        geom_text(aes(label=sprintf("%.2f", val)), size=1.8) +
        geom_vline(xintercept=2.5, color="grey30", linewidth=0.4) +
        scale_fill_gradient2(low="#2166AC", mid="grey97", high="#D6604D", midpoint=0, limits=c(-LFC_LIM, LFC_LIM),
                             na.value="grey85", name="Mean aligned\nlog2FC") +
        theme_minimal(base_size=8) +
        theme(panel.grid=element_blank(), axis.text=element_text(color="black"), axis.title=element_blank(),
              plot.title=element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5, size=6.5, color="grey35")) +
        labs(title=paste0("Top ", cat, " TFs (", if (sx == "male") "Males" else "Females", "): regulon targets in the gene set"),
             subtitle=paste0("Mean log2FC of each TF's gene-set targets (n), each signed by its aging change at ", CAT_AGE[[cat]],
                             " age\naging > 0 by construction; NE < 0 (blue) = NE opposes aging"))
      sv2(p2, paste0("RegulonSummary_", sx, "_", slug(cat)), 5.2, max(2.5, 0.16 * nlevels(s$row) + 1.4))
    }
  }
  cat("Regulon heatmaps in:", OUTR, "\n")
  invisible(NULL)
}

################################################################################
### Run
################################################################################

STEPS$bi <- step_bi
ORDER <- c("robustness", "master", "methodcheck", "atlas", "fold", "bi", "figure5", "regulon")
for (st in ORDER) {
  if (!st %in% c("figure5", "regulon")) {
    status <- CHECKS[[st]]()
    if (status == "ok" && !(st %in% FORCE)) { say(sprintf("STEP %s: up to date -> skipped", st)); next }
    say(sprintf("STEP %s: %s -> running", st, if (status == "ok") "forced" else status))
    if (grepl("out of date", status) || st %in% FORCE) if (!is.null(OWNS[[st]])) supersede(OWNS[[st]])
  } else say("STEP ", st, ": running")
  t0 <- Sys.time()
  STEPS[[st]]()                                  # each step runs in its own environment (a function call)
  invisible(gc())
  say(sprintf("STEP %s finished in %.1f min", st, as.numeric(difftime(Sys.time(), t0, units="mins"))))
  if (!st %in% c("figure5", "regulon")) { status <- CHECKS[[st]](); if (status != "ok") stop(sprintf("STEP %s ran but its output is not right: %s", st, status)) }
  if (!is.null(STOP_AFTER) && st == STOP_AFTER) { say("Stopped after ", st, " (STOP_AFTER)"); break }
}
say("DONE. Figure 5 in ", file.path(NE_BASE, "Paper_Figures", "Figure5", TODAY))
