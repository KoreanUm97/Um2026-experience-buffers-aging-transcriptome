
################################################################################
### Buffering_DREADD_vs_NE.R
###
### Does DREADD activation (old C21 vs old PBS) reproduce NE buffering in old females?
### Approaches B and D (mScarlet >= 50 neurons; glia all cells), female only,
### 3 mice per group (YPBS, OPBS, OC21; YC21 not used).
###
### Held-out fits (k = each OPBS mouse):  B_k = OPBS(-k) vs YPBS,   E_k = OC21 vs OPBS(k)
###   Approach B neurons: ~ mScarlet_median_scaled + grp (as in DREADD_DESeq2.R).
### Every labeling of the 6 old mice (C21 vs PBS; 20 labelings, 1 real) is scored, so each
### result has a rank among its own 20 labelings (1 = most extreme).
###
###   PART 1  DREADD cosine per cell type (held-out; informational)
###   PART 2  DREADD gene calls on DREADD's own aging axis (T1/T3 buffered, A1/A3 amplified)
###           + are there more buffered / amplified genes than in relabeled mice?
###   PART 3  Gene overlap with NE (old female; NE gene calls from Buffering_Pipeline.R step 2):
###           among genes NE calls aging-changed in that cell type,
###             do NE-BUFFERED genes get moved OPPOSITE to NE's aging direction by DREADD
###             more often than other NE aging genes?   (OR, Fisher; rank vs relabelings)
###             do NE-AMPLIFIED genes get moved WITH NE's aging direction by DREADD more often?
###           pooled across cell types (Mantel-Haenszel; relabel rank on summed log OR)
###   PART 4  TF activity (decoupleR ULM, signed CollecTRI, held-out Wald stats):
###             S_own   = mean_k(-a_aging,DREADD * a_C21)   DREADD on its own aging axis
###             S_NEax  = mean_k(-a_aging,NE * a_C21)       DREADD against NE's aging direction
###           per cell type + summed across cell types; ranks vs relabelings;
###           side-by-side with NE TF activity (old female) for the shortlist
###
### Inputs: DREADD pseudobulk caches from DREADD_Buffering_HeldOut.R
###   (pseudobulk_mScarlet_50_Approach_{B,D}.RDS); if missing, df_harmony (DREADD) must be loaded.
### NE inputs: Buffering_Candidates_Master/ (Step2_GeneCalls.RDS, Step5_TFActivity_PerCellType.csv).
### Resumable: per approach x cell type results are cached.
################################################################################

Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages({
  library(conflicted); library(DESeq2); library(dplyr); library(tidyr); library(ggplot2)
  library(parallel); library(decoupleR)
  conflicts_prefer(dplyr::select, .quiet=TRUE); conflicts_prefer(dplyr::filter, .quiet=TRUE)
  conflicts_prefer(dplyr::mutate, .quiet=TRUE); conflicts_prefer(base::intersect, .quiet=TRUE)
  conflicts_prefer(base::setdiff, .quiet=TRUE); conflicts_prefer(base::union, .quiet=TRUE)
  conflicts_prefer(base::unique, .quiet=TRUE); conflicts_prefer(base::Filter, .quiet=TRUE)
  conflicts_prefer(base::setequal, .quiet=TRUE); conflicts_prefer(base::as.data.frame, .quiet=TRUE)
})
if (requireNamespace("RhpcBLASctl", quietly=TRUE)) RhpcBLASctl::blas_set_num_threads(1)

################################################################################
### PARAMETERS
################################################################################

DR_PATH <- "PATH_TO_DREADD_OUTPUT"
NE_PATH <- "PATH_TO_NE_OUTPUT"
NE_MASTER <- file.path(NE_PATH, "Buffering_Candidates_Master")
COLLECTRI_RDS <- "PATH_TO_COLLECTRI/collectri_mouse.rds"
PB_CACHE_DIR  <- file.path(DR_PATH, "DREADD_Buffering_HeldOut")

ARG_FILE   <- NULL              # path to an ARG gene list (.csv with a 'gene' column, or .txt one gene per line); NULL = skip
APPROACHES <- c("B","D")
MSCARLET_THR <- 50
N_CORES    <- 32
P_AGING <- 0.05; P_NE <- 0.05
MIN_REGULON <- 5; ACT_RANK <- 5
CELLTYPES_TO_RUN <- NULL       # NULL = all; e.g. c("CA1 ExNs") for a quick test
TRIM_LO <- 0.10; TRIM_HI <- 0.90; MIN_CELLS_FLOOR <- 5   # Approach D (as in DREADD_Buffering_HeldOut.R)
SHORTLIST <- c("Nr3c1","Creb1","Sirt1","Atf6","Ddit3","Atf4","Foxo1","Foxo3","Nfkb1","Nfkb2","Rela","Relb",
               "Egr1","Egr3","Fos","Jun","Hsf1")

out_root <- file.path(DR_PATH, "Buffering_DREADD_vs_NE")
dir.create(out_root, recursive=TRUE, showWarnings=FALSE)
IN_SUBTYPES <- c("IN1","IN2")
NEURONAL_CELLTYPES <- c("CA1 ExNs","CA2 ExNs","CA3 ExNs","DG ExNs","Sub ExNs","IN")
CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN",
              "Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell","Macrophage","Epithelial","T Cell")
th <- theme_classic(base_size=8) +
  theme(strip.background=element_rect(fill="grey92", color=NA), strip.text=element_text(face="bold"),
        plot.title=element_text(face="bold", hjust=0.5, size=9), plot.subtitle=element_text(color="grey40", hjust=0.5, size=7.5))
sv <- function(p, dir, stem, w, h) {
  ggsave(file.path(dir, paste0(stem, ".pdf")), p, width=min(w, 40), height=min(h, 40), limitsize=FALSE)
  ggsave(file.path(dir, paste0(stem, ".png")), p, width=min(w, 40), height=min(h, 40), dpi=200, bg="white", limitsize=FALSE)
}
wcsv <- function(x, dir, name) write.csv(x, file.path(dir, name), row.names=FALSE)
banner <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep="")

################################################################################
### NE REFERENCE (old female)
################################################################################

banner("NE reference (old female)")
ne_calls <- readRDS(file.path(NE_MASTER, "Step2_GeneCalls.RDS")) %>%
  dplyr::filter(sex == "female", age == "old") %>%
  dplyr::transmute(celltype, gene, ne_aging_lfc=aging_lfc, ne_agechg = n_aging >= 2, ne_T3=T3, ne_A3=A3)
ARG_GENES <- character(0)
if (!is.null(ARG_FILE)) {
  ARG_GENES <- if (grepl("\\.csv$", ARG_FILE, ignore.case=TRUE)) {
    a <- read.csv(ARG_FILE, stringsAsFactors=FALSE); if ("gene" %in% names(a)) a$gene else a[[1]]
  } else readLines(ARG_FILE)
  ARG_GENES <- unique(trimws(ARG_GENES[nzchar(trimws(ARG_GENES))]))
  cat("ARG genes loaded:", length(ARG_GENES), "\n")
}
ne_act <- read.csv(file.path(NE_MASTER, "Step5_TFActivity_PerCellType.csv"), stringsAsFactors=FALSE) %>%
  dplyr::filter(sex == "female", age == "old") %>%
  dplyr::transmute(celltype, source, NE_S=S, NE_aging_act=aging_act, NE_rank_buf=rank_buf, NE_rank_amp=rank_amp,
                   NE_call=dplyr::case_when(buf ~ "buffered", amp ~ "amplified", TRUE ~ "-"))
cat(sprintf("NE gene calls: %d cell types; NE TF activity: %d cell types\n",
            dplyr::n_distinct(ne_calls$celltype), dplyr::n_distinct(ne_act$celltype)))

uniprot <- "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
net_all <- readRDS(COLLECTRI_RDS) %>% dplyr::filter(!grepl(uniprot, source), !grepl(uniprot, target)) %>%
  dplyr::distinct(source, target, .keep_all=TRUE)

################################################################################
### DREADD PSEUDOBULK (cache from DREADD_Buffering_HeldOut.R, else build)
################################################################################

get_pb <- function(ap) {
  f <- file.path(PB_CACHE_DIR, paste0("pseudobulk_mScarlet_", MSCARLET_THR, "_Approach_", ap, ".RDS"))
  if (file.exists(f)) { cat("Pseudobulk cache:", f, "\n"); return(readRDS(f)) }
  stopifnot("Pseudobulk cache missing: load the DREADD df_harmony and rerun" = exists("df_harmony"))
  suppressPackageStartupMessages({ library(Seurat); library(BPCells) })
  md <- df_harmony@meta.data
  inp <- intersect(IN_SUBTYPES, unique(md$celltype)); if (length(inp)) md$celltype[md$celltype %in% inp] <- "IN"
  if (!"mScarlet_raw" %in% colnames(md)) md$mScarlet_raw <- FetchData(df_harmony, vars="mScarlet", layer="counts")$mScarlet
  md <- md %>% dplyr::mutate(age_norm=tolower(age), condition_norm=experimental.condition,
                             age_condition=paste0(age_norm, "_", experimental.condition))
  if ("cell_id" %in% colnames(md)) md <- md %>% dplyr::select(-cell_id)
  meta_full <- md %>% tibble::rownames_to_column("cell_id") %>%
    dplyr::select(cell_id, celltype, orig.ident, age_norm, condition_norm, age_condition, mScarlet_raw)
  base <- meta_full %>% dplyr::filter((celltype %in% NEURONAL_CELLTYPES & mScarlet_raw >= MSCARLET_THR) |
                                      (!celltype %in% NEURONAL_CELLTYPES))
  cells <- if (ap == "B") base else {
    trimmed <- dplyr::bind_rows(lapply(unique(base$celltype), function(ct) {
      cm <- base %>% dplyr::filter(celltype == ct); if (!ct %in% NEURONAL_CELLTYPES) return(cm)
      bnd <- dplyr::bind_rows(lapply(unique(cm$orig.ident), function(m) { v <- cm$mScarlet_raw[cm$orig.ident == m]
        if (length(v) < 5) NULL else data.frame(lo=quantile(v, TRIM_LO), hi=quantile(v, TRIM_HI)) }))
      if (!nrow(bnd)) return(cm); lo <- median(bnd$lo); hi <- median(bnd$hi); if (lo >= hi) return(cm)
      cm %>% dplyr::filter(mScarlet_raw >= lo & mScarlet_raw <= hi) }))
    set.seed(42)
    tgt <- trimmed %>% dplyr::count(celltype, orig.ident, age_condition) %>% dplyr::group_by(celltype) %>%
      dplyr::summarise(target_n=min(n), .groups="drop") %>% dplyr::filter(target_n >= MIN_CELLS_FLOOR)
    dplyr::bind_rows(lapply(tgt$celltype, function(ct) { t <- tgt$target_n[tgt$celltype == ct]
      trimmed %>% dplyr::filter(celltype == ct) %>% dplyr::group_by(orig.ident, age_condition) %>%
        dplyr::group_modify(~dplyr::slice_sample(.x, n=min(t, nrow(.x)))) %>% dplyr::ungroup() }))
  }
  so <- df_harmony[, cells$cell_id]; idx <- match(colnames(so), cells$cell_id)
  so@meta.data$age_norm <- cells$age_norm[idx]; so@meta.data$condition_norm <- cells$condition_norm[idx]
  pb <- Seurat::AggregateExpression(object=so, assays="RNA", return.seurat=FALSE,
          group.by=c("celltype","orig.ident","condition_norm","age_norm"), normalization.method=NULL,
          slot="counts", verbose=FALSE)$RNA
  parse <- function(s) { ctn <- CT_ORDER[sapply(CT_ORDER, function(ct) startsWith(s, paste0(ct, "_")))]
    if (!length(ctn)) return(c(NA, NA, NA, NA)); ctn <- ctn[1]
    p <- strsplit(sub(paste0("^", ctn, "_"), "", s), "_")[[1]]; n <- length(p)
    c(ctn, gsub("-", "_", paste(p[1:(n-2)], collapse="_")), p[n-1], p[n]) }
  pp <- t(sapply(colnames(pb), parse))
  meta <- data.frame(sample=colnames(pb), celltype=pp[,1], orig.ident=pp[,2], condition=pp[,3], age=pp[,4],
                     stringsAsFactors=FALSE)
  keep <- !is.na(meta$celltype)
  PB <- list(counts=as.matrix(pb[, keep, drop=FALSE]), meta=meta[keep, ],
             mcov=cells %>% dplyr::filter(celltype %in% NEURONAL_CELLTYPES) %>% dplyr::group_by(orig.ident, celltype) %>%
               dplyr::summarise(mScarlet_median=median(mScarlet_raw), .groups="drop"))
  saveRDS(PB, f); PB
}

################################################################################
### FITS, GENE CALLS, TF ACTIVITY PER LABELING
################################################################################

fit_one <- function(PB, num, den, ct, ap) {
  cols <- c(den, num)
  cnt  <- round(PB$counts[, cols, drop=FALSE]); cnt <- cnt[rowSums(cnt >= 10) >= 2, , drop=FALSE]
  if (nrow(cnt) < 10) return(NULL)
  cd <- data.frame(grp=factor(rep(c("den","num"), c(length(den), length(num))), levels=c("den","num")), row.names=cols)
  design <- ~grp
  if (ap == "B" && ct %in% NEURONAL_CELLTYPES && !is.null(PB$mcov) && length(cols) >= 4) {
    oid <- PB$meta$orig.ident[match(cols, PB$meta$sample)]
    mv  <- PB$mcov$mScarlet_median[match(paste(oid, ct), paste(PB$mcov$orig.ident, PB$mcov$celltype))]
    if (!any(is.na(mv)) && var(mv) > 0) { cd$mScarlet_median_scaled <- as.numeric(scale(mv)); design <- ~mScarlet_median_scaled + grp }
  }
  tryCatch({
    dds <- suppressMessages(DESeqDataSetFromMatrix(cnt, cd, design=design))
    dds <- tryCatch(DESeq(dds, quiet=TRUE), error=function(e) {
      d <- estimateSizeFactors(dds); d <- estimateDispersionsGeneEst(d); dispersions(d) <- mcols(d)$dispGeneEst
      nbinomWaldTest(d, quiet=TRUE) })
    r <- results(dds, name="grp_num_vs_den")
    data.frame(gene=rownames(r), lfc=r$log2FoldChange, stat=r$stat, pvalue=r$pvalue, stringsAsFactors=FALSE)
  }, error=function(e) NULL)
}

tf_act <- function(B, E) {     # B, E: data.frames with gene, stat
  if (is.null(B) || is.null(E)) return(NULL)
  b <- setNames(B$stat, B$gene); e <- setNames(E$stat, E$gene)
  g <- intersect(names(b)[!is.na(b)], names(e)[!is.na(e)]); if (length(g) < 500) return(NULL)
  mat <- cbind(B=b[g], E=e[g]); rownames(mat) <- g
  r <- tryCatch(decoupleR::run_ulm(mat=mat, net=net_all %>% dplyr::filter(target %in% g), .source="source",
                                   .target="target", .mor="mor", minsize=MIN_REGULON), error=function(e) NULL)
  if (is.null(r)) return(NULL)
  r %>% dplyr::filter(statistic=="ulm") %>% dplyr::select(source, condition, score) %>%
    tidyr::pivot_wider(names_from=condition, values_from=score, names_prefix="a_")
}

score_labeling <- function(PB, ct, ap, ypbs, pbs, c21, ne_g, ne_a) {
  fits <- lapply(seq_along(pbs), function(k) list(B=fit_one(PB, pbs[-k], ypbs, ct, ap), E=fit_one(PB, c21, pbs[k], ct, ap)))
  fits <- Filter(function(f) !is.null(f$B) && !is.null(f$E), fits); if (length(fits) < 2) return(NULL)
  # per-gene, per-fit table
  pf <- dplyr::bind_rows(lapply(seq_along(fits), function(k)
    dplyr::inner_join(fits[[k]]$B %>% dplyr::select(gene, lfc_B=lfc, p_B=pvalue),
                      fits[[k]]$E %>% dplyr::select(gene, lfc_E=lfc, p_E=pvalue), by="gene") %>%
      dplyr::filter(!is.na(lfc_B), !is.na(lfc_E), lfc_B != 0) %>% dplyr::mutate(k=k)))
  calls <- pf %>% dplyr::left_join(ne_g %>% dplyr::select(gene, ne_aging_lfc), by="gene") %>%
    dplyr::mutate(ag=!is.na(p_B) & p_B < P_AGING, nes=!is.na(p_E) & p_E < P_NE,
                  opp=ag & sign(lfc_E) != sign(lfc_B), same=ag & sign(lfc_E) == sign(lfc_B),
                  opp_neax = sign(lfc_E) != sign(ne_aging_lfc), same_neax = sign(lfc_E) == sign(ne_aging_lfc)) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(n_fits=dplyr::n(), n_aging=sum(ag), n_opp=sum(opp), n_same=sum(same),
                     n_opp_nesig=sum(opp & nes), n_same_nesig=sum(same & nes),
                     n_opp_neax=sum(opp_neax, na.rm=TRUE), n_same_neax=sum(same_neax, na.rm=TRUE),
                     aging_lfc=mean(lfc_B), c21_lfc=mean(lfc_E), .groups="drop") %>%
    dplyr::mutate(T1=n_aging >= 3 & n_opp >= 3 & n_opp_nesig >= 2, T3=n_aging >= 2 & n_opp >= 2,
                  A1=n_aging >= 3 & n_same >= 3 & n_same_nesig >= 2, A3=n_aging >= 2 & n_same >= 2,
                  dr_opp_neax = n_opp_neax >= 2, dr_same_neax = n_same_neax >= 2)
  # overlap with NE (genes NE calls aging-changed in this cell type, tested in DREADD)
  ov <- NULL
  if (!is.null(ne_g) && nrow(ne_g)) {
    j <- calls %>% dplyr::inner_join(ne_g %>% dplyr::filter(ne_agechg), by="gene")
    or_of <- function(x, y) { a <- sum(x & y); b <- sum(x & !y); c <- sum(!x & y); d <- sum(!x & !y)
      c(a=a, b=b, c=c, d=d, lor=log(((a + .5) * (d + .5)) / ((b + .5) * (c + .5)))) }
    ob <- or_of(j$ne_T3, j$dr_opp_neax); oa <- or_of(j$ne_A3, j$dr_same_neax)
    ov <- data.frame(test=c("NE-buffered genes reversed by DREADD", "NE-amplified genes pushed further by DREADD"),
                     a=c(ob["a"], oa["a"]), b=c(ob["b"], oa["b"]), c=c(ob["c"], oa["c"]), d=c(ob["d"], oa["d"]),
                     log_or=c(ob["lor"], oa["lor"]), n_ne_aging=nrow(j))
    if (length(ARG_GENES)) {   # same tests with activity-regulated genes removed
      j2 <- j %>% dplyr::filter(!gene %in% ARG_GENES)
      ob2 <- or_of(j2$ne_T3, j2$dr_opp_neax); oa2 <- or_of(j2$ne_A3, j2$dr_same_neax)
      ov <- rbind(ov, data.frame(test=paste(ov$test, "(ARGs removed)"),
                                 a=c(ob2["a"], oa2["a"]), b=c(ob2["b"], oa2["b"]), c=c(ob2["c"], oa2["c"]), d=c(ob2["d"], oa2["d"]),
                                 log_or=c(ob2["lor"], oa2["lor"]), n_ne_aging=nrow(j2)))
    }
  }
  # TF activity
  ta <- dplyr::bind_rows(lapply(seq_along(fits), function(k) { a <- tf_act(fits[[k]]$B, fits[[k]]$E)
    if (is.null(a)) NULL else dplyr::mutate(a, k=k) }))
  tf <- if (nrow(ta)) ta %>% dplyr::left_join(ne_a %>% dplyr::select(source, NE_aging_act), by="source") %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(n_fits=dplyr::n(), S_own=mean(-a_B * a_E), S_NEax=mean(-NE_aging_act * a_E),
                     aging_act=mean(a_B), c21_act=mean(a_E), .groups="drop") %>% dplyr::filter(n_fits >= 2) else NULL
  list(calls=calls, ov=ov, tf=tf,
       counts=data.frame(n_genes=nrow(calls), T3=sum(calls$T3), A3=sum(calls$A3), T1=sum(calls$T1), A1=sum(calls$A1),
                         cos=sum(calls$aging_lfc * calls$c21_lfc) / sqrt(sum(calls$aging_lfc^2) * sum(calls$c21_lfc^2))))
}

################################################################################
### MAIN LOOP (per approach x cell type; resumable)
################################################################################

rank_cut <- function(n_lab) pmax(1, floor(ACT_RANK / 20 * n_lab))   # e.g. 5 of 20 -> 2 of 10
emp_rank <- function(x, real, upper=TRUE) { r <- x[real][1]; if (is.na(r)) return(NA_integer_)
  if (upper) sum(x >= r, na.rm=TRUE) else sum(x <= r, na.rm=TRUE) }

for (ap in APPROACHES) {
  banner(sprintf("Approach %s (mScarlet >= %d)", ap, MSCARLET_THR))
  odir <- file.path(out_root, paste0("Approach_", ap)); udir <- file.path(odir, "units")
  dir.create(udir, recursive=TRUE, showWarnings=FALSE)
  PB <- get_pb(ap)
  cts <- intersect(CT_ORDER, unique(PB$meta$celltype)); if (!is.null(CELLTYPES_TO_RUN)) cts <- intersect(CELLTYPES_TO_RUN, cts)

  for (ct in cts) {
    f_out <- file.path(udir, paste0(gsub("[^A-Za-z0-9]+", "_", ct), ".RDS"))
    if (file.exists(f_out)) { cat("  done already:", ct, "\n"); next }
    m <- PB$meta[PB$meta$celltype == ct, ]
    ypbs <- m$sample[m$age=="young" & m$condition=="PBS"]
    opbs <- m$sample[m$age=="old" & m$condition=="PBS"]; oc21 <- m$sample[m$age=="old" & m$condition=="C21"]
    if (length(ypbs) < 2 || length(opbs) < 2 || length(oc21) < 2 || length(opbs) + length(oc21) < 5) {
      cat(sprintf("  skip (n): %s  YPBS=%d OPBS=%d OC21=%d\n", ct, length(ypbs), length(opbs), length(oc21))); next }
    ne_g <- ne_calls %>% dplyr::filter(celltype == ct); ne_a <- ne_act %>% dplyr::filter(celltype == ct)
    pool <- c(opbs, oc21); mouse <- setNames(m$orig.ident[match(pool, m$sample)], pool)
    labs <- combn(pool, length(oc21), simplify=FALSE)
    cat(sprintf("  %s: YPBS=%d OPBS=%d OC21=%d; NE ref genes=%d, TFs=%d; %d labelings ... ",
                ct, length(ypbs), length(opbs), length(oc21), nrow(ne_g), nrow(ne_a), length(labs)))
    res <- mclapply(labs, function(lab) tryCatch({
      s <- score_labeling(PB, ct, ap, ypbs, setdiff(pool, lab), lab, ne_g, ne_a); if (is.null(s)) return(NULL)
      real <- setequal(lab, oc21); key <- paste(sort(mouse[lab]), collapse="|")
      list(real=real, key=key, calls=if (real) s$calls else NULL,
           ov=if (!is.null(s$ov)) dplyr::mutate(s$ov, real=real, key=key) else NULL,
           tf=if (!is.null(s$tf)) dplyr::mutate(s$tf, real=real, key=key) else NULL,
           counts=dplyr::mutate(s$counts, real=real, key=key))
    }, error=function(e) conditionMessage(e)), mc.cores=N_CORES)
    errs <- unlist(Filter(is.character, res)); if (length(errs)) { cat("\n    errors:", head(errs, 2), "\n") }
    res <- Filter(is.list, res)
    saveRDS(list(celltype=ct, calls=dplyr::bind_rows(lapply(res, `[[`, "calls")),
                 ov=dplyr::bind_rows(lapply(res, `[[`, "ov")), tf=dplyr::bind_rows(lapply(res, `[[`, "tf")),
                 counts=dplyr::bind_rows(lapply(res, `[[`, "counts"))), f_out)
    cat("saved\n")
  }

  ##############################################################################
  ### SUMMARIES (this approach)
  ##############################################################################
  U <- lapply(list.files(udir, pattern="\\.RDS$", full.names=TRUE), readRDS)
  if (!length(U)) next
  ctn <- sapply(U, `[[`, "celltype")
  n_lab <- sapply(U, function(u) nrow(u$counts)); full_cts <- ctn[n_lab == max(n_lab)]
  if (any(n_lab < max(n_lab)))
    cat("Cell types with fewer labelings (one group has 2 mice; ranks out of", paste(unique(n_lab[n_lab < max(n_lab)]), collapse="/"),
        "; excluded from pooled tests):", paste(ctn[n_lab < max(n_lab)], collapse=", "), "\n")

  # PART 1: held-out cosine (real labeling) on DREADD's own aging axis + NE axis
  cosr <- dplyr::bind_rows(lapply(U, function(u) { x <- u$calls; if (is.null(x) || !nrow(x)) return(NULL)
    data.frame(celltype=u$celltype, n_genes=nrow(x),
               cos_mean_lfc = sum(x$aging_lfc * x$c21_lfc) / sqrt(sum(x$aging_lfc^2) * sum(x$c21_lfc^2))) }))
  # PART 2: gene-call counts vs relabelings
  cnt <- dplyr::bind_rows(lapply(U, function(u) dplyr::mutate(u$counts, celltype=u$celltype))) %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise(T3_rank=emp_rank(T3, real), T3_relabel_median=median(T3[!real]),
                     A3_rank=emp_rank(A3, real), A3_relabel_median=median(A3[!real]),
                     T1_rank=emp_rank(T1, real), A1_rank=emp_rank(A1, real),
                     cos_rank=emp_rank(cos, real, FALSE), cos_relabel_median=round(median(cos[!real]), 3),
                     n_genes=n_genes[real][1], T3=T3[real][1], A3=A3[real][1], T1=T1[real][1], A1=A1[real][1],
                     n_labelings=dplyr::n(), .groups="drop") %>%
    dplyr::relocate(n_genes, T3, T3_rank, T3_relabel_median, A3, A3_rank, A3_relabel_median, T1, T1_rank, A1, A1_rank,
                    .after=celltype) %>%
    dplyr::left_join(cosr %>% dplyr::select(celltype, cos_mean_lfc), by="celltype") %>%
    dplyr::arrange(factor(celltype, levels=CT_ORDER))
  wcsv(cnt, odir, "Part1_2_Cosine_GeneCounts_vsRelabel.csv")
  cat("\n-- PART 1-2: DREADD own aging axis. Cosine (mean held-out LFCs) and buffered (T3) / amplified (A3) gene counts;",
      "rank = rank of the real labeling among its 20 (1 = most) --\n")
  print(as.data.frame(cnt %>% dplyr::mutate(cos_mean_lfc=round(cos_mean_lfc, 3))), row.names=FALSE)
  cat("   cos_rank = rank of the real cosine among its labelings (1 = most negative = most buffering-like)\n")
  calls_all <- dplyr::bind_rows(lapply(U, function(u) if (!is.null(u$calls) && nrow(u$calls)) dplyr::mutate(u$calls, celltype=u$celltype)))
  wcsv(calls_all, odir, "Part2_GeneCalls_RealLabeling.csv")

  # PART 1b: do the two cohorts age the same way? (DREADD PBS aging vs NE HC aging, old female)
  ax <- calls_all %>% dplyr::inner_join(ne_calls %>% dplyr::select(celltype, gene, ne_aging_lfc, ne_agechg),
                                        by=c("celltype","gene")) %>%
    dplyr::filter(!is.na(ne_aging_lfc)) %>% dplyr::group_by(celltype) %>%
    dplyr::summarise(n_shared=dplyr::n(),
                     cos_all=sum(aging_lfc * ne_aging_lfc) / sqrt(sum(aging_lfc^2) * sum(ne_aging_lfc^2)),
                     n_NE_aging=sum(ne_agechg),
                     pct_same_sign_NE_aging=100 * mean(sign(aging_lfc[ne_agechg]) == sign(ne_aging_lfc[ne_agechg])),
                     .groups="drop") %>%
    dplyr::mutate(dplyr::across(c(cos_all, pct_same_sign_NE_aging), ~round(.x, 2))) %>%
    dplyr::arrange(factor(celltype, levels=CT_ORDER))
  wcsv(ax, odir, "Part1b_AgingAxis_DREADDcohort_vs_NEcohort.csv")
  cat("\n-- PART 1b: do the two cohorts age the same way? (aging LFCs, DREADD PBS vs NE HC; old female) --\n")
  cat("   cos_all = cosine over shared genes; pct_same_sign_NE_aging = % of NE aging-changed genes changing the same way in DREADD (50 = chance)\n")
  print(as.data.frame(ax), row.names=FALSE)

  # PART 3: overlap with NE
  ov <- dplyr::bind_rows(lapply(U, function(u) if (!is.null(u$ov) && nrow(u$ov)) dplyr::mutate(u$ov, celltype=u$celltype)))
  if (nrow(ov)) {
    ov_ct <- ov %>% dplyr::group_by(celltype, test) %>%
      dplyr::summarise(n_ne_aging=n_ne_aging[real][1], n_ne_set=(a + b)[real][1],
                       pct_in_set=round(100 * a[real][1] / max(1, (a + b)[real][1]), 1),
                       pct_other=round(100 * c[real][1] / max(1, (c + d)[real][1]), 1),
                       odds_ratio=round(exp(log_or[real][1]), 2),
                       fisher_p=tryCatch(stats::fisher.test(matrix(c(a[real][1], c[real][1], b[real][1], d[real][1]), 2),
                                                            alternative="greater")$p.value, error=function(e) NA_real_),
                       rank=emp_rank(log_or, real), n_labelings=dplyr::n(), .groups="drop") %>%
      dplyr::arrange(test, factor(celltype, levels=CT_ORDER))
    wcsv(ov_ct, odir, "Part3_Overlap_withNE_PerCellType.csv")
    ov_pool <- ov %>% dplyr::filter(celltype %in% full_cts) %>% dplyr::group_by(test, key, real) %>% dplyr::summarise(sum_log_or=sum(log_or), n_ct=dplyr::n(), .groups="drop") %>%
      dplyr::group_by(test) %>% dplyr::summarise(rank=emp_rank(sum_log_or, real), n_labelings=dplyr::n(),
                                                 sum_log_or=round(sum_log_or[real][1], 2), n_ct=n_ct[real][1], .groups="drop")
    wcsv(ov_pool, odir, "Part3_Overlap_withNE_Pooled.csv")
    cat("\n-- PART 3: does DREADD move NE's genes the way NE does? (old female NE reference; NE aging direction) --\n")
    cat("   pct_in_set = % of NE-buffered (or NE-amplified) genes that DREADD reverses (or pushes further), in >= 2 of 3 fits\n")
    cat("   pct_other  = same % among the other NE aging-changed genes; rank = real labeling among 20 (1 = strongest)\n")
    print(as.data.frame(ov_ct %>% dplyr::mutate(fisher_p=signif(fisher_p, 2))), row.names=FALSE)
    cat("\n   Pooled across cell types (sum of log odds ratios; rank among 20 labelings):\n")
    print(as.data.frame(ov_pool), row.names=FALSE)
  }

  # PART 4: TF activity
  tf <- dplyr::bind_rows(lapply(U, function(u) if (!is.null(u$tf) && nrow(u$tf)) dplyr::mutate(u$tf, celltype=u$celltype)))
  if (nrow(tf)) {
    tf_ct <- tf %>% dplyr::group_by(celltype, source) %>%
      dplyr::summarise(rank_buf_own=emp_rank(S_own, real), rank_amp_own=emp_rank(S_own, real, FALSE),
                       rank_buf_NEax=emp_rank(S_NEax, real), rank_amp_NEax=emp_rank(S_NEax, real, FALSE),
                       n_labelings=dplyr::n(),
                       S_own=S_own[real][1], S_NEax=S_NEax[real][1], c21_act=c21_act[real][1], .groups="drop") %>%
      dplyr::filter(!is.na(S_own)) %>%
      dplyr::mutate(DR_call_NEax=dplyr::case_when(S_NEax > 0 & rank_buf_NEax <= rank_cut(n_labelings) ~ "buffered",
                                                  S_NEax < 0 & rank_amp_NEax <= rank_cut(n_labelings) ~ "amplified", TRUE ~ "-"),
                    DR_call_own=dplyr::case_when(S_own > 0 & rank_buf_own <= rank_cut(n_labelings) ~ "buffered",
                                                 S_own < 0 & rank_amp_own <= rank_cut(n_labelings) ~ "amplified", TRUE ~ "-")) %>%
      dplyr::left_join(ne_act, by=c("celltype","source"))
    wcsv(tf_ct, odir, "Part4_TFActivity_PerCellType_withNE.csv")
    tf_all <- tf %>% dplyr::filter(celltype %in% full_cts) %>% dplyr::group_by(key, real, source) %>%
      dplyr::summarise(S_NEax_sum=sum(S_NEax, na.rm=TRUE), S_own_sum=sum(S_own), n_ct=dplyr::n(), .groups="drop") %>%
      dplyr::filter(n_ct >= 2) %>% dplyr::group_by(source) %>%
      dplyr::summarise(rank_buf_NEax=emp_rank(S_NEax_sum, real), rank_amp_NEax=emp_rank(S_NEax_sum, real, FALSE),
                       rank_buf_own=emp_rank(S_own_sum, real), rank_amp_own=emp_rank(S_own_sum, real, FALSE),
                       S_NEax_sum=S_NEax_sum[real][1], S_own_sum=S_own_sum[real][1], n_ct=n_ct[real][1], .groups="drop") %>%
      dplyr::relocate(S_NEax_sum, .after=source) %>% dplyr::relocate(S_own_sum, .after=rank_amp_NEax)
    wcsv(tf_all, odir, "Part4_TFActivity_AcrossCellTypes.csv")

    # agreement: NE call vs DREADD call (NE aging axis), per cell type
    agree <- tf_ct %>% dplyr::filter(NE_call != "-") %>%
      dplyr::group_by(celltype) %>%
      dplyr::summarise(n_NE_called=dplyr::n(), same_direction=sum(sign(S_NEax) == sign(NE_S), na.rm=TRUE),
                       DR_called_same=sum(DR_call_NEax == NE_call), .groups="drop") %>%
      dplyr::mutate(pct_same_direction=round(100 * same_direction / n_NE_called, 1))
    # calibrated: per labeling, % of NE-called TFs whose DREADD score has NE's sign; rank of the real labeling
    agree_rl <- tf %>% dplyr::inner_join(ne_act %>% dplyr::filter(NE_call != "-") %>% dplyr::select(celltype, source, NE_S),
                                         by=c("celltype","source")) %>%
      dplyr::group_by(celltype, key, real) %>%
      dplyr::summarise(pct=100 * mean(sign(S_NEax) == sign(NE_S), na.rm=TRUE), .groups="drop") %>%
      dplyr::group_by(celltype) %>%
      dplyr::summarise(rank=emp_rank(pct, real), relabel_median_pct=round(median(pct[!real]), 1), .groups="drop")
    agree <- agree %>% dplyr::left_join(agree_rl, by="celltype")
    wcsv(agree, odir, "Part4_Agreement_NEcalledTFs.csv")
    cat("\n-- PART 4a: TFs NE calls buffered/amplified (old female): does DREADD move them the same way? --\n")
    cat("   same_direction = DREADD score has the same sign as NE's; relabel_median_pct = same % in relabeled mice; rank among 20\n")
    print(as.data.frame(agree %>% dplyr::arrange(factor(celltype, levels=CT_ORDER))), row.names=FALSE)

    cat("\n-- PART 4b: shortlist TFs, per cell type where NE called them (NE vs DREADD, NE aging axis) --\n")
    print(as.data.frame(tf_ct %>% dplyr::filter(source %in% SHORTLIST, NE_call != "-") %>%
      dplyr::mutate(dplyr::across(c(NE_S, S_NEax, S_own), ~round(.x, 2))) %>%
      dplyr::arrange(source, factor(celltype, levels=CT_ORDER)) %>%
      dplyr::select(source, celltype, NE_call, NE_S, DR_call_NEax, S_NEax, rank_buf_NEax, rank_amp_NEax, DR_call_own, S_own)),
      row.names=FALSE)
    cat("\n-- PART 4c: shortlist TFs across cell types (DREADD; rank among 20 labelings) --\n")
    print(as.data.frame(tf_all %>% dplyr::filter(source %in% SHORTLIST) %>%
      dplyr::mutate(dplyr::across(c(S_NEax_sum, S_own_sum), ~round(.x, 1))) %>% dplyr::arrange(source)), row.names=FALSE)

    # figure: NE S vs DREADD S (NE axis), shortlist labelled, per cell type with NE reference
    d4 <- tf_ct %>% dplyr::filter(!is.na(NE_S)) %>%
      dplyr::mutate(lab=ifelse(source %in% SHORTLIST & NE_call != "-", source, NA),
                    celltype=factor(celltype, levels=CT_ORDER))
    if (nrow(d4)) sv(ggplot(d4, aes(NE_S, S_NEax)) +
        geom_hline(yintercept=0, color="grey80") + geom_vline(xintercept=0, color="grey80") +
        geom_point(aes(color=NE_call), size=0.6, alpha=0.6) +
        geom_text(aes(label=lab), size=2, vjust=-0.7, na.rm=TRUE) +
        scale_color_manual(values=c(buffered="#1B5E20", amplified="#6A1B9A", `-`="grey75"), name="NE call") +
        facet_wrap(~celltype, scales="free") + th +
        labs(title=sprintf("TF buffering score: NE (x) vs DREADD (y), old female, Approach %s", ap),
             subtitle="Both against NE's aging direction; S > 0 = reverses aging, S < 0 = pushes further",
             x="NE buffering score S", y="DREADD buffering score (NE aging axis)"),
      odir, "F_TF_NE_vs_DREADD", 14, 10)
  }
  cat("\nOutputs:", odir, "\n")
}
banner("Done")
