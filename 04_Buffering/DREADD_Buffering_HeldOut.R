
################################################################################
### DREADD_Buffering_HeldOut.R
###
### Held-out cosine for the DREADD dataset (female only, young + old),
### Approach B and Approach D. Same held-out logic as the NE robustness step
### (Buffering_Pipeline.R, step 1).
###
###   B (aging)      = OPBS vs YPBS
###   E (activation) = OC21 vs OPBS
###
### Coupled  : all OPBS mice in both B and E
### Held-out : for each OPBS mouse k
###              B_k = OPBS (without k) vs YPBS
###              E_k = OC21 vs OPBS mouse k
###            -> no mouse on both sides; cosines averaged over k
### Gene-level permutation on the averaged held-out cosine + BH FDR (descriptive)
### Mouse-noise null : N_j = OC21 (without j) vs OC21 mouse j (C21 vs C21),
###            cos(B_k, N_j) averaged over k, both signs -> SD, ratio, beyond-all
###
### Cell selection and DESeq2 designs as in DREADD_DESeq2.R (03_DREADD):
###   Approach B : neurons mScarlet_raw >= thr, glia all cells
###                neurons ~ mScarlet_median_scaled + grp, glia ~ grp
###   Approach D : neurons mScarlet_raw >= thr, range trim (10-90th pct,
###                median bounds across mice), binding-constraint downsample
###                (set.seed(42)); all ~ grp
### All genes, unshrunken LFC, gene filter rowSums(counts >= 10) >= 2.
###
### Uses df_harmony if loaded in the session, otherwise reads DREADD_OBJECT.
### Writes the pseudobulk caches pseudobulk_mScarlet_50_Approach_{B,D}.RDS (reused by Buffering_DREADD_vs_NE.R,
### Buffering_Paper_Figures.R and Buffering_Pipeline.R).
################################################################################

suppressPackageStartupMessages({
  library(conflicted)
  library(Seurat)
  library(BPCells)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  conflicts_prefer(dplyr::select,       .quiet=TRUE)
  conflicts_prefer(dplyr::filter,       .quiet=TRUE)
  conflicts_prefer(dplyr::mutate,       .quiet=TRUE)
  conflicts_prefer(base::intersect,     .quiet=TRUE)
  conflicts_prefer(base::setdiff,       .quiet=TRUE)
  conflicts_prefer(base::as.data.frame, .quiet=TRUE)
})

################################################################################
### PARAMETERS
################################################################################

BASE_PATH     <- "PATH_TO_DREADD_OUTPUT"
DREADD_OBJECT <- "PATH_TO_DREADD_OBJECT.rds"   # integrated DREADD Seurat object
fig_date  <- format(Sys.Date(),"%Y-%m-%d")
OUT_DIR   <- file.path(BASE_PATH, "DREADD_Buffering_HeldOut", fig_date)
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
cat("Output:", OUT_DIR, "\n\n")

APPROACHES       <- c("B","D")
MSCARLET_STATES  <- c(50)        # add 0 for the non-thresholded state: c(0, 50)
CELLTYPES_TO_RUN <- NULL         # e.g. c("DG ExNs") for a quick test
N_PERM           <- 2000
MIN_GENES        <- 20
MIN_GENES_REPORT <- 1000
TRIM_LO <- 0.10; TRIM_HI <- 0.90; MIN_CELLS_FLOOR <- 5

IN_SUBTYPES        <- c("IN1","IN2")
NEURONAL_CELLTYPES <- c("CA1 ExNs","CA2 ExNs","CA3 ExNs","DG ExNs","Sub ExNs","IN")
CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN",
              "Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal",
              "Cajal-Retzius Cell","Macrophage","Epithelial","T Cell")
CT_COLS <- c("DG ExNs"="#E53935","CA1 ExNs"="#CC0000","CA2 ExNs"="#8E24AA",
             "CA3 ExNs"="#5E35B1","Sub ExNs"="#80D8FF","IN"="#00897B",
             "Astrocytes"="#43A047","Oligodendrocytes"="#1E88E5","OPC"="#00ACC1",
             "Microglia"="#7CB342","Endothelial"="#FB8C00","Pericyte"="#6D4C41",
             "APC"="#546E7A","Stromal"="#27408B","Macrophage"="#C62828",
             "Cajal-Retzius Cell"="#FDD835","T Cell"="#757575")

if (!exists("df_harmony")) df_harmony <- readRDS(DREADD_OBJECT)

################################################################################
### STEP 1: METADATA
################################################################################

in_present <- intersect(IN_SUBTYPES, unique(df_harmony@meta.data$celltype))
if (length(in_present) > 0)
  df_harmony@meta.data$celltype[df_harmony@meta.data$celltype %in% in_present] <- "IN"
if (!"mScarlet_raw" %in% colnames(df_harmony@meta.data)) {
  raw <- FetchData(df_harmony, vars="mScarlet", layer="counts")
  df_harmony@meta.data$mScarlet_raw <- raw$mScarlet
}
df_harmony@meta.data <- df_harmony@meta.data %>%
  dplyr::mutate(age_norm = dplyr::case_when(tolower(age)=="young" ~ "young",
                                            tolower(age)=="old"   ~ "old",
                                            TRUE ~ tolower(age)),
                condition_norm = experimental.condition,
                age_condition  = paste0(age_norm,"_",experimental.condition))

meta_full <- df_harmony@meta.data
if ("cell_id" %in% colnames(meta_full)) meta_full <- meta_full %>% dplyr::select(-cell_id)
meta_full <- meta_full %>% tibble::rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, celltype, orig.ident, age_norm, condition_norm,
                age_condition, mScarlet_raw)

################################################################################
### STEP 2: CELL SELECTION PER APPROACH
################################################################################

select_cells <- function(approach, thr) {
  base <- meta_full %>%
    dplyr::filter((celltype %in% NEURONAL_CELLTYPES & mScarlet_raw >= thr) |
                  (!celltype %in% NEURONAL_CELLTYPES))
  if (approach == "B") return(base)

  # Approach D: range trim (neurons) ...
  trimmed <- lapply(unique(base$celltype), function(ct) {
    cm <- base %>% dplyr::filter(celltype == ct)
    if (!ct %in% NEURONAL_CELLTYPES) return(cm)
    bnd <- dplyr::bind_rows(lapply(unique(cm$orig.ident), function(m) {
      v <- cm$mScarlet_raw[cm$orig.ident == m]
      if (length(v) < 5) return(NULL)
      data.frame(lo=quantile(v,TRIM_LO), hi=quantile(v,TRIM_HI))
    }))
    if (nrow(bnd) == 0) return(cm)
    lo <- median(bnd$lo); hi <- median(bnd$hi)
    if (lo >= hi) return(cm)
    cm %>% dplyr::filter(mScarlet_raw >= lo & mScarlet_raw <= hi)
  }) %>% dplyr::bind_rows()

  # ... + binding-constraint downsample (same seed as DREADD_DESeq2.R)
  set.seed(42)
  tgt <- trimmed %>% dplyr::count(celltype, orig.ident, age_condition) %>%
    dplyr::group_by(celltype) %>% dplyr::summarise(target_n=min(n), .groups="drop") %>%
    dplyr::filter(target_n >= MIN_CELLS_FLOOR)
  lapply(tgt$celltype, function(ct) {
    t <- tgt$target_n[tgt$celltype == ct]
    trimmed %>% dplyr::filter(celltype == ct) %>%
      dplyr::group_by(orig.ident, age_condition) %>%
      dplyr::group_modify(~dplyr::slice_sample(.x, n=min(t, nrow(.x)))) %>%
      dplyr::ungroup()
  }) %>% dplyr::bind_rows()
}

parse_pb_colname <- function(s) {
  ct_name <- CT_ORDER[sapply(CT_ORDER, function(ct) startsWith(s, paste0(ct,"_")))]
  if (length(ct_name) == 0) return(list(celltype=NA, orig.ident=NA, condition=NA, age=NA))
  ct_name <- ct_name[1]
  parts <- strsplit(sub(paste0("^",ct_name,"_"), "", s), "_")[[1]]; n <- length(parts)
  list(celltype=ct_name,
       orig.ident=gsub("-","_",paste(parts[1:(n-2)],collapse="_")),
       condition=parts[n-1], age=parts[n])
}

make_pseudobulk <- function(cells) {
  so <- df_harmony[, cells$cell_id]
  idx <- match(colnames(so), cells$cell_id)
  so@meta.data$age_norm       <- cells$age_norm[idx]
  so@meta.data$condition_norm <- cells$condition_norm[idx]
  pb <- AggregateExpression(object=so, assays="RNA", return.seurat=FALSE,
          group.by=c("celltype","orig.ident","condition_norm","age_norm"),
          normalization.method=NULL, slot="counts", verbose=FALSE)$RNA
  parsed <- lapply(colnames(pb), parse_pb_colname)
  meta <- data.frame(sample=colnames(pb),
                     celltype  =sapply(parsed,"[[","celltype"),
                     orig.ident=sapply(parsed,"[[","orig.ident"),
                     condition =sapply(parsed,"[[","condition"),
                     age       =sapply(parsed,"[[","age"), stringsAsFactors=FALSE)
  keep <- !is.na(meta$celltype)
  list(counts=as.matrix(pb[, keep, drop=FALSE]), meta=meta[keep, ])
}

################################################################################
### STEP 3: FITS, COSINE, PERMUTATION
################################################################################

# LFC of num vs den (unshrunken). Approach B neurons: ~ mScarlet_median_scaled + grp
fit_lfc <- function(PB, num, den, ct, approach, mcov) {
  cols <- c(den, num)
  cnt  <- round(PB$counts[, cols, drop=FALSE])
  cnt  <- cnt[rowSums(cnt >= 10) >= 2, , drop=FALSE]
  if (nrow(cnt) < 10) return(NULL)
  cd <- data.frame(grp=factor(rep(c("den","num"), c(length(den), length(num))),
                              levels=c("den","num")), row.names=cols)
  design <- ~grp
  if (approach == "B" && ct %in% NEURONAL_CELLTYPES && !is.null(mcov)) {
    oid <- PB$meta$orig.ident[match(cols, PB$meta$sample)]
    mv  <- mcov$mScarlet_median[match(paste(oid, ct), paste(mcov$orig.ident, mcov$celltype))]
    if (!any(is.na(mv)) && var(mv) > 0) {
      cd$mScarlet_median_scaled <- as.numeric(scale(mv))
      design <- ~mScarlet_median_scaled + grp
    }
  }
  tryCatch({
    dds <- suppressMessages(DESeqDataSetFromMatrix(cnt, cd, design=design))
    dds <- tryCatch(DESeq(dds, quiet=TRUE), error=function(e) {
      d <- estimateSizeFactors(dds); d <- estimateDispersionsGeneEst(d)
      dispersions(d) <- mcols(d)$dispGeneEst; nbinomWaldTest(d, quiet=TRUE)
    })
    r <- results(dds, name="grp_num_vs_den"); setNames(r$log2FoldChange, rownames(r))
  }, error=function(e) NULL)
}

align <- function(B, E) {
  if (is.null(B) || is.null(E)) return(NULL)
  g <- intersect(names(B)[!is.na(B)], names(E)[!is.na(E)])
  if (length(g) < MIN_GENES) return(NULL)
  list(b=B[g], e=E[g])
}
cosv <- function(B, E) { x <- align(B, E); if (is.null(x)) NA_real_ else
  sum(x$b*x$e)/(sqrt(sum(x$b^2))*sqrt(sum(x$e^2))) }

# Gene-level permutation on the MEAN of the held-out cosines
perm_mean <- function(pairs) {
  pairs <- Filter(Negate(is.null), pairs)
  if (length(pairs) == 0 || N_PERM == 0) return(c(p=NA, sd=NA))
  pr  <- lapply(pairs, function(x) list(b=x$b, e=x$e, nb=sqrt(sum(x$b^2)), ne=sqrt(sum(x$e^2))))
  obs <- mean(sapply(pr, function(x) sum(x$b*x$e)/(x$nb*x$ne)))
  nul <- replicate(N_PERM, mean(sapply(pr, function(x) sum(x$b*sample(x$e))/(x$nb*x$ne))))
  c(p=(1+sum(nul <= obs))/(N_PERM+1), sd=sd(nul))
}

################################################################################
### STEP 4: MAIN LOOP
################################################################################

rows <- list(); perfit <- list(); nullv <- list()
for (thr in MSCARLET_STATES) for (ap in APPROACHES) {
  state <- paste0("mScarlet_", thr, "_Approach_", ap)
  cat("\n================", state, "================\n")
  cache <- file.path(BASE_PATH, "DREADD_Buffering_HeldOut",
                     paste0("pseudobulk_", state, ".RDS"))
  if (file.exists(cache)) { PB <- readRDS(cache); cells <- NULL } else {
    cells <- select_cells(ap, thr)
    PB <- make_pseudobulk(cells)
    PB$mcov <- cells %>% dplyr::filter(celltype %in% NEURONAL_CELLTYPES) %>%
      dplyr::group_by(orig.ident, celltype) %>%
      dplyr::summarise(mScarlet_median=median(mScarlet_raw), .groups="drop")
    saveRDS(PB, cache)
  }
  print(table(paste(PB$meta$age, PB$meta$condition)) / length(unique(PB$meta$celltype)))

  cts <- intersect(CT_ORDER, unique(PB$meta$celltype))
  if (!is.null(CELLTYPES_TO_RUN)) cts <- intersect(CELLTYPES_TO_RUN, cts)

  for (ct in cts) {
    m    <- PB$meta[PB$meta$celltype == ct, ]
    ypbs <- m$sample[m$age=="young" & m$condition=="PBS"]
    opbs <- m$sample[m$age=="old"   & m$condition=="PBS"]
    oc21 <- m$sample[m$age=="old"   & m$condition=="C21"]
    cat(sprintf("  %-20s YPBS=%d OPBS=%d OC21=%d", ct, length(ypbs), length(opbs), length(oc21)))
    if (length(ypbs) < 2 || length(opbs) < 2 || length(oc21) < 2) { cat("  skip\n"); next }
    f <- function(num, den) fit_lfc(PB, num, den, ct, ap, PB$mcov)

    # Coupled (current cosine)
    cos_coupled <- cosv(f(opbs, ypbs), f(oc21, opbs))

    # Held-out
    Bs <- lapply(seq_along(opbs), function(k) f(opbs[-k], ypbs))
    Es <- lapply(seq_along(opbs), function(k) f(oc21, opbs[k]))
    pairs <- mapply(align, Bs, Es, SIMPLIFY=FALSE)
    cos_k <- mapply(cosv, Bs, Es)
    obs   <- mean(cos_k, na.rm=TRUE)
    pm    <- perm_mean(pairs)
    for (k in seq_along(opbs))
      perfit[[length(perfit)+1]] <- data.frame(state=state, celltype=ct,
        held_out=m$orig.ident[match(opbs[k], m$sample)], cosine=cos_k[k])

    # Mouse-noise null: C21 vs C21 (needs >= 3 OC21 mice)
    null_j <- NA; sd_m <- NA; maxabs <- NA
    if (length(oc21) >= 3) {
      Ns <- lapply(seq_along(oc21), function(j) f(oc21[-j], oc21[j]))
      null_j <- sapply(Ns, function(N) mean(sapply(Bs, cosv, E=N), na.rm=TRUE))
      sd_m   <- sd(c(null_j, -null_j), na.rm=TRUE)
      maxabs <- max(abs(null_j), na.rm=TRUE)
      nullv[[length(nullv)+1]] <- data.frame(state=state, celltype=ct, null_cos=null_j,
        c21_ref=m$orig.ident[match(oc21, m$sample)])
    }
    n_g <- if (!is.null(pairs[[1]])) length(pairs[[1]]$b) else NA
    rows[[length(rows)+1]] <- data.frame(
      state=state, approach=ap, mscarlet_thr=thr, celltype=ct, n_genes=n_g,
      n_YPBS=length(ypbs), n_OPBS=length(opbs), n_OC21=length(oc21),
      cosine_coupled=cos_coupled, cosine_heldout=obs,
      cosine_held_min=min(cos_k,na.rm=TRUE), cosine_held_max=max(cos_k,na.rm=TRUE),
      p_perm_heldout=pm[["p"]], gene_perm_sd=pm[["sd"]],
      mouse_null_sd=sd_m, mouse_null_maxabs=maxabs,
      z_mouse=obs/sd_m, beyond_all_mouse_null=obs < -maxabs,
      stringsAsFactors=FALSE)
    cat(sprintf(" | coupled %.3f | held-out %.3f | mouse SD %.3f\n", cos_coupled, obs, sd_m))
  }
}

summ <- dplyr::bind_rows(rows) %>%
  dplyr::mutate(reportable = n_genes >= MIN_GENES_REPORT) %>%
  dplyr::group_by(state, reportable) %>%
  dplyr::mutate(fdr_perm_heldout = ifelse(reportable, p.adjust(p_perm_heldout,"BH"), NA)) %>%
  dplyr::ungroup()
write.csv(summ, file.path(OUT_DIR, paste0(fig_date,"_DREADD_HeldOut_Summary.csv")), row.names=FALSE)
write.csv(dplyr::bind_rows(perfit), file.path(OUT_DIR, paste0(fig_date,"_DREADD_HeldOut_PerFit.csv")), row.names=FALSE)
if (length(nullv) > 0)
  write.csv(dplyr::bind_rows(nullv), file.path(OUT_DIR, paste0(fig_date,"_DREADD_MouseNoise_Values.csv")), row.names=FALSE)

cat("\n=== SUMMARY ===\n")
print(as.data.frame(summ %>% dplyr::filter(reportable) %>%
  dplyr::select(state, celltype, n_OPBS, n_OC21, cosine_coupled, cosine_heldout,
                fdr_perm_heldout, mouse_null_sd, z_mouse, beyond_all_mouse_null) %>%
  dplyr::mutate(dplyr::across(where(is.double), ~signif(.x,3)))))

################################################################################
### STEP 5: FIGURE
################################################################################

sm <- summ %>% dplyr::filter(reportable)
lv <- rev(intersect(CT_ORDER, unique(sm$celltype)))
pf <- dplyr::bind_rows(perfit) %>% dplyr::semi_join(sm, by=c("state","celltype"))
nv <- if (length(nullv) > 0) dplyr::bind_rows(nullv) %>% dplyr::semi_join(sm, by=c("state","celltype")) else NULL
fx <- function(d) d %>% dplyr::mutate(celltype=factor(celltype, levels=lv))
p <- ggplot() + geom_vline(xintercept=0, linetype="dashed", color="grey30")
if (!is.null(nv)) p <- p +
  geom_point(data=fx(dplyr::bind_rows(nv, dplyr::mutate(nv, null_cos=-null_cos))),
             aes(x=null_cos, y=celltype), color="grey75", size=1.4)
p <- p +
  geom_point(data=fx(pf), aes(x=cosine, y=celltype), color="grey35", size=1.4) +
  geom_point(data=fx(sm), aes(x=cosine_heldout, y=celltype, fill=celltype),
             shape=21, size=3.4, color="black", stroke=0.4) +
  geom_point(data=fx(sm), aes(x=cosine_coupled, y=celltype), shape=5, size=2.4, color="grey20") +
  scale_fill_manual(values=CT_COLS, guide="none") +
  facet_wrap(~state) + theme_classic(base_size=11) +
  labs(title="DREADD: held-out cos(OPBSYPBS, OC21OPBS), all genes, female",
       subtitle=paste0("Open diamond = coupled | filled = held-out mean | ",
                       "dark grey = each held-out OPBS mouse | light grey = C21-vs-C21 mouse noise"),
       x="Cosine similarity", y=NULL)
stem <- file.path(OUT_DIR, paste0(fig_date,"_DREADD_HeldOut_Cosine"))
h <- max(4.5, length(lv)*0.42+3)
ggsave(paste0(stem,".pdf"), p, width=5*length(unique(sm$state))+2, height=h, limitsize=FALSE)
ggsave(paste0(stem,".png"), p, width=5*length(unique(sm$state))+2, height=h, dpi=600, limitsize=FALSE)
cat("\n=== DONE ===\nOutput:", OUT_DIR, "\n")
