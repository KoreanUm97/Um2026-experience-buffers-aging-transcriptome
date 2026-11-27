################################################################################
### Fig4f_GO_ARG_Classes.R
###
### Fig. 4f: GO biological-process enrichment of DREADD-induced ARG classes, all cell types pooled (female only).
###   ARGs per cell type and age = C21 vs PBS padj < 0.05 (no fold-change cutoff), induced (up) and repressed (down)
###   separately, from the DESeq2 tables of DREADD_DESeq2.R (_MergedIN, latest date):
###     Young ARGs   ARG in young            Old ARGs   ARG in old
###     Maintained   ARG at both ages        Lost       ARG in young, not in old       Gained   ARG in old, not in young
###   (sets anchored on one age, as in the eNE pathway analysis: Maintained / Lost need the gene tested in young;
###    Gained needs it tested in old.)
###   Superpooled = union of each set over all cell types.
###   GO (clusterProfiler enrichGO, BP): BH-adjusted p < 0.05 and q < 0.2; sets >= 5 genes; terms of 10-500 genes;
###   simplify(cutoff 0.7, by p.adjust); background = all genes in the anchor comparison's DESeq2 table (YC21YPBS for
###   Young / Maintained / Lost, OC21OPBS for Old / Gained), pooled across cell types. Same settings as the eNE pathways.
###   Plot (Fig. 4f): Maintained Up / Down and Lost Up / Down, top 5 terms each; colour = -log10(padj), size = gene ratio.
###
### Outputs (<OUT_DIR>/Approach<X>/): ARG_SetSizes.csv, ARG_SetMembers.csv, GO_BP_Superpooled_AllResults.csv,
###   GO_BP_Superpooled_Significant.csv, SuperPooled_MaintainedLost_UpDown_GO_BP.pdf/.png/.csv
### APPROACHES: B = covariate, D = trimming + downsampling (both are drawn; pick the one stated in the legend).
### Run:  Rscript Fig4f_GO_ARG_Classes.R     (set DESEQ_ROOT below first; N_CORES > 1 uses forked parallel workers)
################################################################################

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(conflicted); library(dplyr); library(tidyr); library(ggplot2); library(parallel)
  conflicts_prefer(dplyr::select, .quiet = TRUE); conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicts_prefer(dplyr::mutate, .quiet = TRUE); conflicts_prefer(base::intersect, .quiet = TRUE)
  conflicts_prefer(base::setdiff, .quiet = TRUE); conflicts_prefer(base::union, .quiet = TRUE)
  conflicts_prefer(base::unique, .quiet = TRUE); conflicts_prefer(base::Filter, .quiet = TRUE)
})

################################################################################
### SETTINGS
################################################################################

DESEQ_ROOT <- "PATH_TO_DREADD_OUTPUT"                    # same as OUT_ROOT in DREADD_DESeq2.R
OUT_DIR    <- file.path(DESEQ_ROOT, "Paper_Figures", "Figure4", "GO")
APPROACHES <- c("B", "D")
DEG_DIRS   <- c(B = "DESeq2_mScarletPos_ApproachB_MergedIN", D = "DESeq2_mScarletPos_ApproachD_MergedIN")
APPROACH_LABS <- c(B = "Covariate", D = "Downsampling")
ARG_PADJ   <- 0.05
DIRECTIONS <- c("up", "down")
N_CORES    <- 8
CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN","Astrocytes","Oligodendrocytes","OPC",
              "Microglia","Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell")
SETS       <- c("Young ARGs","Old ARGs","Maintained","Lost","Gained")
SET_ANCHOR <- list(`Young ARGs` = "young", `Old ARGs` = "old", Maintained = "young", Lost = "young", Gained = "old")

# ---- GO settings (identical to the eNE pathway analysis) ----------------------------------------------------
GO <- list(
  ont        = "BP",
  fdr        = 0.05,      # enrichGO pvalueCutoff (applied to BH-adjusted p)
  q          = 0.2,       # enrichGO qvalueCutoff
  min_set    = 5,         # minimum genes per gene set
  min_count  = 1,         # no minimum genes per term
  min_gs     = 10,        # enrichGO defaults
  max_gs     = 500,
  simplify   = TRUE,      # clusterProfiler::simplify(cutoff = 0.7, by = "p.adjust")
  sim_cutoff = 0.7,
  universe   = "comparison"   # background = all genes in the anchor comparison's DESeq2 table, pooled across cell types
)

.need <- function() {
  if (!requireNamespace("clusterProfiler", quietly=TRUE) || !requireNamespace("org.Mm.eg.db", quietly=TRUE))
    stop("Install clusterProfiler and org.Mm.eg.db (BiocManager::install(c('clusterProfiler','org.Mm.eg.db','GOSemSim')))")
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Mm.eg.db) })
}

go_enrich <- function(genes, universe, ont=GO$ont) {
  .need()
  universe <- unique(universe); genes <- intersect(unique(genes), universe)
  if (length(genes) < GO$min_set) return(NULL)
  eg <- tryCatch(clusterProfiler::enrichGO(genes, OrgDb=org.Mm.eg.db::org.Mm.eg.db, keyType="SYMBOL", ont=ont,
                                           universe=universe, pAdjustMethod="BH", pvalueCutoff=GO$fdr, qvalueCutoff=GO$q,
                                           minGSSize=GO$min_gs, maxGSSize=GO$max_gs, readable=TRUE), error=function(e) NULL)
  if (is.null(eg) || !nrow(as.data.frame(eg))) return(NULL)                # only significant terms are returned
  if (GO$min_count > 1) { eg@result <- eg@result[eg@result$Count >= GO$min_count, , drop=FALSE]; if (!nrow(eg@result)) return(NULL) }
  if (GO$simplify && requireNamespace("GOSemSim", quietly=TRUE))
    eg <- tryCatch(clusterProfiler::simplify(eg, cutoff=GO$sim_cutoff, by="p.adjust",
                                             semData=if (exists("GO_SEMDATA")) GO_SEMDATA else NULL), error=function(e) eg)
  .tidy(as.data.frame(eg), length(genes), length(universe), paste0("GO_", ont))
}

.tidy <- function(r, n_set, n_universe, db) {
  frac <- function(x) sapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
  r %>% dplyr::mutate(database=db, n_set=n_set, n_universe=n_universe, fold=frac(GeneRatio) / frac(BgRatio),
                      nlp=-log10(p.adjust), sig=p.adjust < GO$fdr, .before=1)
}

go_batch <- function(jobs, n_cores=1) {
  .need()
  # Build the GO annotation and semantic-similarity data ONCE here, so every parallel worker inherits them
  # (otherwise each worker rebuilds them, which is what makes the run slow)
  message("Preparing GO data once for all workers...")
  invisible(tryCatch(clusterProfiler::enrichGO(c("Fos","Arc","Npas4","Egr1","Junb","Nr4a1"), OrgDb=org.Mm.eg.db::org.Mm.eg.db,
                     keyType="SYMBOL", ont=GO$ont, pvalueCutoff=1, qvalueCutoff=1), error=function(e) NULL))
  if (GO$simplify && requireNamespace("GOSemSim", quietly=TRUE))
    GO_SEMDATA <<- GOSemSim::godata(org.Mm.eg.db::org.Mm.eg.db, ont=GO$ont, computeIC=FALSE)   # Wang measure needs no IC
  one <- function(j) {
    r <- go_enrich(j$genes, j$universe)
    if (is.null(r) || !nrow(r)) return(NULL)
    cbind(j$meta[rep(1, nrow(r)), , drop=FALSE], r, row.names=NULL)
  }
  safe <- function(i) tryCatch(one(jobs[[i]]), error=function(e) structure(conditionMessage(e), class="go_error"))
  res <- if (n_cores > 1) parallel::mclapply(seq_along(jobs), safe, mc.cores=n_cores, mc.preschedule=FALSE) else lapply(seq_along(jobs), safe)
  lost <- which(sapply(res, inherits, "try-error"))
  if (length(lost)) { message("Rerunning ", length(lost), " enrichments one at a time"); res[lost] <- lapply(lost, safe) }
  err <- Filter(function(x) inherits(x, "go_error"), res)
  if (length(err)) { message(length(err), " enrichments failed; first: ", unclass(err[[1]])) }
  dplyr::bind_rows(Filter(is.data.frame, res))
}

################################################################################
### ARG CALLS AND SETS
################################################################################

read_deg <- function(ap, comp) {
  f <- list.files(file.path(DESEQ_ROOT, DEG_DIRS[[ap]]), pattern = paste0("_", comp, "_ALL\\.csv$"), recursive = TRUE, full.names = TRUE)
  if (!length(f)) stop("No ", comp, "_ALL.csv for Approach ", ap)
  f <- sort(f, decreasing = TRUE)[1]; cat(sprintf("  Approach %s %s: %s\n", ap, comp, f))
  read.csv(f, stringsAsFactors = FALSE) %>% dplyr::transmute(celltype, gene, lfc = log2FoldChange, pvalue, padj)
}
is_arg <- function(lfc, q, dir) !is.na(q) & !is.na(lfc) & q < ARG_PADJ & (if (dir == "up") lfc > 0 else lfc < 0)

# compareCluster-style dot plot (as in the eNE pathway figures): x = gene set, top_n terms per set
compare_dot <- function(df_all, title, out_stem, top_n = 5, width = 10, height = 11) {
  if (!nrow(df_all) || dplyr::n_distinct(df_all$Cluster) < 2) return(invisible(NULL))
  write.csv(df_all, paste0(out_stem, ".csv"), row.names = FALSE)
  df <- df_all %>% dplyr::group_by(Cluster) %>% dplyr::arrange(p.adjust, .by_group = TRUE) %>% dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(GeneRatio_num = sapply(GeneRatio, function(x) { v <- strsplit(x, "/")[[1]]; as.numeric(v[1]) / as.numeric(v[2]) }),
                  neg_log10_padj = -log10(p.adjust))
  ord <- df %>% dplyr::group_by(Description) %>% dplyr::summarise(min_p = min(p.adjust), .groups = "drop") %>%
    dplyr::arrange(min_p) %>% dplyr::pull(Description)
  df <- df %>% dplyr::mutate(Description = factor(Description, levels = rev(ord)), Cluster = factor(Cluster, levels = unique(df_all$Cluster)))
  p <- ggplot(df, aes(x = Cluster, y = Description, color = neg_log10_padj, size = GeneRatio_num)) + geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#3182BD", high = "#E6550D", name = bquote(-log[10](italic(p)[adj]))) +
    scale_size_continuous(range = c(1.5, 7), name = "Gene ratio") + theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7.5), axis.text.y = element_text(size = 7.5),
          axis.title = element_blank(), legend.position = "right", legend.text = element_text(size = 7),
          legend.title = element_text(size = 7.5), plot.margin = margin(5, 10, 5, 5),
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 8.5, color = "grey40", hjust = 0.5)) +
    labs(title = title, subtitle = paste0("Top ", top_n, " terms per gene set | GO ", GO$ont))
  ggsave(paste0(out_stem, ".pdf"), p, width = width, height = height, device = "pdf", limitsize = FALSE)
  ggsave(paste0(out_stem, ".png"), p, width = width, height = height, device = "png", dpi = 600, bg = "white", limitsize = FALSE)
  cat("  saved:", basename(out_stem), "\n")
}

for (ap in APPROACHES) {
  cat("\n=== Approach", ap, "===\n")
  OUT <- file.path(OUT_DIR, paste0("Approach", ap)); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  dy <- read_deg(ap, "YC21YPBS"); do <- read_deg(ap, "OC21OPBS")
  calls <- dplyr::full_join(dy %>% dplyr::rename(lfc_young = lfc, pvalue_young = pvalue, padj_young = padj),
                            do %>% dplyr::rename(lfc_old = lfc, pvalue_old = pvalue, padj_old = padj), by = c("celltype", "gene")) %>%
    dplyr::filter(celltype %in% CT_ORDER)

  S <- dplyr::bind_rows(lapply(DIRECTIONS, function(dir) {
    d <- calls %>% dplyr::mutate(tested_y = !is.na(pvalue_young), tested_o = !is.na(pvalue_old),
                                 arg_y = is_arg(lfc_young, padj_young, dir), arg_o = is_arg(lfc_old, padj_old, dir))
    dplyr::bind_rows(
      d %>% dplyr::transmute(celltype, gene, set = "Young ARGs", member = arg_y,           bg = tested_y),
      d %>% dplyr::transmute(celltype, gene, set = "Old ARGs",   member = arg_o,           bg = tested_o),
      d %>% dplyr::transmute(celltype, gene, set = "Maintained", member = arg_y & arg_o,   bg = tested_y),
      d %>% dplyr::transmute(celltype, gene, set = "Lost",       member = arg_y & !arg_o,  bg = tested_y),
      d %>% dplyr::transmute(celltype, gene, set = "Gained",     member = arg_o & !arg_y,  bg = tested_o)) %>%
      dplyr::filter(bg) %>% dplyr::select(-bg) %>% dplyr::mutate(direction = dir)
  }))
  sizes <- S %>% dplyr::group_by(direction, celltype, set) %>% dplyr::summarise(n_genes = sum(member), .groups = "drop")
  write.csv(sizes, file.path(OUT, "ARG_SetSizes.csv"), row.names = FALSE)
  write.csv(S %>% dplyr::filter(member) %>% dplyr::select(-member), file.path(OUT, "ARG_SetMembers.csv"), row.names = FALSE)

  # Superpooled GO: union over cell types; background = genes in the anchor table (pooled across cell types)
  f_go <- file.path(OUT, "GO_BP_Superpooled_AllResults.csv")
  if (file.exists(f_go)) { go <- read.csv(f_go, stringsAsFactors = FALSE); cat("Loaded", f_go, "(delete to recompute)\n") } else {
    U <- list(young = unique(calls$gene[!is.na(calls$lfc_young) | !is.na(calls$pvalue_young)]),
              old   = unique(calls$gene[!is.na(calls$lfc_old)   | !is.na(calls$pvalue_old)]))
    keys <- expand.grid(direction = DIRECTIONS, set = SETS, stringsAsFactors = FALSE)
    jobs <- lapply(seq_len(nrow(keys)), function(i) {
      k <- keys[i, ]
      genes <- unique(S$gene[S$direction == k$direction & S$set == k$set & S$member])
      list(meta = data.frame(approach = ap, level = "Superpooled", unit = "All cell types", direction = k$direction, set = k$set),
           genes = genes, universe = U[[SET_ANCHOR[[k$set]]]])
    })
    cat(sprintf("GO: %d gene sets on %d cores...\n", length(jobs), N_CORES))
    go <- go_batch(jobs, n_cores = N_CORES)
    if (nrow(go)) write.csv(go, f_go, row.names = FALSE)
  }
  if (!nrow(go)) { cat("No significant GO terms\n"); next }
  write.csv(go %>% dplyr::filter(p.adjust < GO$fdr) %>%
              dplyr::select(direction, set, ID, Description, GeneRatio, BgRatio, fold, pvalue, p.adjust, geneID),
            file.path(OUT, "GO_BP_Superpooled_Significant.csv"), row.names = FALSE)

  # Fig. 4f: Maintained Up / Down and Lost Up / Down
  ML <- c("Maintained Up", "Maintained Down", "Lost Up", "Lost Down")
  d <- go %>% dplyr::filter(set %in% c("Maintained", "Lost")) %>%
    dplyr::mutate(Cluster = paste(set, tools::toTitleCase(direction))) %>% dplyr::filter(Cluster %in% ML) %>%
    dplyr::arrange(factor(Cluster, levels = ML))
  compare_dot(d, paste0(APPROACH_LABS[[ap]], " | DREADD-induced ARG pathways | Maintained and Lost, Up and Down | Female"),
              file.path(OUT, "SuperPooled_MaintainedLost_UpDown_GO_BP"))
}

cat("\n=== DONE ===\nOutputs:", OUT_DIR, "\n")
