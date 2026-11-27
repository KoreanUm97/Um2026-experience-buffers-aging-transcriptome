################################################################################
### AgeEmergent_ARG_Pathway.R
###
### GO/KEGG pathway enrichment for age-emergent novel EE genes.
###   Novel Up:   sig in MEEMHC or OEEOHC, NOT in YEEYHC, log2FC > 0
###   Novel Down: sig in MEEMHC or OEEOHC, NOT in YEEYHC, log2FC < 0
###
### VISUALIZATION: Gene-level Sankey + Dot plot (SRplot style) for per-CT
###                and super-pooled outputs.
###                compareCluster dot plot for pooled cross-CT comparisons.
###
### Reads DESeq2_Results. GO universe per comparison = genes tested in YEEYHC OR in that aged comparison,
###   so age-emergent genes not tested in young are kept (same gene sets as the ED8a stacked bar).
### Output: Paper_Figures/Figure5/ED8bcd_AgeEmergent_ARG_Pathway/
################################################################################

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

suppressPackageStartupMessages({
  library(conflicted)
  library(ggalluvial)
  library(patchwork)
  conflicts_prefer(dplyr::select,           .quiet=TRUE)
  conflicts_prefer(dplyr::filter,           .quiet=TRUE)
  conflicts_prefer(base::intersect,         .quiet=TRUE)
  conflicts_prefer(base::setdiff,           .quiet=TRUE)
  conflicts_prefer(base::union,             .quiet=TRUE)
  conflicts_prefer(clusterProfiler::select, .quiet=TRUE)
  conflicts_prefer(AnnotationDbi::select,   .quiet=TRUE)
  conflicts_prefer(AnnotationDbi::as.data.frame)
  conflicts_prefer(base::as.data.frame)
})

################################################################################
### PATHS
################################################################################

base_path <- "PATH_TO_NE_OUTPUT/"
deseq_dir <- paste0(base_path,"DESeq2_Results/")
fig_date  <- format(Sys.Date(),"%Y-%m-%d")
out_dir   <- paste0(base_path,"Paper_Figures/Figure5/ED8bcd_AgeEmergent_ARG_Pathway/")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
cat("Output:", out_dir, "\n\n")

PADJ_THR         <- 0.05
Q_THR            <- 0.2
MIN_GENES        <- 5
TOP_N            <- 10
MAX_GENES_PER_PATH <- 6   # cap genes shown per pathway in Sankey (tune: 5=tight, 8=looser)
COMPS            <- c("MEEMHC","OEEOHC")

should_skip <- function(f) {
  if (file.exists(f)) { cat("  Skipping:", basename(f),"\n"); return(TRUE) }
  FALSE
}

################################################################################
### PALETTE
################################################################################

pathway_palette <- function(n) {
  base_cols <- c(
    "#E07B54","#5B8DB8","#6BAF6B","#8F7BB5","#C0B45A",
    "#7BBFBF","#C97BAB","#A07850","#7090B8","#90B870",
    "#B07070","#70A0A0","#A07090","#90A070","#A09060")
  if (n <= length(base_cols)) return(base_cols[seq_len(n)])
  scales::hue_pal()(n)
}

################################################################################
### GENE-LEVEL SANKEY + DOT PLOT (SRplot style)
################################################################################

sankey_dot_plot <- function(enrich_res, title, out_file,
                             top_n              = TOP_N,
                             max_genes_per_path = MAX_GENES_PER_PATH,
                             width    = 14,
                             height   = NULL,
                             low_col  = "#3182BD",
                             high_col = "#E6550D") {

  if (is.null(enrich_res) || nrow(enrich_res) == 0) return(invisible(NULL))
  if (should_skip(out_file)) return(invisible(NULL))

  tryCatch({
    df <- as.data.frame(enrich_res) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::mutate(
        GeneRatio_num = sapply(GeneRatio, function(x) {
          pts <- strsplit(x, "/")[[1]]
          as.numeric(pts[1]) / as.numeric(pts[2])
        }),
        neg_log10_padj = -log10(p.adjust),
        Desc_short = ifelse(
          nchar(Description) > 45,
          paste0(substr(Description, 1, 42), "..."),
          Description))

    pathway_levels <- rev(df$Desc_short)
    df$Desc_short  <- factor(df$Desc_short, levels = pathway_levels)

    # Gene-to-pathway edge table — capped to max_genes_per_path per pathway
    edge_df <- df %>%
      dplyr::select(Desc_short, geneID, neg_log10_padj) %>%
      dplyr::mutate(genes = strsplit(as.character(geneID), "/")) %>%
      tidyr::unnest(genes) %>%
      dplyr::rename(Gene = genes, Pathway = Desc_short) %>%
      dplyr::mutate(Pathway = factor(Pathway, levels = levels(df$Desc_short))) %>%
      # ── cap representative genes per pathway ────────────────────────────────
      dplyr::group_by(Pathway) %>%
      dplyr::slice_head(n = max_genes_per_path) %>%
      dplyr::ungroup()
      # ────────────────────────────────────────────────────────────────────────

    n_paths     <- nlevels(df$Desc_short)
    path_colors <- setNames(pathway_palette(n_paths), levels(df$Desc_short))

    # Gene ordering: by rank of most significant pathway
    gene_rank <- edge_df %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(top_rank = min(as.integer(Pathway)), .groups="drop") %>%
      dplyr::arrange(top_rank, Gene)
    gene_levels  <- unique(gene_rank$Gene)
    edge_df$Gene <- factor(edge_df$Gene, levels = gene_levels)
    n_genes      <- length(gene_levels)

    # LEFT: Gene-level Sankey
    # fill=Pathway in both geom_alluvium and geom_stratum:
    #   flows colored by pathway; gene strata split-colored by pathway mix
    p_sankey <- ggplot(
      edge_df,
      aes(axis1 = Gene, axis2 = Pathway, y = 1)) +
      ggalluvial::geom_alluvium(
        aes(fill = Pathway),
        width=0.25, alpha=0.5, curve_type="sigmoid", knot.pos=0.4) +
      ggalluvial::geom_stratum(
        aes(fill = Pathway),
        width=0.25, color="grey70", linewidth=0.25) +
      geom_text(
        stat=    "stratum",
        aes(label=after_stat(stratum)),
        size=2.6, hjust=1, nudge_x=-0.15,
        color="grey15", fontface="italic") +
      scale_fill_manual(values=path_colors, guide="none") +
      scale_x_discrete(limits=c("Gene","Pathway"), expand=c(0.35,0.05)) +
      theme_void(base_size=9) +
      theme(
        axis.text.x = element_text(size=9, face="bold", color="grey30",
                                   margin=margin(t=4)),
        plot.margin = margin(10,0,10,65))

    # RIGHT: Dot plot with color strip
    dot_xmin <- min(df$GeneRatio_num) * 0.75
    strip_w  <- (max(df$GeneRatio_num) - min(df$GeneRatio_num)) * 0.04

    p_dot <- ggplot(df, aes(x=GeneRatio_num, y=Desc_short)) +
      geom_tile(
        aes(x=dot_xmin - strip_w/2, fill=Desc_short),
        width=strip_w, height=0.75, color=NA) +
      geom_point(aes(color=neg_log10_padj, size=Count), alpha=0.9) +
      scale_color_gradient(
        low=low_col, high=high_col,
        name=bquote(-log[10](italic(p)[adj]))) +
      scale_fill_manual(values=path_colors, guide="none") +
      scale_size_continuous(range=c(2,9), name="Count") +
      scale_x_continuous(
        expand=expansion(mult=c(0.15,0.18)),
        labels=scales::number_format(accuracy=0.01)) +
      theme_classic(base_size=9) +
      theme(
        axis.text.y     = element_text(size=8, color="grey15"),
        axis.ticks.y    = element_blank(),
        axis.line.y     = element_blank(),
        axis.text.x     = element_text(size=7.5),
        axis.title.x    = element_text(size=8.5),
        plot.margin     = margin(10,12,10,4),
        legend.position = "right",
        legend.text     = element_text(size=7),
        legend.title    = element_text(size=7.5)) +
      labs(x="Gene Ratio", y=NULL)

    auto_h <- if (is.null(height)) max(5, n_genes*0.22+2) else height

    combined <- p_sankey + p_dot +
      patchwork::plot_layout(widths=c(2.2,1)) +
      patchwork::plot_annotation(
        title    = title,
        subtitle = paste0("Top ",nrow(df)," terms | ordered by p.adjust",
                          " | max ",max_genes_per_path," genes/pathway shown"),
        theme    = theme(
          plot.title      = element_text(face="bold",size=12,hjust=0.5),
          plot.subtitle   = element_text(size=8.5,color="grey40",hjust=0.5),
          plot.background = element_rect(fill="white",color=NA)))

    ggsave(out_file, combined,
           width=width, height=auto_h,
           device="pdf", limitsize=FALSE)
    cat("  Saved:", basename(out_file), "\n")

  }, error=function(e) cat("  sankey_dot error:", e$message, "\n"))
}

################################################################################
### COMPARECLUSTER DOT PLOT (multi-cluster: keep standard dotplot)
################################################################################

run_compare <- function(gene_lists, universe, ont="BP",
                         title, out_file, top_n=TOP_N,
                         width=14, height=10) {
  gene_lists <- gene_lists[lengths(gene_lists) >= MIN_GENES]
  if (length(gene_lists) < 2) {
    cat("  Too few gene sets\n"); return(invisible(NULL))
  }
  if (should_skip(out_file)) return(invisible(NULL))
  tryCatch({
    res <- compareCluster(
      geneClusters=gene_lists, fun="enrichGO",
      OrgDb=org.Mm.eg.db, keyType="SYMBOL", ont=ont,
      pAdjustMethod="BH", pvalueCutoff=0.05,
      universe=universe, readable=TRUE)
    if (is.null(res)||nrow(res)==0) {
      cat("  No significant terms\n"); return(invisible(NULL))
    }
    write.csv(as.data.frame(res),
              gsub(".pdf",".csv",out_file,fixed=TRUE), row.names=FALSE)
    n_cls <- length(unique(res@compareClusterResult$Cluster))
    p <- dotplot(res, showCategory=top_n) +
      labs(title=title) +
      theme(plot.title  = element_text(face="bold",size=12),
            axis.text.x = element_text(angle=35,hjust=1,size=8),
            axis.text.y = element_text(size=8))
    ggsave(out_file, p,
           width=max(width, n_cls*1.0+3.5), height=height, device="pdf")
    cat("  Saved compareCluster:", basename(out_file), "\n")
    return(res)
  }, error=function(e) {
    cat("  compareCluster error:", e$message, "\n"); NULL
  })
}

################################################################################
### ENRICHMENT HELPERS
################################################################################

run_go <- function(genes, universe, ont="BP") {
  genes <- intersect(genes, universe)
  if (length(genes) < MIN_GENES) return(NULL)
  tryCatch({
    res <- enrichGO(
      gene=genes, universe=universe, OrgDb=org.Mm.eg.db,
      keyType="SYMBOL", ont=ont, pAdjustMethod="BH",
      pvalueCutoff=0.05, qvalueCutoff=Q_THR, readable=TRUE)
    if (is.null(res)||nrow(res)==0) return(NULL)
    clusterProfiler::simplify(res, cutoff=0.7, by="p.adjust")
  }, error=function(e) NULL)
}

sym_to_entrez <- function(genes) {
  res <- AnnotationDbi::mapIds(org.Mm.eg.db, keys=genes,
                                column="ENTREZID", keytype="SYMBOL",
                                multiVals="first")
  as.character(res[!is.na(res)])
}

run_kegg <- function(genes, universe) {
  genes <- intersect(genes, universe)
  if (length(genes) < MIN_GENES) return(NULL)
  eg <- sym_to_entrez(genes); eu <- sym_to_entrez(universe)
  if (length(eg) < MIN_GENES) return(NULL)
  tryCatch({
    res <- enrichKEGG(gene=eg, universe=eu, organism="mmu",
                       pAdjustMethod="BH", pvalueCutoff=0.05,
                       qvalueCutoff=Q_THR)
    if (is.null(res)||nrow(res)==0) return(NULL)
    setReadable(res, OrgDb=org.Mm.eg.db, keyType="ENTREZID")
  }, error=function(e) NULL)
}

################################################################################
### GENE SET BUILDER
################################################################################

build_novel <- function(yee_df, comp_df, ct) {
  y     <- yee_df  %>% dplyr::filter(celltype==ct)
  d     <- comp_df %>% dplyr::filter(celltype==ct)
  y_sig <- y$gene[!is.na(y$padj) & y$padj < PADJ_THR]
  d_up  <- d$gene[!is.na(d$padj) & d$padj < PADJ_THR & d$log2FoldChange > 0]
  d_dn  <- d$gene[!is.na(d$padj) & d$padj < PADJ_THR & d$log2FoldChange < 0]
  list(
    novel_up   = setdiff(d_up, y_sig),
    novel_down = setdiff(d_dn, y_sig),
    novel_all  = setdiff(union(d_up, d_dn), y_sig)
  )
}

################################################################################
### MAIN LOOP
################################################################################

for (current_sex in c("male","female")) {
  cat("\n########################################\n")
  cat("### SEX:", toupper(current_sex), "\n")
  cat("########################################\n\n")

  sex_out <- file.path(out_dir, current_sex)
  dir.create(sex_out, recursive=TRUE, showWarnings=FALSE)

  yee_file <- file.path(deseq_dir,
                         paste0(current_sex,"_ALL_YEEYHC.csv"))
  if (!file.exists(yee_file)) {
    cat("  YEEYHC not found -- skip\n"); next
  }
  yee_all   <- read.csv(yee_file, stringsAsFactors=FALSE)
  universe  <- unique(yee_all$gene[!is.na(yee_all$gene)])
  celltypes <- sort(unique(yee_all$celltype))
  cat(sprintf("Universe: %d genes | Cell types: %d\n",
              length(universe), length(celltypes)))

  for (comp in COMPS) {
    cat(sprintf("\n--- Comparison: %s ---\n", comp))

    comp_file <- file.path(deseq_dir,
                            paste0(current_sex,"_ALL_",comp,".csv"))
    if (!file.exists(comp_file)) {
      cat("  Not found:", basename(comp_file), "\n"); next
    }
    comp_all <- read.csv(comp_file, stringsAsFactors=FALSE)
    # universe for this comparison: genes tested in YEEYHC or in the aged comparison
    # (age-emergent genes absent from YEEYHC would otherwise be dropped by intersect(genes, universe))
    universe <- union(unique(yee_all$gene[!is.na(yee_all$gene)]), unique(comp_all$gene[!is.na(comp_all$gene)]))
    cat(sprintf("  Universe (%s): %d genes\n", comp, length(universe)))
    comp_out <- file.path(sex_out, comp)
    dir.create(comp_out, recursive=TRUE, showWarnings=FALSE)

    pool_up   <- list()
    pool_dn   <- list()
    pool_comb <- list()

    ##########################################################
    ### PER CELL TYPE
    ##########################################################

    cat("\n  Per cell type:\n")

    for (ct in celltypes) {
      gs      <- build_novel(yee_all, comp_all, ct)
      ct_slug <- gsub("[^A-Za-z0-9]+","_",ct)
      ct_out  <- file.path(comp_out, ct_slug)
      dir.create(ct_out, recursive=TRUE, showWarnings=FALSE)

      cat(sprintf("  %-22s  novel_up=%d  novel_dn=%d\n",
                  ct, length(gs$novel_up), length(gs$novel_down)))

      # Accumulate for pooled
      if (length(gs$novel_up)   >= MIN_GENES)
        pool_up[[paste0(ct,"\nNovel Up")]]   <- gs$novel_up
      if (length(gs$novel_down) >= MIN_GENES)
        pool_dn[[paste0(ct,"\nNovel Down")]] <- gs$novel_down
      if (length(gs$novel_all)  >= MIN_GENES)
        pool_comb[[paste0(ct,"\nNovel All")]] <- gs$novel_all

      sets <- list(
        list(genes=gs$novel_up,   slug="NovelUp",   lbl="Novel Up"),
        list(genes=gs$novel_down, slug="NovelDown",  lbl="Novel Down"),
        list(genes=gs$novel_all,  slug="NovelAll",   lbl="Novel (Up+Down)")
      )

      for (s in sets) {
        if (length(s$genes) < MIN_GENES) next

        # GO BP — Sankey+Dot
        out_go <- file.path(ct_out,
          paste0(fig_date,"_",ct_slug,"_",s$slug,"_GO_BP_SankeyDot.pdf"))
        if (!should_skip(out_go)) {
          r <- run_go(s$genes, universe, "BP")
          if (!is.null(r)) {
            write.csv(as.data.frame(r),
                      gsub(".pdf",".csv",out_go,fixed=TRUE), row.names=FALSE)
            sankey_dot_plot(r,
              paste0(ct," | ",s$lbl," | GO BP | ",comp),
              out_go, top_n=TOP_N, max_genes_per_path=MAX_GENES_PER_PATH)
          }
        }

        # KEGG — Sankey+Dot
        out_kk <- file.path(ct_out,
          paste0(fig_date,"_",ct_slug,"_",s$slug,"_KEGG_SankeyDot.pdf"))
        if (!should_skip(out_kk)) {
          rk <- run_kegg(s$genes, universe)
          if (!is.null(rk)) {
            write.csv(as.data.frame(rk),
                      gsub(".pdf",".csv",out_kk,fixed=TRUE), row.names=FALSE)
            sankey_dot_plot(rk,
              paste0(ct," | ",s$lbl," | KEGG | ",comp),
              out_kk, top_n=TOP_N, max_genes_per_path=MAX_GENES_PER_PATH)
          }
        }
      }
    }

    ##########################################################
    ### POOLED COMPARECLUSTER (standard dotplot — multi-cluster)
    ##########################################################

    cat(sprintf("\n  Pooled compareCluster (%s):\n", comp))

    run_compare(
      pool_up, universe, ont="BP",
      title    = paste0("GO BP -- Novel Up per cell type | ",
                         tools::toTitleCase(current_sex)," | ",comp),
      out_file = file.path(sex_out,
        paste0(fig_date,"_",current_sex,"_",comp,
               "_UP_compareCluster_GO_BP.pdf")),
      top_n=10, width=max(14,length(pool_up)*1.0+3.5), height=10)

    run_compare(
      pool_dn, universe, ont="BP",
      title    = paste0("GO BP -- Novel Down per cell type | ",
                         tools::toTitleCase(current_sex)," | ",comp),
      out_file = file.path(sex_out,
        paste0(fig_date,"_",current_sex,"_",comp,
               "_DOWN_compareCluster_GO_BP.pdf")),
      top_n=10, width=max(14,length(pool_dn)*1.0+3.5), height=10)

    run_compare(
      pool_comb, universe, ont="BP",
      title    = paste0("GO BP -- Novel All (Up+Down) per cell type | ",
                         tools::toTitleCase(current_sex)," | ",comp),
      out_file = file.path(sex_out,
        paste0(fig_date,"_",current_sex,"_",comp,
               "_Combined_compareCluster_GO_BP.pdf")),
      top_n=10, width=max(14,length(pool_comb)*1.0+3.5), height=10)

    ##########################################################
    ### SUPER-POOLED — Sankey+Dot
    ##########################################################

    cat(sprintf("\n  Super-pooled (%s):\n", comp))

    super <- list(
      "Novel Up"   = unique(unlist(lapply(celltypes, function(ct)
                       build_novel(yee_all,comp_all,ct)$novel_up))),
      "Novel Down" = unique(unlist(lapply(celltypes, function(ct)
                       build_novel(yee_all,comp_all,ct)$novel_down))),
      "Novel All"  = unique(unlist(lapply(celltypes, function(ct)
                       build_novel(yee_all,comp_all,ct)$novel_all)))
    )

    for (nm in names(super)) {
      genes   <- super[[nm]]
      nm_slug <- gsub("[^A-Za-z0-9]+","_",nm)
      if (length(genes) < MIN_GENES) next
      cat(sprintf("    %s: %d genes\n", nm, length(genes)))

      # GO BP Sankey+Dot
      out_go <- file.path(sex_out,
        paste0(fig_date,"_",current_sex,"_",comp,
               "_SuperPooled_",nm_slug,"_GO_BP_SankeyDot.pdf"))
      if (!should_skip(out_go)) {
        r <- run_go(genes, universe, "BP")
        if (!is.null(r)) {
          write.csv(as.data.frame(r),
                    gsub(".pdf",".csv",out_go,fixed=TRUE), row.names=FALSE)
          sankey_dot_plot(r,
            paste0("Super-pooled: ",nm," | ",current_sex," | ",comp),
            out_go, top_n=TOP_N, max_genes_per_path=MAX_GENES_PER_PATH,
            width=13)
        }
      }

      # KEGG Sankey+Dot
      out_kk <- file.path(sex_out,
        paste0(fig_date,"_",current_sex,"_",comp,
               "_SuperPooled_",nm_slug,"_KEGG_SankeyDot.pdf"))
      if (!should_skip(out_kk)) {
        rk <- run_kegg(genes, universe)
        if (!is.null(rk)) {
          write.csv(as.data.frame(rk),
                    gsub(".pdf",".csv",out_kk,fixed=TRUE), row.names=FALSE)
          sankey_dot_plot(rk,
            paste0("Super-pooled: ",nm," | ",current_sex," | ",comp," | KEGG"),
            out_kk, top_n=TOP_N, max_genes_per_path=MAX_GENES_PER_PATH,
            width=13)
        }
      }
    }

    # Super-pooled compareCluster: Novel Up vs Down vs All
    super_filt <- super[lengths(super) >= MIN_GENES]
    if (length(super_filt) >= 2) {
      run_compare(
        super_filt, universe, ont="BP",
        title    = paste0("GO BP -- Super-pooled Novel ARGs | ",
                           tools::toTitleCase(current_sex)," | ",comp),
        out_file = file.path(sex_out,
          paste0(fig_date,"_",current_sex,"_",comp,
                 "_SuperPooled_compareCluster_GO_BP.pdf")),
        top_n=12, width=12, height=11)
    }
  }
}

cat("\n=== DONE ===\nOutputs:", out_dir, "\n\n")
cat("Per CT outputs (_SankeyDot.pdf):\n")
cat("  NovelUp_GO_BP_SankeyDot    / NovelUp_KEGG_SankeyDot\n")
cat("  NovelDown_GO_BP_SankeyDot  / NovelDown_KEGG_SankeyDot\n")
cat("  NovelAll_GO_BP_SankeyDot   / NovelAll_KEGG_SankeyDot\n")
cat("Super-pooled outputs (_SankeyDot.pdf):\n")
cat("  SuperPooled_Novel_Up/Down/All_GO_BP + KEGG\n")
cat("CompareCluster outputs (standard dotplot):\n")
cat("  _UP/DOWN/Combined/SuperPooled_compareCluster_GO_BP.pdf\n")
