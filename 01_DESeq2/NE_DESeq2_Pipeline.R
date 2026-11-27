
################################################################################
### NE_DESeq2_Pipeline.R
### Pseudobulk DESeq2 for the NE dataset (36 mice; LogNorm integrated object), each sex separately.
### Interneurons merged into one cluster (PV IN + SST IN + VIP IN + Lamp5 IN = "IN").
### Pseudobulk = counts summed per cell type x mouse; DESeq2 ~ group; genes with >= 10 counts in >= 2 samples.
### Comparisons: XEEXHC (eNE vs HC within age), age comparisons within HC and within eNE, and
###   OEEYHC / MEEYHC (aged eNE vs young HC, used for the buffering index).
### Skip-if-exists throughout for safe reruns.
### Output: <base_path>/DESeq2_Results/<sex>_ALL_<comparison>.csv (+ per cell type CSVs, DEG summary, volcano plots)
################################################################################

library(conflicted)
library(BPCells)
library(Seurat)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyr)
library(tibble)

conflicts_prefer(base::intersect,     .quiet=TRUE)
conflicts_prefer(base::setdiff,       .quiet=TRUE)
conflicts_prefer(base::as.data.frame, .quiet=TRUE)
conflicts_prefer(dplyr::select,       .quiet=TRUE)
conflicts_prefer(dplyr::filter,       .quiet=TRUE)
conflicts_prefer(dplyr::arrange,      .quiet=TRUE)
conflicts_prefer(dplyr::mutate,       .quiet=TRUE)
conflicts_prefer(dplyr::pull,         .quiet=TRUE)
conflicts_prefer(DESeq2::vst,         .quiet=TRUE)
conflicts_prefer(dplyr::desc,         .quiet=TRUE)
conflicts_prefer(dplyr::rename,       .quiet=TRUE)

################################################################################
### SETUP
################################################################################

NE_OBJECT <- "PATH_TO_NE_OBJECT.rds"        # integrated NE Seurat object (celltype, orig.ident, experimental.condition, age, sex)
if (!exists("df_harmony")) df_harmony <- readRDS(NE_OBJECT)

cat("Cell types:"); print(table(df_harmony$celltype))
cat("Sex:");        print(table(df_harmony$sex))
cat("Age:");        print(table(df_harmony$age))
cat("Condition:");  print(table(df_harmony$experimental.condition))

if (any(c("PV IN","SST IN","VIP IN","Lamp5 IN") %in% df_harmony$celltype))
  warning("Separate IN subtype labels still present -- merge them into \"IN\" first")

base_path    <- "PATH_TO_NE_OUTPUT/"
deseq_path   <- paste0(base_path,"DESeq2_Results/")
volcano_path <- paste0(base_path,"Volcano_Plots/")
dir.create(deseq_path,   recursive=TRUE, showWarnings=FALSE)
dir.create(volcano_path, recursive=TRUE, showWarnings=FALSE)

PADJ_THRESHOLD <- 0.05
YOUNG_LABEL    <- "young"
OLD_LABEL      <- "old"

age_vals     <- unique(tolower(df_harmony@meta.data$age))
cat("Age labels in object:", paste(sort(age_vals),collapse=", "), "\n")
MIDDLE_LABEL <- age_vals[!age_vals %in% c("young","old")][1]
cat("Using middle age label:", MIDDLE_LABEL, "\n")
if (is.na(MIDDLE_LABEL)) stop("Could not detect middle-age label.")

################################################################################
### PSEUDOBULK AGGREGATION
################################################################################

cat("\n=== Aggregating pseudobulk counts ===\n")

pb_all <- suppressWarnings(
  AggregateExpression(
    object               = df_harmony,
    assays               = "RNA",
    return.seurat        = FALSE,
    group.by             = c("celltype","orig.ident",
                              "experimental.condition","age","sex"),
    normalization.method = NULL,
    slot                 = "counts",
    verbose              = TRUE
  )$RNA
)

split_names <- strsplit(colnames(pb_all),"_")
parse_right <- function(parts,from_end) parts[length(parts)-from_end+1]
pb_meta <- data.frame(
  sample     = colnames(pb_all),
  sex        = sapply(split_names,parse_right,from_end=1),
  age        = sapply(split_names,parse_right,from_end=2),
  condition  = sapply(split_names,parse_right,from_end=3),
  orig.ident = sapply(split_names,parse_right,from_end=4),
  celltype   = sapply(split_names,function(x)
                 paste(x[seq_len(length(x)-4)],collapse=" ")),
  stringsAsFactors=FALSE
)
pb_meta$age_norm <- tolower(pb_meta$age)
pb_meta$sex_norm <- tolower(pb_meta$sex)

cat("Pseudobulk samples:",nrow(pb_meta),"\n")
cat("Age x Condition x Sex:\n")
print(table(pb_meta$age_norm,pb_meta$condition,pb_meta$sex_norm))
celltypes_all <- sort(unique(pb_meta$celltype))
cat("Cell types:",paste(celltypes_all,collapse=", "),"\n")

################################################################################
### DESeq2 HELPER (with skip-if-exists)
################################################################################

run_deseq2 <- function(counts_mat, meta_df,
                        group_col, numerator, denominator,
                        celltype, comp_name,
                        out_dir=deseq_path, sex_label=current_sex) {
  # Skip if already done
  out_f <- file.path(out_dir,paste0(sex_label,"_",celltype,"_",comp_name,".csv"))
  if (file.exists(out_f)) {
    cat("  [skip existing]\n")
    return(read.csv(out_f,stringsAsFactors=FALSE))
  }

  ct_keep    <- meta_df$celltype==celltype
  sub_meta   <- meta_df[ct_keep,]
  sub_counts <- counts_mat[,ct_keep,drop=FALSE]

  grp_keep   <- sub_meta[[group_col]] %in% c(numerator,denominator)
  sub_meta   <- sub_meta[grp_keep,]
  sub_counts <- sub_counts[,grp_keep,drop=FALSE]
  if (ncol(sub_counts)==0) return(NULL)

  sub_meta$grp <- factor(sub_meta[[group_col]],
                          levels=c(denominator,numerator))
  grp_tab <- table(sub_meta$grp)
  cat("  Groups:",paste(names(grp_tab),grp_tab,sep="=",collapse=", "),"\n")

  if (length(grp_tab)<2 || any(grp_tab<2)) {
    cat("  Skipping -- need >= 2 samples per group\n"); return(NULL)
  }

  keep_genes <- rowSums(sub_counts>=10)>=2
  sub_counts <- sub_counts[keep_genes,,drop=FALSE]
  if (nrow(sub_counts)<10) {
    cat("  Skipping -- too few genes\n"); return(NULL)
  }

  tryCatch({
    dds    <- DESeqDataSetFromMatrix(sub_counts,sub_meta,design=~grp)
    dds    <- DESeq(dds,quiet=TRUE)
    res    <- results(dds,alpha=PADJ_THRESHOLD)
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::mutate(celltype=celltype,comparison=comp_name)
    n_sig <- sum(!is.na(res_df$padj)&res_df$padj<PADJ_THRESHOLD)
    cat("  DEGs:",n_sig,"\n")
    write.csv(res_df,out_f,row.names=FALSE)
    res_df
  }, error=function(e) { cat("  Error:",e$message,"\n"); NULL })
}

################################################################################
### CROSS-AGE HELPER (for OEEYHC, MEEYHC -- crosses condition and age)
################################################################################

run_cross_age <- function(pb_sex, pb_meta_sex, mask,
                           grp_col, ref_label, num_label,
                           comp_name, sex_label, celltypes) {
  all_f <- file.path(deseq_path,paste0(sex_label,"_ALL_",comp_name,".csv"))
  if (file.exists(all_f)) {
    cat("  ALL file exists -- loading\n")
    return(split(read.csv(all_f,stringsAsFactors=FALSE),
                 read.csv(all_f,stringsAsFactors=FALSE)$celltype))
  }

  sub_pb   <- pb_sex[,mask,drop=FALSE]
  sub_meta <- pb_meta_sex[mask,]
  sub_meta[[grp_col]] <- factor(sub_meta[[grp_col]],
                                  levels=c(ref_label,num_label))

  comp_results <- list()
  for (ct in celltypes) {
    ct_f <- file.path(deseq_path,paste0(sex_label,"_",ct,"_",comp_name,".csv"))
    if (file.exists(ct_f)) {
      cat(sprintf("  %-22s [skip]\n",ct))
      comp_results[[ct]] <- read.csv(ct_f,stringsAsFactors=FALSE)
      next
    }
    cat(sprintf("  %-22s: ",ct))
    ct_k  <- sub_meta$celltype==ct
    sm    <- sub_meta[ct_k,]
    sc    <- sub_pb[,ct_k,drop=FALSE]
    if (ncol(sc)==0) { cat("no samples\n"); next }
    gt    <- table(sm[[grp_col]])
    cat(paste(names(gt),gt,sep="=",collapse=" | ")," ")
    if (length(gt)<2||any(gt<2)) { cat("insufficient\n"); next }
    kg <- rowSums(sc>=10)>=2; sc <- sc[kg,,drop=FALSE]
    if (nrow(sc)<10) { cat("few genes\n"); next }
    res <- tryCatch({
      dds <- DESeqDataSetFromMatrix(sc,sm,design=as.formula(paste0("~",grp_col)))
      dds <- DESeq(dds,quiet=TRUE)
      r   <- results(dds,alpha=PADJ_THRESHOLD)
      as.data.frame(r) %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::mutate(celltype=ct,comparison=comp_name)
    },error=function(e){cat("Error:",e$message,"\n");NULL})
    if (!is.null(res)) {
      n_sig <- sum(!is.na(res$padj)&res$padj<PADJ_THRESHOLD)
      cat(sprintf("%d genes, %d sig\n",nrow(res),n_sig))
      write.csv(res,ct_f,row.names=FALSE)
      comp_results[[ct]] <- res
    } else cat("\n")
  }
  if (length(comp_results)>0) {
    all_df <- dplyr::bind_rows(comp_results)
    write.csv(all_df,all_f,row.names=FALSE)
    cat(sprintf("  Saved ALL_%s: %d rows\n",comp_name,nrow(all_df)))
  }
  comp_results
}

################################################################################
### VOLCANO PLOT FUNCTIONS
################################################################################

create_volcano <- function(deg_df, comparison_name, celltype,
                            out_dir, sex_label) {
  out_f <- file.path(out_dir,paste0(sex_label,"_",celltype,"_",
                                     comparison_name,"_Volcano.pdf"))
  if (file.exists(out_f)) return(invisible(NULL))
  df <- deg_df %>%
    dplyr::mutate(
      negLog10padj = -log10(padj),
      significant  = !is.na(padj)&padj<PADJ_THRESHOLD&abs(log2FoldChange)>0.5,
      color = dplyr::case_when(
        significant&log2FoldChange>0.5  ~ "Upregulated",
        significant&log2FoldChange< -0.5 ~ "Downregulated",
        TRUE ~ "Not Significant"))
  top_genes <- df %>% dplyr::filter(significant) %>%
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
    head(50) %>% dplyr::pull(gene)
  label_df <- df[df$gene %in% top_genes,]
  p <- ggplot(df,aes(x=log2FoldChange,y=negLog10padj)) +
    geom_point(aes(color=color),alpha=0.6,size=1.5) +
    scale_color_manual(values=c(Upregulated="red",Downregulated="blue",
                                 "Not Significant"="grey70"),name=NULL) +
    geom_vline(xintercept=c(-0.5,0.5),linetype="dashed",color="grey40") +
    geom_hline(yintercept=-log10(PADJ_THRESHOLD),linetype="dashed",color="grey40") +
    geom_text_repel(data=label_df,aes(label=gene),
                     size=2.5,box.padding=0.5,max.overlaps=50,seed=42) +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(face="bold",size=13)) +
    labs(title=paste0(celltype," (",sex_label,"): ",comparison_name),
         subtitle=paste0("Top 50 labeled | n=",nrow(df)),
         x="log2 Fold Change",y="-log10(adjusted p-value)")
  ggsave(out_f,p,width=10,height=8,device="pdf")
}

create_combined_volcano <- function(deg_list,comparison_name,out_dir,sex_label) {
  out_f <- file.path(out_dir,paste0("CombinedVolcano_",sex_label,"_",
                                     comparison_name,".pdf"))
  if (file.exists(out_f)) return(invisible(NULL))
  if (length(deg_list)==0) return(NULL)
  all_data <- dplyr::bind_rows(lapply(names(deg_list),function(ct)
    deg_list[[ct]] %>% dplyr::mutate(
      celltype=ct,
      sig=dplyr::case_when(
        !is.na(padj)&padj<PADJ_THRESHOLD&log2FoldChange>0 ~ "Up",
        !is.na(padj)&padj<PADJ_THRESHOLD&log2FoldChange<0 ~ "Down",
        TRUE ~ "NS"))))
  counts_df <- all_data %>% dplyr::filter(sig!="NS") %>%
    dplyr::group_by(celltype,sig) %>%
    dplyr::summarise(n=n(),.groups="drop") %>%
    tidyr::pivot_wider(names_from=sig,values_from=n,values_fill=0)
  if (!"Up"   %in% names(counts_df)) counts_df$Up   <- 0L
  if (!"Down" %in% names(counts_df)) counts_df$Down <- 0L
  ymax <- max(abs(all_data$log2FoldChange),na.rm=TRUE)
  y_pad <- ymax*0.20
  p <- ggplot() +
    geom_jitter(data=dplyr::filter(all_data,sig=="NS"),
                aes(x=celltype,y=log2FoldChange),
                width=0.25,alpha=0.4,size=1.2,color="gray70") +
    geom_jitter(data=dplyr::filter(all_data,sig=="Up"),
                aes(x=celltype,y=log2FoldChange),
                width=0.25,alpha=0.85,size=1.8,color="#D6604D") +
    geom_jitter(data=dplyr::filter(all_data,sig=="Down"),
                aes(x=celltype,y=log2FoldChange),
                width=0.25,alpha=0.85,size=1.8,color="#4393C3") +
    geom_hline(yintercept=0,linetype="dashed",color="black",linewidth=0.8) +
    geom_text(data=counts_df,
              aes(x=celltype,y=ymax+y_pad,label=paste0("",Up)),
              color="#D6604D",size=3.5,fontface="bold") +
    geom_text(data=counts_df,
              aes(x=celltype,y=-ymax-y_pad,label=paste0("",Down)),
              color="#4393C3",size=3.5,fontface="bold") +
    coord_cartesian(ylim=c(-ymax-y_pad*2.5,ymax+y_pad*2.5),clip="off") +
    theme_classic(base_size=12) +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=9),
          plot.title=element_text(face="bold",size=13),
          plot.subtitle=element_text(size=10,color="grey40"),
          panel.grid.major.y=element_line(color="gray90",linewidth=0.3)) +
    labs(title=paste0(comparison_name," - ",tools::toTitleCase(sex_label)),
         subtitle="Red=Up (padj<0.05) | Blue=Down | Gray=NS",
         x="Cell Type",y="log2 Fold Change")
  ggsave(out_f,p,
         width=max(8,length(unique(all_data$celltype))*0.5+4),
         height=7,device="pdf")
  cat("  Saved combined volcano:",comparison_name,"\n")
  p
}

################################################################################
### COMPARISONS
################################################################################

EEvHC_comparisons <- list(
  list(name="YEEYHC",  age=YOUNG_LABEL,  num="EE", den="HC"),
  list(name="OEEOHC",  age=OLD_LABEL,    num="EE", den="HC"),
  list(name="MEEMHC",  age=MIDDLE_LABEL, num="EE", den="HC")
)
HC_age_comparisons <- list(
  list(name="OHCYHC", num=OLD_LABEL,    den=YOUNG_LABEL),
  list(name="OHCMHC", num=OLD_LABEL,    den=MIDDLE_LABEL),
  list(name="MHCYHC", num=MIDDLE_LABEL, den=YOUNG_LABEL)
)
EE_age_comparisons <- list(
  list(name="OEEYEE", num=OLD_LABEL,    den=YOUNG_LABEL),
  list(name="OEEMEE", num=OLD_LABEL,    den=MIDDLE_LABEL),
  list(name="MEEYEE", num=MIDDLE_LABEL, den=YOUNG_LABEL)
)

################################################################################
### MAIN LOOP
################################################################################

for (current_sex in c("male","female")) {

  cat(sprintf("\n%s\n### SEX: %s\n%s\n",
              strrep("#",80),toupper(current_sex),strrep("#",80)))

  sex_keep    <- pb_meta$sex_norm==current_sex
  pb_sex      <- pb_all[,sex_keep,drop=FALSE]
  pb_meta_sex <- pb_meta[sex_keep,]
  cat("Samples:",nrow(pb_meta_sex),"\n")
  celltypes <- sort(unique(pb_meta_sex$celltype))

  all_results <- list()

  # -- EE vs HC (within age) -------------------------------------------------
  for (comp in EEvHC_comparisons) {
    comp_name <- comp$name
    all_f     <- file.path(deseq_path,paste0(current_sex,"_ALL_",comp_name,".csv"))
    cat("\n===",comp_name,"===\n")
    if (file.exists(all_f)) {
      cat("  ALL file exists -- skip\n")
      all_results[[comp_name]] <- split(
        read.csv(all_f,stringsAsFactors=FALSE),
        read.csv(all_f,stringsAsFactors=FALSE)$celltype)
      next
    }
    age_keep    <- pb_meta_sex$age_norm==tolower(comp$age)
    pb_age      <- pb_sex[,age_keep,drop=FALSE]
    pb_meta_age <- pb_meta_sex[age_keep,]
    comp_results <- list()
    for (ct in celltypes) {
      cat("  Celltype:",ct,"\n")
      res <- run_deseq2(pb_age,pb_meta_age,"condition",
                         comp$num,comp$den,ct,comp_name)
      if (!is.null(res)) comp_results[[ct]] <- res
    }
    if (length(comp_results)>0) {
      all_results[[comp_name]] <- comp_results
      write.csv(dplyr::bind_rows(comp_results),all_f,row.names=FALSE)
    }
  }

  # -- HC age comparisons ----------------------------------------------------
  hc_keep    <- pb_meta_sex$condition=="HC"
  pb_hc      <- pb_sex[,hc_keep,drop=FALSE]
  pb_meta_hc <- pb_meta_sex[hc_keep,]
  for (comp in HC_age_comparisons) {
    comp_name <- comp$name
    all_f     <- file.path(deseq_path,paste0(current_sex,"_ALL_",comp_name,".csv"))
    cat("\n===",comp_name,"===\n")
    if (file.exists(all_f)) {
      cat("  ALL file exists -- skip\n")
      all_results[[comp_name]] <- split(
        read.csv(all_f,stringsAsFactors=FALSE),
        read.csv(all_f,stringsAsFactors=FALSE)$celltype)
      next
    }
    comp_results <- list()
    for (ct in celltypes) {
      cat("  Celltype:",ct,"\n")
      res <- run_deseq2(pb_hc,pb_meta_hc,"age_norm",
                         tolower(comp$num),tolower(comp$den),ct,comp_name)
      if (!is.null(res)) comp_results[[ct]] <- res
    }
    if (length(comp_results)>0) {
      all_results[[comp_name]] <- comp_results
      write.csv(dplyr::bind_rows(comp_results),all_f,row.names=FALSE)
    }
  }

  # -- EE age comparisons ----------------------------------------------------
  ee_keep    <- pb_meta_sex$condition=="EE"
  pb_ee      <- pb_sex[,ee_keep,drop=FALSE]
  pb_meta_ee <- pb_meta_sex[ee_keep,]
  for (comp in EE_age_comparisons) {
    comp_name <- comp$name
    all_f     <- file.path(deseq_path,paste0(current_sex,"_ALL_",comp_name,".csv"))
    cat("\n===",comp_name,"===\n")
    if (file.exists(all_f)) {
      cat("  ALL file exists -- skip\n")
      all_results[[comp_name]] <- split(
        read.csv(all_f,stringsAsFactors=FALSE),
        read.csv(all_f,stringsAsFactors=FALSE)$celltype)
      next
    }
    comp_results <- list()
    for (ct in celltypes) {
      cat("  Celltype:",ct,"\n")
      res <- run_deseq2(pb_ee,pb_meta_ee,"age_norm",
                         tolower(comp$num),tolower(comp$den),ct,comp_name)
      if (!is.null(res)) comp_results[[ct]] <- res
    }
    if (length(comp_results)>0) {
      all_results[[comp_name]] <- comp_results
      write.csv(dplyr::bind_rows(comp_results),all_f,row.names=FALSE)
    }
  }

  # -- OEEYHC: Old EE vs Young HC (cross-age, cross-condition) -------------
  cat("\n=== OEEYHC: Old EE vs Young HC ===\n")
  oeeyhc_mask <- (pb_meta_sex$age_norm==YOUNG_LABEL & pb_meta_sex$condition=="HC") |
                  (pb_meta_sex$age_norm==OLD_LABEL   & pb_meta_sex$condition=="EE")
  pb_meta_sex_oeeyhc <- pb_meta_sex
  pb_meta_sex_oeeyhc$grp_oeeyhc <- ifelse(
    pb_meta_sex_oeeyhc$age_norm==OLD_LABEL,"OEE","YHC")
  res_oeeyhc <- run_cross_age(
    pb_sex, pb_meta_sex_oeeyhc, oeeyhc_mask,
    "grp_oeeyhc","YHC","OEE","OEEYHC",current_sex,celltypes)
  if (length(res_oeeyhc)>0) all_results[["OEEYHC"]] <- res_oeeyhc

  # -- MEEYHC: Middle EE vs Young HC (cross-age, cross-condition) ----------
  # same structure as OEEYHC, middle age instead of old
  cat("\n=== MEEYHC: Middle EE vs Young HC ===\n")
  meeyhc_mask <- (pb_meta_sex$age_norm==YOUNG_LABEL   & pb_meta_sex$condition=="HC") |
                  (pb_meta_sex$age_norm==MIDDLE_LABEL & pb_meta_sex$condition=="EE")
  pb_meta_sex_meeyhc <- pb_meta_sex
  pb_meta_sex_meeyhc$grp_meeyhc <- ifelse(
    pb_meta_sex_meeyhc$age_norm==MIDDLE_LABEL,"MEE","YHC")
  res_meeyhc <- run_cross_age(
    pb_sex, pb_meta_sex_meeyhc, meeyhc_mask,
    "grp_meeyhc","YHC","MEE","MEEYHC",current_sex,celltypes)
  if (length(res_meeyhc)>0) all_results[["MEEYHC"]] <- res_meeyhc

  # -- VOLCANOS --------------------------------------------------------------
  cat("\n=== Generating volcano plots for",toupper(current_sex),"===\n")
  all_comp_names <- c(sapply(EEvHC_comparisons,"[[","name"),
                       sapply(HC_age_comparisons,"[[","name"),
                       sapply(EE_age_comparisons,"[[","name"))
  volcano_plots <- list()
  for (comp_name in all_comp_names) {
    if (!comp_name %in% names(all_results)) next
    for (ct in names(all_results[[comp_name]]))
      create_volcano(all_results[[comp_name]][[ct]],comp_name,ct,
                      volcano_path,current_sex)
    p <- create_combined_volcano(all_results[[comp_name]],comp_name,
                                  volcano_path,current_sex)
    if (!is.null(p)) volcano_plots[[comp_name]] <- p
  }
  if (length(volcano_plots)>0) {
    pdf_f <- file.path(volcano_path,
      paste0(current_sex,"_AllComparisons_Combined_Volcano.pdf"))
    if (!file.exists(pdf_f)) {
      pdf(pdf_f,width=max(10,length(celltypes)*0.5+4),
          height=7*length(volcano_plots))
      for (p in volcano_plots) print(p)
      dev.off()
    }
  }

  # -- DEG SUMMARY -----------------------------------------------------------
  all_sex_df <- dplyr::bind_rows(lapply(all_results,dplyr::bind_rows))
  if (nrow(all_sex_df)>0) {
    summary_wide <- all_sex_df %>%
      dplyr::filter(!is.na(padj)&padj<PADJ_THRESHOLD) %>%
      dplyr::group_by(celltype,comparison) %>%
      dplyr::summarise(n_total=n(),n_up=sum(log2FoldChange>0),
                        n_down=sum(log2FoldChange<0),.groups="drop") %>%
      tidyr::pivot_wider(names_from=comparison,values_from=n_total,values_fill=0)
    cat("\nDEG Summary --",toupper(current_sex),":\n")
    print(summary_wide)
    write.csv(summary_wide,
              file.path(deseq_path,paste0(current_sex,"_DEG_Summary.csv")),
              row.names=FALSE)
  }
  cat("\n--- Completed",toupper(current_sex),"---\n")
}

cat("\n=== DONE ===\n")
cat("Results:", deseq_path, "\n")
cat("\nComparisons now available:\n")
print(list.files(deseq_path,pattern="^male_ALL_",full.names=FALSE))
