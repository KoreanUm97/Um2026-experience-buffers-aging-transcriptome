################################################################################
### ED10e_BufferingIndex.R
###
### Buffering index per cell type at old age, males and females in separate panels:
###   BI = sum(log2FC OEEYHC^2) / sum(log2FC OHCYHC^2) over genes with a log2 fold change in both
###   (< 1 = old eNE closer to young HC than old HC is; dot size = genes tested).
### Plot: BI Style 1 (lollipop) of Buffering_VizTest.R, as used for the panel.
### Reads the buffering index written by Buffering_Pipeline.R step "bi"
###   (<NE_BASE>/Figures/<date>_Buffering_Index/<sex>_BI_Summary.csv; newest folder is used).
### Output: <NE_BASE>/Paper_Figures/ED10/ED10e_BufferingIndex_<sex>.pdf/.png + source data
################################################################################

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

NE_BASE  <- "PATH_TO_NE_OUTPUT"
BI_GLOB  <- "_Buffering_Index$"            # folder written by the pipeline's bi step
BI_FILE  <- "%s_BI_Summary.csv"            # %s = sex
OUT      <- file.path(NE_BASE, "Paper_Figures", "ED10")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN","Astrocytes","Oligodendrocytes","OPC",
              "Microglia","Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell","Macrophage","Epithelial","T Cell")
CT_COLS <- c("DG ExNs"="#E53935","CA1 ExNs"="#CC0000","CA2 ExNs"="#8E24AA","CA3 ExNs"="#5E35B1","Sub ExNs"="#80D8FF",
             "IN"="#00897B","Astrocytes"="#43A047","Oligodendrocytes"="#1E88E5","OPC"="#00ACC1","Microglia"="#7CB342",
             "Endothelial"="#FB8C00","Pericyte"="#6D4C41","APC"="#546E7A","Stromal"="#27408B","Macrophage"="#C62828",
             "Cajal-Retzius Cell"="#FDD835","T Cell"="#757575")
SEX_LAB <- c(male="Males", female="Females")

dirs <- list.dirs(file.path(NE_BASE, "Figures"), recursive=FALSE)
dirs <- dirs[grepl(BI_GLOB, basename(dirs))]
if (!length(dirs)) stop("No buffering index folder under ", file.path(NE_BASE, "Figures"), ": run Buffering_Pipeline.R step bi")
bi_dir <- dirs[which.max(file.mtime(dirs))]
cat("Buffering index from:", bi_dir, "\n")

base_theme <- theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold", size=13, hjust=0.5),
        plot.subtitle=element_text(size=9, color="grey45", hjust=0.5),
        legend.position="top")

for (sx in names(SEX_LAB)) {
  f <- file.path(bi_dir, sprintf(BI_FILE, sx)); if (!file.exists(f)) { cat("  missing:", f, "\n"); next }
  d <- read.csv(f, stringsAsFactors=FALSE) %>% dplyr::filter(tolower(label) == "old") %>%
    dplyr::mutate(celltype=factor(celltype, levels=rev(c(intersect(CT_ORDER, celltype), sort(setdiff(celltype, CT_ORDER))))))
  write.csv(d, file.path(OUT, paste0("SourceData_ED10e_BufferingIndex_", sx, ".csv")), row.names=FALSE)
  cts <- levels(d$celltype); miss <- setdiff(cts, names(CT_COLS))
  pal <- c(CT_COLS[intersect(names(CT_COLS), cts)], setNames(rep("grey60", length(miss)), miss))
  # BI Style 1 (lollipop) from Buffering_VizTest.R
  p <- ggplot(d, aes(x=bi, y=celltype, color=celltype)) +
    geom_vline(xintercept=1, linetype="dashed", color="grey40", linewidth=0.8) +
    annotate("rect", xmin=-Inf, xmax=1, ymin=-Inf, ymax=Inf, fill="#2166AC", alpha=0.04) +
    geom_segment(aes(x=1, xend=bi, y=celltype, yend=celltype), linewidth=0.7, color="grey70") +
    geom_point(aes(size=n_genes), alpha=0.9) +
    geom_text(aes(label=sprintf("%.3f", bi), hjust=ifelse(bi<1, 1.2, -0.2)), size=2.8, color="grey20") +
    scale_color_manual(values=pal, guide="none") +
    scale_size_continuous(range=c(3,8), name="Genes tested") +
    scale_x_continuous(expand=expansion(mult=c(0.15,0.15))) +
    base_theme +
    labs(title=paste0("Buffering Index - Lollipop  |  ", SEX_LAB[[sx]], " | Old"),
         subtitle="BI < 1 = eNE compresses aging transcriptional distance\nDot size = genes tested",
         x=expression("Buffering Index  (sum OEEYHC"^2*" / sum OHCYHC"^2*")"), y=NULL)
  h <- max(5, nrow(d) * 0.45 + 3)
  ggsave(file.path(OUT, paste0("ED10e_BufferingIndex_", sx, ".pdf")), p, width=9, height=h)
  ggsave(file.path(OUT, paste0("ED10e_BufferingIndex_", sx, ".png")), p, width=9, height=h, dpi=600, bg="white")
  cat("  saved: ED10e_BufferingIndex_", sx, "\n", sep="")
}
