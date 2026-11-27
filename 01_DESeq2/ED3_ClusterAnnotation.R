############################################################
# ED3_ClusterAnnotation.R
#
# Extended Data Fig. 3: stacked violin plot of the marker genes used for
# manual cell-type annotation of the eNE dataset.
#
# Needs the integrated eNE Seurat object loaded as df_harmony (deposited in GEO).
#
# Marker-group order:
#   Pan-Excitatory
#   Pan-Neuronal
#   Neural Stem
#   DG ExNs
#   CA3 ExNs
#   CA1 ExNs
#   CA2 ExNs
#   Sub ExNs
#   IN
#   Astrocytes
#   Oligodendrocytes
#   OPC
#   Microglia
#   Endothelial
#   Pericyte
#   APC
#   Stromal
#   Cajal-Retzius Cell
#   Macrophage
#   Epithelial
#   T Cell
#
# Y-axis cell-type order remains separate.
############################################################


############################
# 1. Libraries
############################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)


############################
# 2. Input object
############################

# df_harmony should already be loaded, e.g.
#   df_harmony <- readRDS("PATH_TO_eNE_OBJECT.rds")

BASE <- "PATH_TO_OUTPUT"     # same folder as BASE in Fig2_ED4_Pipeline.R
if (grepl("^PATH_TO_", BASE)) stop("Set BASE at the top of this script")
if (!exists("df_harmony")) stop("Load the integrated eNE object as df_harmony first")
if (!"data" %in% Layers(df_harmony[["RNA"]])) df_harmony <- NormalizeData(df_harmony, verbose = FALSE)   # deposited object has counts only

DefaultAssay(df_harmony) <- "RNA"


############################
# 3. Y-axis cell-type order
############################

# This controls the order of the cell types
# on the LEFT side of the violin plot.
#
# This is kept independent from the marker
# group order below.

CT_ORDER_Y <- c(
  
  "DG ExNs",
  "CA3 ExNs",
  "CA1 ExNs",
  "CA2 ExNs",
  "Sub ExNs",
  "IN",
  "Astrocytes",
  "Oligodendrocytes",
  "OPC",
  "Microglia",
  "Endothelial",
  "Pericyte",
  "APC",
  "Stromal",
  "Cajal-Retzius Cell",
  "Macrophage",
  "Epithelial",
  "T Cell",
  "Pan-Excitatory",
  "Pan-Neuronal",
  "Neural Stem"
)


############################
# 4. Marker-group order
############################

# This controls the order of marker groups
# and genes along the x-axis.

CT_ORDER_MARKERS <- c(
  
  # Broad markers FIRST
  "Pan-Excitatory",
  "Pan-Neuronal",
  "Neural Stem",
  
  # Neuronal
  "DG ExNs",
  "CA3 ExNs",
  "CA1 ExNs",
  "CA2 ExNs",
  "Sub ExNs",
  "IN",
  
  # Glial
  "Astrocytes",
  "Oligodendrocytes",
  "OPC",
  "Microglia",
  
  # Vascular / other
  "Endothelial",
  "Pericyte",
  "APC",
  "Stromal",
  "Cajal-Retzius Cell",
  "Macrophage",
  "Epithelial",
  "T Cell"
)


############################
# 5. Cell-type colors
############################

# Cell-type colors used in every figure

CT_COLORS <- c(
  
  "DG ExNs"            = "#E53935",
  "CA1 ExNs"           = "#CC0000",
  "CA2 ExNs"           = "#8E24AA",
  "CA3 ExNs"           = "#5E35B1",
  "Sub ExNs"           = "#80D8FF",
  "IN"                 = "#00897B",
  "Astrocytes"         = "#43A047",
  "Oligodendrocytes"   = "#1E88E5",
  "OPC"                = "#00ACC1",
  "Microglia"          = "#7CB342",
  "Endothelial"        = "#FB8C00",
  "Pericyte"           = "#6D4C41",
  "APC"                = "#546E7A",
  "Stromal"            = "#27408B",
  "Macrophage"         = "#C62828",
  "Cajal-Retzius Cell" = "#FDD835",
  "T Cell"             = "#757575",
  "Ambiguous"          = "#BDBDBD",
  
  # Additional colors ONLY for marker-only groups
  "Pan-Excitatory"     = "#FF7043",
  "Pan-Neuronal"       = "#26A69A",
  "Neural Stem"        = "#AB47BC",
  "Epithelial"         = "#8D6E63"
)


############################
# 6. Marker genes
############################

marker_genes <- list(
  
  ##########################################################
  # Broad
  ##########################################################
  
  "Pan-Excitatory" = c(
    "Slc17a7",
    "Slc17a6",
    "Cux2"
  ),
  
  "Pan-Neuronal" = c(
    "Rbfox3",
    "Map2",
    "Sv2b"
  ),
  
  "Neural Stem" = c(
    "Sox2",
    "Nestin",
    "Fabp7"
  ),
  
  
  ##########################################################
  # Neuronal
  ##########################################################
  
  "DG ExNs" = c(
    "Prox1",
    "Glis3",
    "Egr3"
  ),
  
  "CA3 ExNs" = c(
    "Cpne4",
    "Spock1"
  ),
  
  "CA1 ExNs" = c(
    "Spink8",
    "Mpped1"
  ),
  
  "CA2 ExNs" = c(
    "Tiam2",
    "Pcp4",
    "Rgs14",
    "Amigo2"
  ),
  
  "Sub ExNs" = c(
    "Tshz2"
  ),
  
  "IN" = c(
    "Gad1",
    "Gad2",
    "Slc32a1"
  ),
  
  
  ##########################################################
  # Glial
  ##########################################################
  
  "Astrocytes" = c(
    "Gfap",
    "Aldh1l1",
    "Aqp4"
  ),
  
  "Oligodendrocytes" = c(
    "Olig1",
    "Olig2",
    "Ptgds"
  ),
  
  "OPC" = c(
    "Pdgfra"
  ),
  
  "Microglia" = c(
    "C1qc",
    "Itgam",
    "Trem2"
  ),
  
  
  ##########################################################
  # Vascular / other
  ##########################################################
  
  "Endothelial" = c(
    "Pecam1",
    "Flt1",
    "Cldn5"
  ),
  
  "Pericyte" = c(
    "Pdgfrb",
    "Kcnj8",
    "Abcc9"
  ),
  
  "APC" = c(
    "S100a6"
  ),
  
  "Stromal" = c(
    "Col1a1",
    "Col1a2",
    "Dcn"
  ),
  
  "Cajal-Retzius Cell" = c(
    "Reln",
    "Trp73",
    "Lhx1"
  ),
  
  "Macrophage" = c(
    "Cd68",
    "Csf1r",
    "Ms4a7"
  ),
  
  "Epithelial" = c(
    "Epcam",
    "Krt8",
    "Krt18"
  ),
  
  "T Cell" = c(
    "Cd3e",
    "Cd3d",
    "Cd247"
  )
)


############################
# 7. Check marker groups
############################

stopifnot(
  all(
    names(marker_genes) %in%
      CT_ORDER_MARKERS
  )
)


############################
# 8. Check duplicate genes
############################

all_marker_genes <- unlist(
  marker_genes,
  use.names = FALSE
)

duplicate_genes <- unique(
  all_marker_genes[
    duplicated(all_marker_genes)
  ]
)

if (length(duplicate_genes) > 0) {
  
  warning(
    paste0(
      "Duplicate marker gene(s) detected: ",
      paste(
        duplicate_genes,
        collapse = ", "
      )
    )
  )
}


############################
# 9. Determine cell types
#    actually present
############################

ct_present <- CT_ORDER_Y[
  CT_ORDER_Y %in%
    unique(df_harmony$celltype)
]

cat("\nCell types present in object:\n")
print(ct_present)


############################
# 10. Build ordered marker table
############################

gene_group_raw <- bind_rows(
  
  lapply(
    
    names(marker_genes),
    
    function(ct) {
      
      data.frame(
        group = ct,
        gene = marker_genes[[ct]],
        stringsAsFactors = FALSE
      )
      
    }
  )
)


# Apply desired marker-group order

gene_group_raw$group <- factor(
  gene_group_raw$group,
  levels = CT_ORDER_MARKERS
)


# Sort by marker-group order

gene_group_raw <- gene_group_raw %>%
  arrange(group)


# Ordered genes

genes <- gene_group_raw$gene


############################
# 11. Check genes in RNA assay
############################

genes_present <- genes[
  genes %in%
    rownames(df_harmony[["RNA"]])
]

genes_missing <- setdiff(
  genes,
  genes_present
)

if (length(genes_missing) > 0) {
  
  warning(
    paste0(
      "The following marker genes were not found ",
      "in the RNA assay:\n",
      paste(
        genes_missing,
        collapse = ", "
      )
    )
  )
}


############################
# 12. Keep only genes present
############################

gene_group <- gene_group_raw %>%
  filter(
    gene %in% genes_present
  )

genes <- gene_group$gene


############################
# 13. Fetch normalized expression
############################

fetch_vars <- c(
  "celltype",
  genes
)

plot_data <- FetchData(
  df_harmony,
  vars = fetch_vars,
  layer = "data"
)


############################
# 14. Restrict to cell types
#     actually present
############################

plot_data <- plot_data %>%
  filter(
    celltype %in% ct_present
  )


############################
# 15. Preserve original
#     Y-axis order
############################

# IMPORTANT:
#
# Do NOT use rev() here.
#
# This preserves the original order
# defined in CT_ORDER_Y.

plot_data$celltype <- factor(
  plot_data$celltype,
  levels = ct_present
)


############################
# 16. Convert to long format
############################

plot_long <- plot_data %>%
  
  pivot_longer(
    cols = all_of(genes),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  
  left_join(
    
    gene_group %>%
      mutate(
        gene = as.character(gene),
        group = as.character(group)
      ),
    
    by = "gene"
  )


############################
# 17. Set gene order
############################

plot_long$gene <- factor(
  plot_long$gene,
  levels = genes
)


plot_long$group <- factor(
  plot_long$group,
  levels = CT_ORDER_MARKERS
)


############################
# 18. Main stacked violin plot
############################

p_violin <- ggplot(
  
  plot_long,
  
  aes(
    x = gene,
    y = expression,
    fill = celltype
  )
) +
  
  geom_violin(
    
    scale = "width",
    
    trim = TRUE,
    
    linewidth = 0.15,
    
    color = "grey30"
  ) +
  
  
  facet_grid(
    
    celltype ~ .,
    
    scales = "free_y",
    
    switch = "y"
  ) +
  
  
  scale_x_discrete(
    drop = FALSE
  ) +
  
  
  # Cell-type colors
  scale_fill_manual(
    
    values = CT_COLORS,
    
    drop = FALSE
  ) +
  
  
  labs(
    
    x = NULL,
    
    y = "Normalized expression"
  ) +
  
  
  theme_classic(
    base_size = 11
  ) +
  
  
  theme(
    
    # Cell-type labels
    strip.placement = "outside",
    
    strip.background = element_blank(),
    
    strip.text.y.left = element_text(
      
      angle = 0,
      
      hjust = 1,
      
      face = "bold",
      
      size = 10
    ),
    
    
    # Gene labels
    axis.text.x = element_text(
      
      angle = 90,
      
      vjust = 0.5,
      
      hjust = 1,
      
      size = 8
    ),
    
    
    axis.title.y = element_text(
      size = 11
    ),
    
    
    # Space between cell-type rows
    panel.spacing.y = unit(
      0.15,
      "lines"
    ),
    
    
    panel.grid = element_blank(),
    
    
    # No legend because cell-type labels
    # are already displayed on the left
    legend.position = "none",
    
    
    plot.margin = margin(
      5,
      10,
      5,
      10
    )
  )


############################
# 19. Gene x-position geometry
############################

gene_geom <- data.frame(
  
  gene = genes,
  
  stringsAsFactors = FALSE
  
) %>%
  
  left_join(
    
    gene_group %>%
      mutate(
        gene = as.character(gene),
        group = as.character(group)
      ),
    
    by = "gene"
  ) %>%
  
  mutate(
    x = seq_along(gene)
  )


############################
# 20. One box per marker group
############################

group_geom <- gene_geom %>%
  
  group_by(group) %>%
  
  summarise(
    
    xmin = min(x) - 0.5,
    
    xmax = max(x) + 0.5,
    
    xmid = mean(
      c(
        min(x) - 0.5,
        max(x) + 0.5
      )
    ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    # Use the same CT color for
    # the corresponding marker group
    color = CT_COLORS[group],
    
    
    # White text for darker colors
    txt_col = ifelse(
      
      group %in% c(
        "DG ExNs",
        "CA1 ExNs",
        "CA3 ExNs",
        "Oligodendrocytes",
        "Microglia",
        "Macrophage",
        "Epithelial",
        "T Cell",
        "APC",
        "Stromal"
      ),
      
      "white",
      "black"
    )
  )


############################
# 21. Grouped marker strip
############################

p_strip <- ggplot() +
  
  
  # ONE rectangle for EACH marker group
  #
  # For example:
  #
  # DG ExNs:
  #   Prox1 | Glis3 | Egr3
  #
  # all get ONE box.

  geom_rect(
    
    data = group_geom,
    
    aes(
      
      xmin = xmin,
      
      xmax = xmax,
      
      ymin = 0,
      
      ymax = 1,
      
      fill = color
    ),
    
    color = "white",
    
    linewidth = 0.5
  ) +
  
  
  # One centered group label
  
  geom_text(
    
    data = group_geom,
    
    aes(
      
      x = xmid,
      
      y = 0.5,
      
      label = group,
      
      color = txt_col
    ),
    
    angle = 90,
    
    size = 3,
    
    fontface = "bold"
  ) +
  
  
  scale_fill_identity() +
  
  scale_color_identity() +
  
  
  scale_x_continuous(
    
    breaks = gene_geom$x,
    
    labels = NULL,
    
    expand = c(
      0,
      0
    )
  ) +
  
  
  scale_y_continuous(
    
    limits = c(
      0,
      1
    ),
    
    expand = c(
      0,
      0
    )
  ) +
  
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  
  theme_void() +
  
  
  theme(
    
    plot.margin = margin(
      0,
      10,
      0,
      10
    )
  )


############################
# 22. Gene labels
############################

p_gene_labels <- ggplot(
  
  gene_geom,
  
  aes(
    x = x,
    y = 0,
    label = gene
  )
) +
  
  
  geom_text(
    
    angle = 90,
    
    hjust = 1,
    
    vjust = 0.5,
    
    size = 2.8
  ) +
  
  
  scale_x_continuous(
    
    breaks = gene_geom$x,
    
    labels = NULL,
    
    expand = c(
      0,
      0
    )
  ) +
  
  
  scale_y_continuous(
    
    limits = c(
      -1,
      0
    ),
    
    expand = c(
      0,
      0
    )
  ) +
  
  
  theme_void() +
  
  
  theme(
    
    plot.margin = margin(
      0,
      10,
      0,
      10
    )
  )


############################
# 23. Combine plots
############################

final_plot <-
  
  p_violin /
  
  p_strip /
  
  p_gene_labels +
  
  plot_layout(
    
    heights = c(
      1,
      0.20,
      0.35
    )
  )


############################
# 24. Display
############################

print(final_plot)


############################
# 25. Output directory
############################

out_dir <- file.path(BASE, "Figures", "ED3_ClusterAnnotation")


if (!dir.exists(out_dir)) {
  
  dir.create(
    out_dir,
    recursive = TRUE
  )
}


############################
# 26. Save PDF
############################

ggsave(
  
  filename = file.path(
    out_dir,
    "ClusterAnnotation_StackedViolin.pdf"
  ),
  
  plot = final_plot,
  
  width = 18,
  
  height = 16,
  
  units = "in",
  
  limitsize = FALSE
)


############################
# 27. Save PNG
############################

ggsave(
  
  filename = file.path(
    out_dir,
    "ClusterAnnotation_StackedViolin.png"
  ),
  
  plot = final_plot,
  
  width = 18,
  
  height = 16,
  
  units = "in",
  
  dpi = 300,
  
  limitsize = FALSE
)


############################
# 28. Confirmation
############################

cat(
  
  "\n\n========================================\n",
  
  "Cluster annotation plot complete!\n",
  
  "========================================\n\n",
  
  "Output directory:\n",
  
  out_dir,
  
  "\n\n",
  
  "Files:\n",
  
  "  ClusterAnnotation_StackedViolin.pdf\n",
  
  "  ClusterAnnotation_StackedViolin.png\n\n",
  
  "Marker group order:\n",
  
  paste(
    seq_along(names(marker_genes)),
    names(marker_genes),
    sep = ". ",
    collapse = "\n"
  ),
  
  "\n\n========================================\n"
)
