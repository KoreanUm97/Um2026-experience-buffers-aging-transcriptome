################################################################################
### Fig4b_ED8c_UMAP.R
###
### UMAPs of the DREADD dataset (female only; IN1 + IN2 merged into "IN").
###   Fig. 4b   cell types; mScarlet (raw UMI, cells < 50 UMI grey); age; stimulation (C21 vs PBS)
###   ED8c      mScarlet split by age x stimulation
### mScarlet colour scale: raw UMI from MSCARLET_THRESH (50) to MSCARLET_CAP (600); cells above the cap are drawn in the
### darkest colour, cells below the threshold in grey.
###
### Run:  Rscript Fig4b_ED8c_UMAP.R     (set OBJ_RDS and OUT_DIR below first)
################################################################################

suppressPackageStartupMessages({
  library(conflicted)
  library(Seurat)
  library(BPCells)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  conflicts_prefer(dplyr::select, .quiet = TRUE)
  conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicts_prefer(dplyr::mutate, .quiet = TRUE)
})

################################################################################
### SETTINGS
################################################################################

OBJ_RDS <- "PATH_TO_DREADD_OBJECT.rds"          # DREADD Seurat object (GEO)
OUT_DIR <- file.path("PATH_TO_OUTPUT", "Paper_Figures", "Figure4", "UMAP")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PT_SIZE         <- 0.5
EXPORT_DPI      <- 300
MSCARLET_THRESH <- 50    # cells below this are shown as grey background
MSCARLET_CAP    <- 600   # top of the colour scale (cells above are drawn in the darkest colour)

should_skip <- function(f) {
  if (file.exists(f)) { cat("  Skipping:", basename(f), "\n"); return(TRUE) }
  FALSE
}
save_plot <- function(p, pdf_path, width, height) {
  ggsave(pdf_path, p, width = width, height = height, device = "pdf")
  ggsave(sub("\\.pdf$", ".png", pdf_path), p, width = width, height = height,
         device = "png", dpi = EXPORT_DPI, bg = "white")
  cat("  Saved:", basename(pdf_path), "\n")
}

################################################################################
### LOAD + METADATA
################################################################################

cat("Loading", OBJ_RDS, "\n")
obj <- readRDS(OBJ_RDS)

all_reds <- names(obj@reductions)
umap_key <- if ("full.umap.rpca.sketch" %in% all_reds) "full.umap.rpca.sketch" else
            if ("umap.rpca.sketch"       %in% all_reds) "umap.rpca.sketch"      else
            if ("umap"                   %in% all_reds) "umap"                  else all_reds[1]
cat(sprintf("Using reduction: %s (%d cells)\n\n", umap_key, nrow(obj@reductions[[umap_key]]@cell.embeddings)))

# Merge IN subtypes
in_present <- intersect(c("IN1", "IN2"), unique(obj@meta.data$celltype))
if (length(in_present) > 0) obj@meta.data$celltype[obj@meta.data$celltype %in% in_present] <- "IN"

obj@meta.data <- obj@meta.data %>%
  dplyr::mutate(
    age_norm       = dplyr::case_when(tolower(age) == "young" ~ "Young",
                                      tolower(age) == "old"   ~ "Old",
                                      TRUE ~ tools::toTitleCase(tolower(age))),
    condition_norm = experimental.condition,
    age_condition  = paste(age_norm, condition_norm))

# mScarlet raw UMI; the display column is NA below the threshold (drawn grey).
# Named "mScarlet" so the feature label on the plots reads "mScarlet".
if (!"mScarlet_raw" %in% colnames(obj@meta.data))
  obj@meta.data$mScarlet_raw <- FetchData(obj, vars = "mScarlet", layer = "counts")$mScarlet
obj@meta.data[["mScarlet"]] <- ifelse(obj@meta.data$mScarlet_raw >= MSCARLET_THRESH,
                                      obj@meta.data$mScarlet_raw, NA_real_)

cat("Cells:", ncol(obj), "\n"); print(table(obj$celltype)); print(table(obj$age_condition))
cat(sprintf("mScarlet >= %d UMI: %d cells\n", MSCARLET_THRESH, sum(!is.na(obj[["mScarlet"]]))))

################################################################################
### COLOURS + THEME
################################################################################

CT_PALETTE <- c(
  "DG ExNs" = "#E53935", "CA1 ExNs" = "#CC0000", "CA2 ExNs" = "#8E24AA", "CA3 ExNs" = "#5E35B1",
  "Sub ExNs" = "#80D8FF", "IN" = "#00897B", "Astrocytes" = "#43A047", "Oligodendrocytes" = "#1E88E5",
  "OPC" = "#00ACC1", "Microglia" = "#7CB342", "Endothelial" = "#FB8C00", "Pericyte" = "#6D4C41",
  "APC" = "#546E7A", "Stromal" = "#27408B", "Macrophage" = "#C62828", "Cajal-Retzius Cell" = "#FDD835",
  "T Cell" = "#757575", "Ambiguous" = "#BDBDBD")
celltypes <- sort(unique(obj$celltype))
CT_PALETTE[setdiff(celltypes, names(CT_PALETTE))] <- "#999999"
ct_palette <- CT_PALETTE[celltypes]
age_cols  <- c(Young = "#1565C0", Old = "#FF1744")
cond_cols <- c(PBS = "#1565C0", C21 = "#FF1744")

MSCARLET_SCALE <- scale_color_gradientn(
  colors   = c("#FFF5F0", "#FCBBA1", "#FC6B50", "#EF3B2C", "#CB181D"),
  limits   = c(MSCARLET_THRESH, MSCARLET_CAP),
  oob      = scales::squish,
  na.value = "grey88",
  name     = sprintf("mScarlet\n(raw UMI\n≥%d)", MSCARLET_THRESH))

umap_theme <- theme_classic(base_size = 11) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_text(size = 9, color = "grey50"),
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8.5, color = "grey40", hjust = 0.5),
        legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))

################################################################################
### Fig. 4b: CELL TYPES, AGE, STIMULATION
################################################################################

out_f <- file.path(OUT_DIR, "UMAP_CellType.pdf")
if (!should_skip(out_f)) {
  p <- DimPlot(obj, reduction = umap_key, group.by = "celltype", label = TRUE, label.size = 3.5, repel = TRUE,
               pt.size = PT_SIZE, raster = FALSE) +
    scale_color_manual(values = ct_palette) + umap_theme +
    labs(title = "DREADD dataset: cell types", subtitle = paste0("n=", ncol(obj), " cells | Female only | IN1+IN2 merged"),
         x = "UMAP 1", y = "UMAP 2", color = "Cell type")
  save_plot(p, out_f, 10, 8)
}

out_f <- file.path(OUT_DIR, "UMAP_CellType_NoLabel.pdf")
if (!should_skip(out_f)) {
  p <- DimPlot(obj, reduction = umap_key, group.by = "celltype", label = FALSE, pt.size = PT_SIZE, raster = FALSE) +
    scale_color_manual(values = ct_palette) + umap_theme +
    labs(title = "DREADD dataset: cell types", subtitle = "Female only | IN1+IN2 merged",
         x = "UMAP 1", y = "UMAP 2", color = "Cell type")
  save_plot(p, out_f, 10, 8)
}

out_f <- file.path(OUT_DIR, "UMAP_Age.pdf")
if (!should_skip(out_f)) {
  p <- DimPlot(obj, reduction = umap_key, group.by = "age_norm", pt.size = PT_SIZE, shuffle = TRUE, raster = FALSE) +
    scale_color_manual(values = age_cols) + umap_theme +
    labs(title = "DREADD dataset: age", subtitle = "Female only", x = "UMAP 1", y = "UMAP 2", color = "Age")
  save_plot(p, out_f, 8, 7)
}

out_f <- file.path(OUT_DIR, "UMAP_Condition.pdf")
if (!should_skip(out_f)) {
  p <- DimPlot(obj, reduction = umap_key, group.by = "condition_norm", pt.size = PT_SIZE, shuffle = TRUE, raster = FALSE) +
    scale_color_manual(values = cond_cols) + umap_theme +
    labs(title = "DREADD dataset: condition", subtitle = "C21 = DREADD activated | PBS = vehicle | Female only",
         x = "UMAP 1", y = "UMAP 2", color = "Condition")
  save_plot(p, out_f, 8, 7)
}

################################################################################
### Fig. 4b: mScarlet;  ED8c: mScarlet split by age x stimulation
################################################################################

out_f <- file.path(OUT_DIR, "UMAP_mScarlet.pdf")
if (!should_skip(out_f)) {
  p <- FeaturePlot(obj, features = "mScarlet", reduction = umap_key, pt.size = PT_SIZE,
                   order = TRUE, raster = FALSE) +                 # expressing cells drawn on top
    MSCARLET_SCALE + umap_theme +
    labs(title = sprintf("DREADD dataset: mScarlet (threshold=%d UMI)", MSCARLET_THRESH),
         subtitle = "Grey = below threshold (ambient RNA) | Colored = expressing cells | Female only",
         x = "UMAP 1", y = "UMAP 2")
  save_plot(p, out_f, 8, 7)
}

out_f <- file.path(OUT_DIR, "UMAP_mScarlet_ByAgeCondition.pdf")
if (!should_skip(out_f)) {
  p <- FeaturePlot(obj, features = "mScarlet", reduction = umap_key, split.by = "age_condition", ncol = 2,
                   pt.size = PT_SIZE, order = TRUE, raster = FALSE) &
    MSCARLET_SCALE & umap_theme & labs(x = "UMAP 1", y = "UMAP 2")
  save_plot(p, out_f, 14, 14)
}

cat("\n=== DONE ===\nOutputs:", OUT_DIR, "\n")
