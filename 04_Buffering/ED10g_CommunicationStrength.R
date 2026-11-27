
################################################################################
### ED10g_CommunicationStrength.R
###
### Total communication strength (CellChat, NeuronChat) per sex x age x condition (6 conditions x 2 sexes),
### built from the integrated NE Seurat object with these settings:
###   CellChat   : all cell types; cell types with < 10 cells dropped; CellChatDB.mouse; triMean;
###                filterCommunication(min.cells = 10); computeCommunProbPathway; aggregateNet
###   NeuronChat : neurons only (DG, CA3, CA1, CA2, Sub ExNs, IN); >= 50 neurons per condition; cell types with
###                < 10 cells dropped; DB = "mouse"; run_NeuronChat(M = 100)
### Total strength = sum of the CellChat weight matrix, and sum of net_aggregation(net, "weight") for NeuronChat.
###
### Objects are cached as <date>_CellChat_EE/<sex>/<date>_<sex>_cc_<cond>.rds (and _nc_ for NeuronChat), so a
### rerun skips finished objects.
################################################################################

suppressPackageStartupMessages({
  library(CellChat); library(NeuronChat); library(Seurat)
  library(dplyr); library(tidyr); library(ggplot2); library(future)
  library(conflicted)
  conflicts_prefer(dplyr::select, .quiet = TRUE); conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicts_prefer(base::intersect, .quiet = TRUE); conflicts_prefer(base::union, .quiet = TRUE)
})
future::plan("sequential")

################################################################################
### SETTINGS
################################################################################
base_path <- "PATH_TO_NE_OUTPUT/"
NE_OBJECT <- "PATH_TO_NE_OBJECT.rds"            # integrated NE Seurat object
fig_date <- format(Sys.Date(), "%Y-%m-%d")      # set to the first run's date when resuming a run on another day
cc_dir   <- paste0(base_path, "CellChat/",   fig_date, "_CellChat_EE/")
nc_dir   <- paste0(base_path, "NeuronChat/", fig_date, "_NeuronChat_EE/")
out_dir  <- paste0(base_path, "CellChat_NeuronChat/", fig_date, "_TotalStrength/")
for (d in c(cc_dir, nc_dir, out_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

AGE_KEY <- "age"; COND_KEY <- "experimental.condition"; SEX_KEY <- "sex"
SEXES <- c("female", "male")
CONDS <- list(YHC = list(age = "young", cond = "HC"),     YEE = list(age = "young", cond = "EE"),
              MHC = list(age = "Middleage", cond = "HC"), MEE = list(age = "Middleage", cond = "EE"),
              OHC = list(age = "old", cond = "HC"),       OEE = list(age = "old", cond = "EE"))
CT_ORDER <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN","Astrocytes","Oligodendrocytes","OPC","Microglia",
              "Endothelial","Pericyte","APC","Stromal","Cajal-Retzius Cell","Macrophage","Epithelial","T Cell")
NEURON_CTS <- c("DG ExNs","CA3 ExNs","CA1 ExNs","CA2 ExNs","Sub ExNs","IN")
CC_MIN_CELLS <- 10
NC_MIN_CELLS <- 50
NC_M_PERMS   <- 100

################################################################################
### HELPERS
################################################################################
to_dgc <- function(expr) if (inherits(expr, "dgCMatrix")) expr else
  tryCatch(as(expr, "dgCMatrix"), error = function(e) as(BPCells::write_matrix_memory(expr), "dgCMatrix"))

make_cc <- function(sub) {
  expr <- to_dgc(GetAssayData(sub, assay = "RNA", layer = "data"))
  ct_labels <- sub$celltype; names(ct_labels) <- colnames(sub)
  ct_counts <- table(ct_labels); keep_cts <- names(ct_counts[ct_counts >= CC_MIN_CELLS])
  if (length(keep_cts) < 3) { cat("    too few cell types - skip\n"); return(NULL) }
  keep_cells <- names(ct_labels)[ct_labels %in% keep_cts]
  expr <- expr[, keep_cells]; ct_labels <- ct_labels[keep_cells]
  ct_labels <- factor(ct_labels, levels = intersect(CT_ORDER, unique(ct_labels)))
  tryCatch({
    cc <- createCellChat(object = expr, meta = data.frame(labels = ct_labels, row.names = names(ct_labels)), group.by = "labels")
    cc@DB <- CellChatDB.mouse
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, type = "triMean")
    cc <- filterCommunication(cc, min.cells = CC_MIN_CELLS)
    cc <- computeCommunProbPathway(cc)
    aggregateNet(cc)
  }, error = function(e) { cat("    CellChat error:", e$message, "\n"); NULL })
}

make_nc <- function(sub) {
  if (is.null(sub) || ncol(sub) < NC_MIN_CELLS) return(NULL)
  nc_cells <- colnames(sub)[sub$celltype %in% NEURON_CTS]
  if (length(nc_cells) < NC_MIN_CELLS) return(NULL)
  sub <- sub[, nc_cells]
  ct_counts <- table(sub$celltype); keep_cts <- names(ct_counts[ct_counts >= 10])
  if (length(keep_cts) < 2) return(NULL)
  sub <- sub[, sub$celltype %in% keep_cts]
  expr <- to_dgc(GetAssayData(sub, assay = "RNA", layer = "data"))
  cell_groups <- sub$celltype; names(cell_groups) <- colnames(sub)
  tryCatch({
    nc <- createNeuronChat(object = expr, DB = "mouse", group.by = cell_groups)
    run_NeuronChat(nc, M = NC_M_PERMS)
  }, error = function(e) { cat("    NeuronChat error:", e$message, "\n"); NULL })
}

################################################################################
### BUILD (or load) ALL 12 + 12 OBJECTS
################################################################################
cat("Reading", NE_OBJECT, "\n")
df_harmony <- readRDS(NE_OBJECT)
cat("Sex x age x condition (cells):\n")
print(table(df_harmony@meta.data[[SEX_KEY]], paste(df_harmony@meta.data[[AGE_KEY]], df_harmony@meta.data[[COND_KEY]])))

cc_all <- list(); nc_all <- list()
for (sx in SEXES) {
  cat(sprintf("\n--- %s ---\n", toupper(sx)))
  for (d in c(cc_dir, nc_dir)) dir.create(file.path(d, sx), recursive = TRUE, showWarnings = FALSE)
  cc_all[[sx]] <- list(); nc_all[[sx]] <- list()
  sx_obj <- df_harmony[, colnames(df_harmony)[tolower(df_harmony@meta.data[[SEX_KEY]]) == sx]]
  for (nm in names(CONDS)) {
    f_cc <- file.path(cc_dir, sx, paste0(fig_date, "_", sx, "_cc_", nm, ".rds"))
    f_nc <- file.path(nc_dir, sx, paste0(fig_date, "_", sx, "_nc_", nm, ".rds"))
    cells <- colnames(sx_obj)[sx_obj@meta.data[[AGE_KEY]] == CONDS[[nm]]$age & sx_obj@meta.data[[COND_KEY]] == CONDS[[nm]]$cond]
    cat(sprintf("  %s: %d cells\n", nm, length(cells)))
    sub <- if (length(cells)) sx_obj[, cells] else NULL
    if (file.exists(f_cc)) { cc_all[[sx]][[nm]] <- readRDS(f_cc); cat("    CellChat: cached\n") }
    else if (length(cells) >= CC_MIN_CELLS * 3) {
      cc <- make_cc(sub); if (!is.null(cc)) { saveRDS(cc, f_cc); cc_all[[sx]][[nm]] <- cc; cat("    CellChat: built\n") } }
    if (file.exists(f_nc)) { nc_all[[sx]][[nm]] <- readRDS(f_nc); cat("    NeuronChat: cached\n") }
    else { nc <- make_nc(sub); if (!is.null(nc)) { saveRDS(nc, f_nc); nc_all[[sx]][[nm]] <- nc; cat("    NeuronChat: built\n") } }
  }
  rm(sx_obj); invisible(gc())
}
rm(df_harmony); invisible(gc())

################################################################################
### TOTAL STRENGTH + PLOTS
################################################################################
tot_cc <- function(cc) if (is.null(cc)) NA_real_ else sum(cc@net$weight, na.rm = TRUE)
tot_nc <- function(nc) if (is.null(nc)) NA_real_ else
  tryCatch(sum(net_aggregation(nc@net, method = "weight"), na.rm = TRUE), error = function(e) NA_real_)
all_str_df <- dplyr::bind_rows(lapply(SEXES, function(sx) dplyr::bind_rows(lapply(names(CONDS), function(nm) data.frame(
  sex = sx, condition = nm, cc_total = tot_cc(cc_all[[sx]][[nm]]), nc_total = tot_nc(nc_all[[sx]][[nm]]),
  age = dplyr::case_when(grepl("^Y", nm) ~ "Young", grepl("^M", nm) ~ "Middle", TRUE ~ "Old"),
  cond = ifelse(grepl("EE$", nm), "EE", "HC")))))) %>%
  dplyr::mutate(age = factor(age, levels = c("Young", "Middle", "Old")), sex_cond = paste0(sex, "_", cond))
write.csv(all_str_df %>% dplyr::select(sex, condition, age, cond, cc_total, nc_total),
          file.path(out_dir, paste0(fig_date, "_TotalStrength.csv")), row.names = FALSE)
print(all_str_df %>% dplyr::select(sex, condition, cc_total, nc_total), row.names = FALSE)

save_both <- function(p, stem, w, h) {
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h)
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 300, bg = "white"); cat("  Saved:", basename(stem), "\n")
}
line_cols <- c(EE = "#FF1744", HC = "#1565C0")
for (metric in c("cc_total", "nc_total")) {
  lbl <- ifelse(metric == "cc_total", "CellChat (all cell types)", "NeuronChat (neurons)")
  d <- all_str_df %>% dplyr::filter(!is.na(.data[[metric]]))
  # ED10g: both sexes drawn on one panel (solid = female, dashed = male)
  p_ov <- ggplot(d, aes(x = age, y = .data[[metric]], color = cond, linetype = sex, group = sex_cond)) +
    geom_line(linewidth = 1.2) + geom_point(aes(shape = sex), size = 3) +
    scale_color_manual(values = line_cols, name = "Condition") +
    scale_linetype_manual(values = c(female = "solid", male = "dashed"), name = "Sex") +
    scale_shape_manual(values = c(female = 16, male = 17), name = "Sex") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.subtitle = element_text(size = 8.5, color = "grey40", hjust = 0.5), legend.position = "top") +
    labs(title = paste0("Total communication strength\n", lbl),
         subtitle = "Solid=Female | Dashed=Male | Red=EE | Blue=HC", x = "Age", y = "Total communication strength")
  save_both(p_ov, file.path(out_dir, paste0(fig_date, "_TotalStrength_", metric)), 7, 6)
  # one panel per sex (sexes not overlaid)
  p_sx <- ggplot(d, aes(x = age, y = .data[[metric]], color = cond, group = cond)) +
    geom_line(linewidth = 1.2) + geom_point(size = 3) +
    scale_color_manual(values = line_cols, name = NULL) +
    facet_wrap(~sex, nrow = 1, scales = "free_y", labeller = labeller(sex = tools::toTitleCase)) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), strip.text = element_text(face = "bold"),
          strip.background = element_blank(), legend.position = "top") +
    labs(title = paste0("Total communication strength\n", lbl), x = "Age", y = "Total communication strength")
  save_both(p_sx, file.path(out_dir, paste0(fig_date, "_TotalStrength_BySex_", metric)), 8, 4.5)
}
cat("\nDone:", out_dir, "\n")
