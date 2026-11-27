################################################################################
### run_demo.R  -  small simulated demo of the core analysis steps
###
### Simulates pseudobulk counts (2 cell types x 3 ages x 2 conditions x 3 mice, one sex)
### with (i) activity-regulated genes (ARGs) whose induction shrinks with age and
### (ii) aging genes that stimulation partially reverses, then runs the same core steps
### as the full pipelines:
###   - pseudobulk DESeq2, design ~ grp, genes with >= 10 counts in >= 2 samples
###     (01_DESeq2/NE_DESeq2_Pipeline.R)
###   - ARG calls (padj < 0.05, |log2FC| >= 0.25) and dynamic range of the young ARG
###     program at each age (02_Dynamic_Range_Compression)
###   - buffering: cosine between the aging (OHCYHC) and stimulation (OEEOHC) fold-change
###     vectors, and buffering index BI = sum(OEEYHC^2) / sum(OHCYHC^2) (04_Buffering)
###
### Usage:  Rscript demo/run_demo.R          (from the repository root)
### Output: demo/demo_output/  (DESeq2 tables, demo_summary.csv, demo_plots.pdf)
### Runtime: about 1 minute on a standard desktop.
################################################################################

suppressPackageStartupMessages({ library(DESeq2) })
set.seed(1)
t0 <- Sys.time()

out_dir <- file.path("demo", "demo_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. simulate pseudobulk counts ----------------------------------------
n_genes <- 3000
genes   <- sprintf("Gene%04d", seq_len(n_genes))
ages    <- c("Y", "M", "O")
conds   <- c("HC", "EE")
n_mice  <- 3
cts     <- c("CA1_ExNs", "Astrocytes")

arg_idx   <- 1:150                       # activity-regulated genes (induced by EE)
aging_idx <- 151:450                     # genes changed with age
arg_lfc   <- c(Y = 1.5, M = 1.0, O = 0.5)                   # induction contracts with age
aging_lfc <- sample(c(-0.8, 0.8), length(aging_idx), TRUE)  # aging direction per gene
age_frac  <- c(Y = 0, M = 0.5, O = 1)                       # fraction of the aging change present
buffer    <- 0.4                                            # EE reverses 40% of the aging change

meta <- expand.grid(mouse = seq_len(n_mice), cond = conds, age = ages, stringsAsFactors = FALSE)
meta$grp    <- paste0(meta$age, meta$cond)
meta$sample <- paste0(meta$grp, "_m", meta$mouse)

simulate_ct <- function() {
  base <- rgamma(n_genes, shape = 2, rate = 2 / 200)          # mean expression
  sapply(seq_len(nrow(meta)), function(j) {
    a <- meta$age[j]; ee <- meta$cond[j] == "EE"
    lfc <- numeric(n_genes)
    lfc[aging_idx] <- aging_lfc * age_frac[a] * (if (ee) 1 - buffer else 1)
    if (ee) lfc[arg_idx] <- lfc[arg_idx] + arg_lfc[a]
    mu <- base * 2^lfc * exp(rnorm(1, 0, 0.1))               # per-mouse size factor
    rnbinom(n_genes, mu = mu, size = 1 / 0.05)
  })
}

## ---- 2. pseudobulk DESeq2 per cell type ------------------------------------
comps <- list(YEEYHC = c("YEE", "YHC"), MEEMHC = c("MEE", "MHC"), OEEOHC = c("OEE", "OHC"),
              OHCYHC = c("OHC", "YHC"), OEEYHC = c("OEE", "YHC"))
res <- list()
for (ct in cts) {
  cnt <- simulate_ct(); dimnames(cnt) <- list(genes, meta$sample)
  keep <- rowSums(cnt >= 10) >= 2
  dds <- DESeqDataSetFromMatrix(cnt[keep, ], data.frame(grp = factor(meta$grp), row.names = meta$sample), ~ grp)
  dds <- DESeq(dds, quiet = TRUE)
  for (cp in names(comps)) {
    r <- as.data.frame(results(dds, contrast = c("grp", comps[[cp]])))
    r$gene <- rownames(r)
    write.csv(r, file.path(out_dir, sprintf("DESeq2_%s_%s.csv", ct, cp)), row.names = FALSE)
    res[[ct]][[cp]] <- r
  }
}

## ---- 3. ARGs, dynamic range, buffering -------------------------------------
sig_up <- function(r) r$gene[!is.na(r$padj) & r$padj < 0.05 & r$log2FoldChange >= 0.25]
cosine <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))

summ <- do.call(rbind, lapply(cts, function(ct) {
  R <- res[[ct]]
  young_args <- sig_up(R$YEEYHC)
  med_lfc <- function(cp) median(R[[cp]]$log2FoldChange[match(young_args, R[[cp]]$gene)], na.rm = TRUE)
  aging_g <- R$OHCYHC$gene[!is.na(R$OHCYHC$padj) & R$OHCYHC$padj < 0.05]
  lfc <- function(cp) R[[cp]]$log2FoldChange[match(aging_g, R[[cp]]$gene)]
  data.frame(celltype = ct,
             n_ARG_young = length(young_args),
             n_ARG_middle = length(sig_up(R$MEEMHC)),
             n_ARG_old = length(sig_up(R$OEEOHC)),
             medLFC_youngARGs_young = round(med_lfc("YEEYHC"), 2),
             medLFC_youngARGs_middle = round(med_lfc("MEEMHC"), 2),
             medLFC_youngARGs_old = round(med_lfc("OEEOHC"), 2),
             n_aging_genes = length(aging_g),
             cosine_aging_vs_EE = round(cosine(lfc("OHCYHC"), lfc("OEEOHC")), 2),
             buffering_index = round(sum(lfc("OEEYHC")^2) / sum(lfc("OHCYHC")^2), 2))
}))
write.csv(summ, file.path(out_dir, "demo_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)

## ---- 4. plots ----------------------------------------------------------------
pdf(file.path(out_dir, "demo_plots.pdf"), width = 9, height = 4.5)
par(mfrow = c(1, 2))
m <- t(as.matrix(summ[, c("medLFC_youngARGs_young", "medLFC_youngARGs_middle", "medLFC_youngARGs_old")]))
barplot(m, beside = TRUE, names.arg = summ$celltype, legend.text = c("young", "middle", "old"),
        ylab = "median log2FC of young ARGs (EE vs HC)", main = "Dynamic range compression")
R <- res[[cts[1]]]; g <- R$OHCYHC$gene[!is.na(R$OHCYHC$padj) & R$OHCYHC$padj < 0.05]
plot(R$OHCYHC$log2FoldChange[match(g, R$OHCYHC$gene)], R$OEEOHC$log2FoldChange[match(g, R$OEEOHC$gene)],
     pch = 16, cex = 0.6, xlab = "aging log2FC (OHC vs YHC)", ylab = "stimulation log2FC (OEE vs OHC)",
     main = paste0("Buffering, ", cts[1])); abline(h = 0, v = 0, lty = 2)
invisible(dev.off())

cat(sprintf("\nDone in %.1f s. Output in %s\n", as.numeric(difftime(Sys.time(), t0, units = "secs")), out_dir))
