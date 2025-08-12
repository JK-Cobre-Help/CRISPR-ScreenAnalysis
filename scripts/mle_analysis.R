#!/usr/bin/env Rscript

# -------------------------
# Load Libraries
# -------------------------
suppressPackageStartupMessages({
  library(MAGeCKFlute)
  library(cowplot)
  library(ggplot2)
  library(ggrepel)
  library(rlang)
})

# -------------------------
# Capture command-line arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
gene_summary <- args[1]
count_summary <- args[2]
design_file <- args[3]
output_dir <- args[4]
proj_name <- args[5]
organism <- args[6]
norm_method <- args[7]

# -------------------------
# Ensure output directory exists
# -------------------------
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -------------------------
# Load Gene Summary table and Count Summary
# -------------------------
replicates <- read.table(gene_summary, header = TRUE, sep = "\t", check.names = FALSE)
countsummary <- read.table(count_summary, header = TRUE, sep = "\t", check.names = FALSE)
design <- read.table(design_file, header = TRUE, sep = "\t", check.names = FALSE)

# -------------------------
# Derive treat/ctrl from design matrix
# (exclude 'group' and 'baseline'; use remaining factor columns in order)
# -------------------------
cols <- setdiff(colnames(design), c("group","baseline"))
if (length(cols) < 2) {
  stop("Design matrix must contain at least two factor columns (excluding 'baseline').")
}
# Convention: first factor column = control, last factor column = treatment
ctrlname <- cols[1]
treatname <- cols[length(cols)]

# -------------------------
# Run MAGeCK Flute MLE
# -------------------------
# setwd is used to ensure artifacts are written inside output_dir
oldwd <- getwd(); setwd(output_dir)
on.exit(setwd(oldwd), add = TRUE)

FluteMLE(
  replicates,
  treatname = treatname,
  ctrlname = ctrlname,
  proj = proj_name,
  organism = organism,
  norm_method = norm_method
)

# -------------------------
# Generate Custom Plots
# -------------------------

# pull out beta scores
gdata <- ReadBeta(replicates)

# sanity check: columns exist
stopifnot(ctrlname %in% names(gdata), treatname %in% names(gdata))

# Simple QC: histogram of beta scores from the MLE gene summary.
p1 <- ggplot(gdata, aes(x = .data[[ctrlname]])) +
  geom_histogram(bins = 50) +
  theme_bw() +
  labs(title = sprintf("Beta score distribution (%s vs %s)", treatname, ctrlname),
       x = "Beta", y = "Count")

p2 <- ggplot(gdata, aes(x = .data[[treatname]])) +
  geom_histogram(bins = 50) +
  theme_bw() +
  labs(title = sprintf("Beta score distribution (%s vs %s)", treatname, ctrlname),
       x = "Beta", y = "Count")

# plot histograms them side by side
out_hist_png <- file.path(output_dir, "beta_hist.png")
hists <- cowplot::plot_grid(p1, p2, labels = c("Control Beta","Treatment Beta"), nrow = 1)
ggsave(out_hist_png, plot = hists, width = 12, height = 6, dpi = 150)

################################################################################

# methods for normalization
gdata_cc = NormalizeBeta(gdata, samples=c(ctrlname, treatname), method="cell_cycle")
gdata_loe = NormalizeBeta(gdata, samples=c(ctrlname, treatname), method="loess")

#to compare density of the non-normalized betas
p_raw <- DensityView(gdata, samples=c(ctrlname, treatname))

#to compare density of the normalized betas
p_cc <- DensityView(gdata_cc, samples=c(ctrlname, treatname))

#to compare density of the normalized betas
p_loess <- DensityView(gdata_loe, samples=c(ctrlname, treatname))

# plot them side by side
norm_grid <- cowplot::plot_grid(p_raw, p_cc, p_loess, labels = c("NonNormalized", "CellCycle", "Loess"), nrow = 1)
out_norm_png <- file.path(output_dir, "beta_norm.png")
ggsave(out_norm_png, plot = norm_grid, width = 12, height = 6, dpi = 150)

