#!/usr/bin/env Rscript

# -------------------------
# Load Libraries
# -------------------------
suppressPackageStartupMessages({
  library(MAGeCKFlute)
  library(cowplot)
  library(ggplot2)
  library(ggrepel)
})

# -------------------------
# Capture command-line arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
gene_summary <- strsplit(args[1], " ")[[1]]
count_summary <- strsplit(args[2], " ")[[1]]
design_file <- args[3]
output_dir <- args[4]

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

# Project name defaults to the directory name
proj <- basename(normalizePath(output_dir, mustWork = FALSE))

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
  proj = proj,
  organism = "hsa",
  norm_method = "none"
)

# -------------------------
# Generate Custom Plots
# -------------------------
# Simple QC: histogram of beta scores from the MLE gene summary.
p <- ggplot(replicates, aes(x = beta)) +
  geom_histogram(bins = 60) +
  theme_bw() +
  labs(title = sprintf("Beta score distribution (%s vs %s)", treatname, ctrlname),
       x = "Beta", y = "Count")

ggsave(out_png, plot = p, width = 8, height = 6, dpi = 150)
