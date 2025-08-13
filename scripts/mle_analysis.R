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
  library(dplyr)
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
oldwd <- getwd()
setwd(output_dir)

FluteMLE(
  replicates,
  treatname = treatname,
  ctrlname = ctrlname,
  proj = proj_name,
  organism = organism,
  norm_method = norm_method
)

setwd(oldwd)
# -------------------------
# Generate Beta Score Histograms
# -------------------------
# pull out beta scores
gdata <- ReadBeta(replicates)

# sanity check: columns exist
stopifnot(ctrlname %in% names(gdata), treatname %in% names(gdata))

# Simple QC: histogram of beta scores from the MLE gene summary.
p1 <- ggplot(gdata, aes(x = .data[[ctrlname]])) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  labs(title = sprintf("Beta score distribution (%s vs %s)", treatname, ctrlname),
       x = "Beta", y = "Count")

p2 <- ggplot(gdata, aes(x = .data[[treatname]])) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  labs(title = sprintf("Beta score distribution (%s vs %s)", treatname, ctrlname),
       x = "Beta", y = "Count")

# plot histograms them side by side
out_hist_png <- file.path(output_dir, "beta_hist.png")
hists <- cowplot::plot_grid(p1, p2, nrow = 1)
ggsave(out_hist_png, plot = hists, width = 12, height = 6, dpi = 150)

# -------------------------
# Plot normalization (none | cell_cycle | loess)
# -------------------------
# methods for normalization
gdata_cc = NormalizeBeta(gdata, samples=c(ctrlname, treatname), method="cell_cycle")
gdata_loe = NormalizeBeta(gdata, samples=c(ctrlname, treatname), method="loess")

#to compare density of the non-normalized betas
p_raw <- DensityView(gdata, samples=c(ctrlname, treatname))

#to compare density of the normalized betas
p_cc <- DensityView(gdata_cc, samples=c(ctrlname, treatname))

#to compare density of the normalized betas
p_loe <- DensityView(gdata_loe, samples=c(ctrlname, treatname))

# plot them side by side
norm_grid <- cowplot::plot_grid(p_raw, p_cc, p_loe, labels = c("NonNormalized", "CellCycle", "Loess"), nrow = 1)
out_norm_png <- file.path(output_dir, "beta_norm.png")
ggsave(out_norm_png, plot = norm_grid, width = 12, height = 6, dpi = 150)

# -------------------------
# Choose normalization from config (none | cell_cycle | loess)
# -------------------------
nm <- tolower(norm_method)
if (nm == "cell_cycle") {
  gsel <- NormalizeBeta(gdata, samples = c(ctrlname, treatname), method = "cell_cycle")
} else if (nm == "loess") {
  gsel <- NormalizeBeta(gdata, samples = c(ctrlname, treatname), method = "loess")
} else {
  gsel <- gdata  # "none" or anything else falls back to raw
}

# -------------------------
# Build positive/negative selection table
# -------------------------
# Dynamic FDR column names in the gene_summary:
fdr_ctrl_col  <- paste0(ctrlname,  ".fdr")
fdr_treat_col <- paste0(treatname, ".fdr")

# Start with betas (possibly normalized) and diff (treat - ctrl)
sel_df <- data.frame(
  Gene = gsel$Gene,
  !!ctrlname := as.numeric(gsel[[ctrlname]]),    # control beta
  !!treatname := as.numeric(gsel[[treatname]]),   # treatment beta
  diff = as.numeric(gsel[[treatname]]) - as.numeric(gsel[[ctrlname]])
)

# Add FDRs if present in the MLE gene summary (replicates). Align by Gene to be safe.
have_ctrl_fdr  <- fdr_ctrl_col  %in% names(replicates)
have_treat_fdr <- fdr_treat_col %in% names(replicates)

fdr_keep <- c("Gene",
              if (have_ctrl_fdr)  fdr_ctrl_col,
              if (have_treat_fdr) fdr_treat_col)

if (length(fdr_keep) > 1) {
  fdr_df <- replicates[, fdr_keep, drop = FALSE]
  sel_df <- merge(sel_df, fdr_df, by = "Gene", all.x = TRUE)
} else {
  # If no FDR cols, add NAs with informative names
  if (!have_ctrl_fdr)  sel_df[[paste0(ctrlname,  "_FDR")]]  <- NA_real_
  if (!have_treat_fdr) sel_df[[paste0(treatname, "_FDR")]] <- NA_real_
}

# If original FDR columns exist with ".fdr" suffix, also add friendly aliases
if (fdr_ctrl_col %in% names(sel_df))  sel_df[[paste0(ctrlname,  "_FDR")]]  <- sel_df[[fdr_ctrl_col]]
if (fdr_treat_col %in% names(sel_df)) sel_df[[paste0(treatname, "_FDR")]] <- sel_df[[fdr_treat_col]]

# -------------------------
# Beta vs Beta scatter (color = treatment FDR) with top-50 labels
# -------------------------
treat_fdr_col <- paste0(treatname, ".fdr")

df_scatter <- merge(
  gsel[, c("Gene", ctrlname, treatname)],
  replicates[, c("Gene", treat_fdr_col)],
  by = "Gene",
  all.x = TRUE
)

p_beta_scatter <- ggplot(
  df_scatter,
  aes(
    x = .data[[ctrlname]],
    y = .data[[treatname]],
    color = .data[[treat_fdr_col]]
  )
) +
  geom_point(alpha = 0.8, size = 1.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = sprintf("Beta vs Beta (%s): %s vs %s (color = %s FDR)",
                    nm, treatname, ctrlname, treatname),
    x = sprintf("%s beta", ctrlname),
    y = sprintf("%s beta", treatname),
    color = sprintf("%s FDR", treatname)
  ) +
  scale_colour_gradient(low = "#ff7f00", high = "#7570b3")

# Top 50 by lowest treatment FDR (ignore NAs)
top_lab <- df_scatter |>
  dplyr::filter(!is.na(.data[[treat_fdr_col]])) |>
  dplyr::arrange(.data[[treat_fdr_col]]) |>
  dplyr::slice_head(n = 50)

p_beta_scatter_labeled <- p_beta_scatter +
  ggrepel::geom_text_repel(
    data = top_lab,
    aes(label = Gene),
    size = 2.5, max.overlaps = 100,
    box.padding = 0.3, point.padding = 0.2, segment.size = 0.2
  )

ggsave(file.path(output_dir, "beta_scatter.png"),
       plot = p_beta_scatter, width = 8, height = 8, dpi = 150)

ggsave(file.path(output_dir, "beta_scatter_labeled.png"),
       plot = p_beta_scatter_labeled, width = 8, height = 8, dpi = 150)
              
