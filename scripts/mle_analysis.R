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
  library(tibble)
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
fdr_threshold <- as.numeric(args[8])
y_cap <- as.numeric(args[9])
effect_thr <- as.numeric(args[10])

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
# Figure out ctrl/treat from design
# -------------------------
cols <- setdiff(colnames(design), c("group","baseline"))
if (length(cols) < 2) stop("Design matrix must contain at least two factor columns (excluding 'baseline').")
ctrlname  <- cols[1]
treatname <- cols[length(cols)]

# -------------------------
# Helpers for tolerant column matching (.|_|-||)
# -------------------------
escape_regex <- function(x) gsub("([][{}()+*^$.|\\?\\-])", "\\\\\\1", x)
guess_metric_col <- function(prefix, metric, cols) {
  pat <- paste0("^", escape_regex(prefix), "([._\\-|])?", escape_regex(metric), "$")
  cand <- grep(pat, cols, ignore.case = TRUE, value = TRUE)
  if (length(cand)) return(cand[1])
  tail_pat <- paste0("([._\\-|])", escape_regex(metric), "$")
  cand2 <- grep(tail_pat, cols, ignore.case = TRUE, value = TRUE)
  if (length(cand2)) return(cand2[1])
  NA_character_
}

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
  theme_bw() +
  labs(title = sprintf("Beta score distribution (%s)", ctrlname),
       x = "Beta", y = "Count")

p2 <- ggplot(gdata, aes(x = .data[[treatname]])) +
  geom_histogram(bins = 50) +
  theme_bw() +
  labs(title = sprintf("Beta score distribution (%s)", treatname),
       x = "Beta", y = "Count")

# plot histograms them side by side
out_hist_png <- file.path(output_dir, "beta_hist.png")
hists <- cowplot::plot_grid(p1, p2, nrow = 1)
ggsave(out_hist_png, plot = hists, width = 12, height = 6, dpi = 150)

# -------------------------
# Plot normalization (none | cell_cycle | loess)
# -------------------------
# compare methods for normalization
gdata_cc = NormalizeBeta(gdata, samples=c(ctrlname, treatname), method="cell_cycle")
gdata_loe = NormalizeBeta(gdata, samples=c(ctrlname, treatname), method="loess")
p_raw <- DensityView(gdata, samples=c(ctrlname, treatname))
p_cc <- DensityView(gdata_cc, samples=c(ctrlname, treatname))
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
fdr_ctrl_col  <- guess_metric_col(ctrlname,  "fdr", names(replicates))
fdr_treat_col <- guess_metric_col(treatname, "fdr", names(replicates))

# sanity check the matched columns; stop early if we couldn't find them
if (any(is.na(c(fdr_ctrl_col, fdr_treat_col)))) {
  stop(sprintf(
    "Could not locate FDR columns. Control match: %s | Treatment match: %s\nAvailable columns:\n  - %s",
    ifelse(is.na(fdr_ctrl_col),  "NA", fdr_ctrl_col),
    ifelse(is.na(fdr_treat_col), "NA", fdr_treat_col),
    paste(names(replicates), collapse = "\n  - ")
  ))
}

# (optional) print which columns we’ll use — super helpful for debugging
message("Using FDR columns: ctrl=", fdr_ctrl_col, " | treat=", fdr_treat_col)

# Start with betas (possibly normalized) and diff (treat - ctrl)
sel_df <- tibble::tibble(
  Gene = gsel$Gene,
  !!ctrlname := as.numeric(gsel[[ctrlname]]),
  !!treatname := as.numeric(gsel[[treatname]]),
  diff =  as.numeric(gsel[[treatname]]) - as.numeric(gsel[[ctrlname]])
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

# Order selection table by absolute effect size (|Δβ|)
if (paste0(treatname, "_FDR") %in% names(sel_df) && "diff" %in% names(sel_df)) {
  sel_df <- sel_df |>
    dplyr::mutate(effect_abs = abs(diff)) |>
    dplyr::arrange(dplyr::desc(effect_abs), .data[[paste0(treatname, "_FDR")]])
} else if ("diff" %in% names(sel_df)) {
  sel_df <- sel_df |>
    dplyr::mutate(effect_abs = abs(diff)) |>
    dplyr::arrange(dplyr::desc(effect_abs))
}

write.table(sel_df, file.path(output_dir, sprintf("selection_table.norm_%s.tsv", nm)), sep="\t", quote=FALSE, row.names=FALSE)

# -------------------------
# Counts QC plots
# -------------------------
cs <- countsummary
if ("Label" %in% names(cs)) cs$Label <- gsub("_R[12]_001$", "", cs$Label)

# gini index
p_gini <- BarView(cs, x = "Label", y = "GiniIndex", ylab="Gini index", main="Evenness of sgRNA reads") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
ggsave(file.path(output_dir, "qc_gini.png"), plot = p_gini, width = 10, height = 8, dpi = 150)

# zero counts
p_zero <- BarView(cs, x = "Label", y = "Zerocounts", fill = "#394E80",
                  ylab="Zero Count sgRNAs", main="Missed sgRNAs") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
ggsave(file.path(output_dir, "qc_zero_counts.png"), plot = p_zero, width = 10, height = 8, dpi = 150)

# mapping rates
p_map <- MapRatesView(cs) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "qc_maprates.png"), plot = p_map, width = 10, height = 8, dpi = 150)

# -------------------------
# Volcano (diff vs -log10 FDR) with cap + effect-size guides
# -------------------------
treat_fdr_col <- fdr_treat_col

# tuning knobs
eps        <- 1e-16       # protects against FDR==0

# Base dataframe: betas + FDR
vol_df <- merge(
  gsel[, c("Gene", ctrlname, treatname)],
  replicates[, c("Gene", treat_fdr_col)],
  by = "Gene", all.x = TRUE
)

# Effect size and robust FDR transforms
vol_df$diff     <- as.numeric(vol_df[[treatname]]) - as.numeric(vol_df[[ctrlname]])
vol_df$FDR_num  <- suppressWarnings(as.numeric(vol_df[[treat_fdr_col]]))
vol_df$LogFDR   <- -log10(pmax(vol_df$FDR_num, eps))
vol_df$LogFDR_c <- pmin(vol_df$LogFDR, y_cap)  # cap for presentation

# Keep rows with finite values
vol_df <- vol_df[is.finite(vol_df$diff) & is.finite(vol_df$LogFDR_c), , drop = FALSE]

# Threshold line for FDR (respect cap so it's always visible)
thr_line <- min(-log10(fdr_threshold), y_cap)

## --- MAGeCKFlute volcano (capped) ---
p_vol_flute <- ScatterView(vol_df, x = "diff", y = "LogFDR_c", label = "Gene",
                           model = "volcano", top = 20) +
  ggtitle(sprintf("%s vs %s (norm=%s, FDR≤%.2f, y≤%d)",
                  treatname, ctrlname, nm, fdr_threshold, y_cap)) +
  geom_hline(yintercept = thr_line, linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = 0,         linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = c(-effect_thr, effect_thr),
             linetype = "dotted", color = "grey50")
ggsave(file.path(output_dir, "volcano_diff_vs_neglog10FDR.png"),
       plot = p_vol_flute, width = 8, height = 6, dpi = 150)

# -------------------------
# Custom Volcano (capped) with threshold lines and labels
# -------------------------
# Top genes to label (best FDR, then largest |Δβ|)
top_lab <- vol_df |>
  dplyr::arrange(FDR_num, dplyr::desc(abs(diff))) |>
  dplyr::slice_head(n = 20)

p_vol_custom <- ggplot(vol_df, aes(x = diff, y = LogFDR_c, color = LogFDR_c)) +
  geom_point(alpha = 0.85, size = 1.6, na.rm = TRUE) +
  # threshold lines
  geom_hline(yintercept = thr_line, linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = c(-effect_thr, effect_thr),
             linetype = "dotted", color = "grey50") +
  # labels
  ggrepel::geom_text_repel(
    data = top_lab,
    aes(label = Gene),
    size = 2.5, max.overlaps = 100, box.padding = 0.3,
    point.padding = 0.2, segment.size = 0.2, show.legend = FALSE
  ) +
  theme_bw() +
  labs(
    title = sprintf("%s vs %s (norm=%s, FDR≤%.2f, y≤%d, |Δβ|≥%.2f guides)",
                    treatname, ctrlname, nm, fdr_threshold, y_cap, effect_thr),
    x = "Δβ (treat − ctrl)", y = "-log10(FDR)", color = "-log10(FDR, capped)"
  ) +
  scale_colour_gradient(low = "#7570b3", high = "#ff7f00")

ggsave(file.path(output_dir, "custom_volcano_diff_vs_neglog10FDR.png"),
       plot = p_vol_custom, width = 8, height = 6, dpi = 150)

# -------------------------
# Beta vs Beta scatter (color = treatment FDR) top-30 labeled
# -------------------------
df_scatter <- merge(
  gsel[, c("Gene", ctrlname, treatname)],
  replicates[, c("Gene", treat_fdr_col)],
  by = "Gene", all.x = TRUE
)
p_beta_scatter <- ggplot(df_scatter, aes(
  x = .data[[ctrlname]],
  y = .data[[treatname]],
  color = .data[[treat_fdr_col]]
)) +
  geom_point(alpha = 0.8, size = 1.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(title = sprintf("Beta vs Beta (%s): %s vs %s (color = %s FDR)", nm, treatname, ctrlname, treatname),
       x = sprintf("%s beta", ctrlname), y = sprintf("%s beta", treatname),
       color = sprintf("%s FDR", treatname)) +
  scale_colour_gradient(low = "#ff7f00", high = "#7570b3")

top_lab <- df_scatter |>
  dplyr::filter(!is.na(.data[[treat_fdr_col]])) |>
  dplyr::mutate(delta_beta = .data[[treatname]] - .data[[ctrlname]],
                effect_abs = abs(delta_beta)) |>
  dplyr::arrange(.data[[treat_fdr_col]], dplyr::desc(effect_abs)) |>
  dplyr::slice_head(n = 30)

p_beta_scatter_labeled <- p_beta_scatter +
  ggrepel::geom_text_repel(
    data = top_lab,
    aes(label = Gene),
    size = 2.5, max.overlaps = 100, box.padding = 0.3,
    point.padding = 0.2, segment.size = 0.2
  )

ggsave(file.path(output_dir, "beta_scatter.png"), plot = p_beta_scatter, width = 7, height = 6, dpi = 150)
ggsave(file.path(output_dir, "beta_scatter_labeled.png"), plot = p_beta_scatter_labeled, width = 7, height = 6, dpi = 150)

# -------------------------
# Significant genes (filter by treatment FDR), include control FDR if present
# -------------------------
# Add control FDR to vol_df if available
if (!is.na(fdr_ctrl_col) && fdr_ctrl_col %in% names(replicates)) {
  vol_df <- merge(
    vol_df,
    replicates[, c("Gene", fdr_ctrl_col), drop = FALSE],
    by = "Gene",
    all.x = TRUE
  )
}

# Add |Δβ| and sort by effect then FDR
sig_list <- vol_df |>
  dplyr::filter(!is.na(.data[[treat_fdr_col]]), .data[[treat_fdr_col]] <= fdr_threshold) |>
  dplyr::mutate(effect_abs = abs(diff)) |>
  dplyr::arrange(dplyr::desc(effect_abs), .data[[treat_fdr_col]])

# Friendly column names for output
treat_fdr_out <- paste0(treatname, "_FDR")
ctrl_fdr_out  <- paste0(ctrlname,  "_FDR")

names(sig_list)[names(sig_list) == ctrlname] <- paste0(ctrlname,  "_beta")
names(sig_list)[names(sig_list) == treatname] <- paste0(treatname, "_beta")
names(sig_list)[names(sig_list) == treat_fdr_col] <- treat_fdr_out
if (!is.na(fdr_ctrl_col) && fdr_ctrl_col %in% names(sig_list)) {
  names(sig_list)[names(sig_list) == fdr_ctrl_col] <- ctrl_fdr_out
}

# Final column order (includes |Δβ|)
keep_cols <- c(
  "Gene",
  paste0(ctrlname,  "_beta"),
  paste0(treatname, "_beta"),
  "diff",
  "effect_abs",
  treat_fdr_out
)
if (ctrl_fdr_out %in% names(sig_list)) keep_cols <- c(keep_cols, ctrl_fdr_out)

sig_list <- sig_list[, keep_cols, drop = FALSE]

sig_path <- file.path(output_dir, sprintf("%s_sig_hits_treatmentFDR_%.2f.tsv", proj_name, fdr_threshold))
write.table(sig_list, sig_path, sep = "\t", quote = FALSE, row.names = FALSE)
