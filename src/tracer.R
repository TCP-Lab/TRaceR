#!/usr/bin/env Rscript

# --- Packages -----------------------------------------------------------------

# library(tidyverse)

library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(zoo, warn.conflicts = FALSE) |> suppressWarnings()     # rollapply()
library(pracma, warn.conflicts = FALSE) |> suppressWarnings()  # trapz()
#library(r4tcpl)

# Function loading
source("./src/util_functions.R")

# --- Input Parsing ------------------------------------------------------------

# Extract command-line arguments.
in_dir <- commandArgs(trailingOnly = TRUE)[1]
out_dir <- commandArgs(trailingOnly = TRUE)[2]

# # Interactive debug (from the project root directory)
# in_dir <- "./data/in/KU_lessUV_NAC/KU"
# out_dir <- "./data/out/KU_lessUV_NAC/KU"

# b_type <- "rel"
# b_param <- 0.10   # fraction of initial samples to compute F0 (first 10%)
b_type <- "abs"
b_param <- 10

# --- Data Loading -------------------------------------------------------------

# List CSVs
file_pattern <- "\\.csv$"
files <- list.files(in_dir, pattern = file_pattern, full.names = TRUE)
if (length(files) == 0) {
  warning(paste("No CSV files found in", in_dir, "input directory."))
  quit(status = 1)
}

# To collect final stats
all_summary_list <- list()

# Loop over experiments
for (fpath in files) {
  
  message("Processing: ", fpath)
  fpath |> basename() |> sub(file_pattern, "", x=_) -> exp_id
  
  # Read the entire CSV. Mind the special format:
  # - , for decimal point
  # - ; as field separator
  # - skip the first row
  # - using the second row for heading
  # - fileEncoding "UTF-16LE" (with BOM)
  fpath |> read.table(header = TRUE,
                      sep = ";",
                      dec = ",",
                      skip = 1,
                      check.names = TRUE,
                      stringsAsFactors = FALSE,
                      fileEncoding = "UTF-16LE") |>
    lapply(as_num) |> as.data.frame() -> raw_traces
  # Sanitize heading labels
  colnames(raw_traces)[1] <- "Time"
  colnames(raw_traces) |> sub("\\.*$", "", x=_) -> colnames(raw_traces)
  
  # Expect first column to be the vector of time samples
  time_vec <- raw_traces$Time
  if (any(is.na(time_vec))) {
    warning("NAs in Time vector!")
  }
  # Signals are the remaining columns
  raw_traces <- raw_traces[,-1]
  # Check for and remove constantly-zero traces
  raw_traces |> sapply(\(x)all(x==0)) |> which() -> null_traces
  if (length(null_traces) > 0) {
    raw_traces <- raw_traces[,-null_traces]
  }
  # Store ROI names
  ROIs <- colnames(raw_traces)
  
  # --- Plot Raw Traces --------------------------------------------------------
  
  p_raw <- plot_traces(raw_traces, time_vec,
                       title = paste("Raw traces:", exp_id),
                       axis_labels = c("Time (s)", "Raw fluorescence (a.u.)"))
  # Save the Plot
  r4tcpl::savePlots(
    \(){print(p_raw)},
    width_px = 2000,
    figure_Name = paste0(exp_id, "_Raw_traces"),
    figure_Folder = out_dir)
  
  # --- Normalize Traces -------------------------------------------------------
  
  raw_traces |> lapply(normalize, b_type, b_param) |> as.data.frame() -> norm_traces
  
  p_norm <- plot_traces(norm_traces, time_vec,
                        title = paste("Normalized traces (dF/F0):", exp_id),
                        axis_labels = c("Time (s)", "dF/F0"))
  # Save the Plot
  r4tcpl::savePlots(
    \(){print(p_norm)},
    width_px = 2000,
    figure_Name = paste0(exp_id, "_Norm_traces"),
    figure_Folder = out_dir)
  
  # --- Collapse Detection -----------------------------------------------------
  
  norm_traces |> sapply(extract) -> collapse_idx
  
  collapse_idx |> sapply(\(x)ifelse(is.na(x), NA_real_, time_vec[x])) -> collapse_times
  
  # Check for Debug and thr fine-tuning
  # norm_traces |>
  #   lapply(\(x)convolve(x, rev(c(1, 1, 1, -1, -1, -1)/6), type = "filter")) |>
  #   as.data.frame() |> plot_traces()
  
  # --- Descriptive Stats ------------------------------------------------------
  
  # Compute the Area Under the Curve (AUC) up to collapse
  ROIs |> sapply(\(x)auc(norm_traces[[x]], collapse_idx[x], time_vec)) -> AUCs
  
  # Find the maximum of the normalized (denoised) signal
  norm_traces |> sapply(\(x) x |> denoise() |> max()) -> max_norm
  
  # 2nd-order Kinetics:
  # ~~~~~~~~~~~~~~~~~~
  # Onset Time
  # 10-90% Rise Time
  # Half-height width
  
  # Onset Time, defined as the time before the (denoised) signal crosses the
  # 5-times noise-level of the baseline (five-sigma criterion)
  
  # Baseline of the normalize traces
  F0_Norm <- norm_traces |> sapply(baseline, b_type, b_param)
  # Baseline noise (5x)
  baseline_length <- ifelse(b_type == "abs", b_param,
                            ceiling(length(time_vec)*b_param))
  norm_traces[1:baseline_length,] |> apply(2,sd) |> {\(x)5*x}() -> t5s_vals
  
  ROIs |> sapply(\(x){
    w <- 5
    norm_traces[[x]] |> denoise(win_width = w) -> smooth_norm
    time_smooth <- time_vec[w:(length(time_vec)-w+1)]
    t5s <- cross_time(time_vec, smooth_norm, t5s_vals[x])
  }) -> onset_time
  
  # Find the 10-90% Rise Time (after gentle denoising)
  t10_vals <- F0_Norm + 0.10*(max_norm - F0_Norm)
  t90_vals <- F0_Norm + 0.90*(max_norm - F0_Norm)
  
  ROIs |> sapply(\(x){
    w <- 5
    norm_traces[[x]] |> denoise(win_width = w) -> smooth_norm
    time_smooth <- time_vec[w:(length(time_vec)-w+1)]
    t10 <- cross_time(time_vec, smooth_norm, t10_vals[x])
    t90 <- cross_time(time_vec, smooth_norm, t90_vals[x])
    ifelse(is.na(t10) || is.na(t90), NA_real_, t90 - t10)
    }) -> rise_time_10_90
  
  # store Descriptive Stats
  summary_tbl <- tibble(cell = ROIs)
  summary_tbl$F0 <- raw_traces |> sapply(baseline, b_type, b_param)
  summary_tbl$collapse_idx <- collapse_idx
  summary_tbl$collapse_times <- collapse_times
  summary_tbl$AUCs <- AUCs
  summary_tbl$max_norm <- max_norm
  summary_tbl$rise_time_10_90 <- rise_time_10_90
  
  # --- Save per-cell summary --------------------------------------------------
  
  out_summary_csv <- file.path(out_dir, paste0(exp_id, "_perCell_stats.csv"))
  write.csv(summary_tbl, out_summary_csv, row.names = FALSE)
  message("Per-cell stats written to: ", out_summary_csv)
  
}









