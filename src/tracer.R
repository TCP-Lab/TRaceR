#!/usr/bin/env Rscript

# --- Packages -----------------------------------------------------------------

#library(r4tcpl)
#library(tidyverse)
library(dplyr, warn.conflicts = FALSE) |> suppressWarnings()
library(ggplot2, warn.conflicts = FALSE) |> suppressWarnings()
library(patchwork) # wrap_plots()
library(zoo, warn.conflicts = FALSE) |> suppressWarnings()    # rollapply()
library(pracma, warn.conflicts = FALSE) |> suppressWarnings() # trapz(), findpeaks()
#library(signal) # Butterworth filter

# Function loading
source("./src/util_functions.R")

# --- Input Parsing ------------------------------------------------------------

# Extract command-line arguments.
in_dir <- commandArgs(trailingOnly = TRUE)[1]
out_dir <- commandArgs(trailingOnly = TRUE)[2]

# # Interactive debug (from the project root directory)
# in_dir <- "./data/in/KU_lessUV_NAC/KU"
# out_dir <- "./data/out/KU_lessUV_NAC/KU"
#
# in_dir <- "./data/in/KU_lessUV_NAC/DMSO_NAC"
# out_dir <- "./data/out/KU_lessUV_NAC/DMSO_NAC"
#
# in_dir <- "./data/in/KU_Gd_CBX/KU"
# out_dir <- "./data/out/KU_Gd_CBX/KU"
#
# in_dir <- "./data/in/KU_OMO/KU"
# out_dir <- "./data/out/KU_OMO/KU"

# b_type <- "rel"
# b_param <- 0.10   # Fraction of initial samples to compute F0 (first 10%)
b_type <- "abs"
b_param <- 10     # Number of initial samples to be considered for F0

# --- Data Loading -------------------------------------------------------------

message(" \nOpening: ", in_dir)

# List CSVs
file_pattern <- "\\.csv$"
files <- list.files(in_dir, pattern = file_pattern,
                    full.names = TRUE, ignore.case = TRUE)
if (length(files) == 0) {
  warning("  > No CSV files found in ", in_dir, " input directory.")
  quit(status = 1)
} else {
  message("  > ", length(files), " CSV files found.")
}

# To collect final stats
all_summary_list <- list()

# Loop over experiments
for (fpath in files) {
  
  message("  >> Processing: ", basename(fpath))
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
    warning("  >>> NAs in Time vector!")
  }
  time_vec |> diff() |> median() -> dt # Sampling time
  # Signals are the remaining columns
  raw_traces <- raw_traces[,-1]
  # Search and destroy constantly-zero traces
  raw_traces |> sapply(\(x)all(x==0)) |> which() -> null_traces
  if (length(null_traces) > 0) {
    raw_traces <- raw_traces[,-null_traces]
  }
  # Store ROI names
  ROIs <- colnames(raw_traces)
  
  message("  >>> ", ncol(raw_traces), " ROIs x ",
          nrow(raw_traces), " time samples",
          " (dt = ", round(dt, digits = 3), " s) - fluoIQR (a.u.): ",
          round(IQR(unlist(raw_traces), na.rm = TRUE), digits = 2))
  
  # --- Plot Raw Traces --------------------------------------------------------
  
  p_raw <- plot_traces(raw_traces, time_vec,
                       title = paste("Raw traces:", exp_id),
                       axis_labels = c("Time (s)", "Raw fluorescence (a.u.)"))
  
  # Save the Plot
  r4tcpl::savePlots(
    \(){print(p_raw)},
    width_px = 2000,
    figure_Name = paste0(exp_id, "_1RawTraces"),
    figure_Folder = out_dir,
    pdf_out = TRUE,
    png_out = TRUE) # This (default) option sometimes (!) triggers the following
                    # annoying warning on stderr: "Fontconfig warning: using
                    # without calling FcInit()" which I cannot suppress in any way.
  
  # --- Normalize Traces -------------------------------------------------------
  
  raw_traces |> lapply(normalize, b_type, b_param) |> as.data.frame() -> norm_traces
  
  p_norm <- plot_traces(norm_traces, time_vec,
                        title = paste("Normalized traces (dF/F0):", exp_id),
                        axis_labels = c("Time (s)", "dF/F0"))
  # Save the Plot
  r4tcpl::savePlots(
    \(){print(p_norm)},
    width_px = 2000,
    figure_Name = paste0(exp_id, "_2NormTraces"),
    figure_Folder = out_dir)
  
  # --- Collapse Detection -----------------------------------------------------
  
  mw <- 5 # width of the median pre-filter (median_width)
  step_win <- c(1, 1, 1, -1, -1, -1)
  protect <- 200
  thr <- 10
  
  
  raw_traces |> sapply(\(x) x |> step_convolve() |> IQR()) |> median() -> average_noise
  message("  >>> Average noise estimate: ", round(average_noise, digits = 3))
  thr <- 7*average_noise
  
  # Note here the use of RAW traces to allow for absolute thr values!!
  # It is virtually impossible to define a universal threshold for normalized
  # traces because the range of the y-axis--being a function of the baseline
  # values (F0)--is much more variable. However, the drop points do not change
  # as a result of normalization.
  raw_traces |> sapply(extract, mw, step_win, protect, thr) -> collapse_idx
  
  collapse_idx |> sapply(\(x)ifelse(is.na(x), NA_real_, time_vec[x])) -> collapse_times
  
  # Checkpoint for debugging and threshold fine-tuning
  raw_traces |> lapply(step_convolve, mw, step_win) |> as.data.frame() -> conv_traces
  conv_traces |> select(which(!is.na(collapse_idx))) |> mutate(threshold = -thr) |>
    plot_traces(title = paste0("Kept: ", sum(!is.na(collapse_idx))),
                axis_labels = NULL) -> p_conv_kept
  conv_traces |> select(which(is.na(collapse_idx))) |> mutate(threshold = -thr) |>
    plot_traces(title = paste0("Discarded: ", sum(is.na(collapse_idx))),
                axis_labels = NULL) -> p_conv_disc
  
  arrange_plots_grid(list(p_conv_kept, p_conv_disc),
                     nrow = 2, ncol = 1, tag_levels = NULL,
                     title = paste0("Total Traces: ", length(collapse_idx)),
                     blank_plot = NULL) -> p_conv
  
  # Save the Plot
  r4tcpl::savePlots(
    \(){print(p_conv)},
    width_px = 2000,
    figure_Name = paste0(exp_id, "_3ConvTraces"),
    figure_Folder = out_dir)
  
  # --- Descriptive Stats ------------------------------------------------------
  
  # Compute the Area Under the Curve (AUC) up to collapse
  ROIs |> sapply(\(x)auc(norm_traces[[x]], collapse_idx[x], time_vec)) -> AUCs
  
  # Find the maximum of the normalized (and denoised) signal
  norm_traces |> sapply(\(x) x |> denoise() |> max()) -> max_norm
  
  # 2nd-order Kinetics:
  # ~~~~~~~~~~~~~~~~~~
  # Onset Time
  # 10-90% Rise Time
  # Half-height width
  
  # Gentle denoising of the normalized signals
  w <- 5
  ROIs |> sapply(\(x)norm_traces[[x]] |> denoise(win_width = w)) -> smooth_norm
  time_vec[w:(length(time_vec)-w+1)] -> time_smooth
  
  # Onset Time, defined as the time before the (denoised) signal crosses the
  # 5-times noise-level of the baseline (five-sigma criterion)
  baseline_length <- ifelse(b_type == "abs", b_param,
                            ceiling(length(time_vec)*b_param))
  norm_traces[1:baseline_length,] |> apply(2,sd) |> {\(x)5*x}() -> t5s_vals
  
  ROIs |> sapply(\(x){
    cross_time(time_smooth, smooth_norm[,x], t5s_vals[x])}) -> onset_time
  
  # Find the 10-90% Rise Time
  t10_vals <- 0.10 * max_norm
  t90_vals <- 0.90 * max_norm
  ROIs |> sapply(\(x){
    t10 <- cross_time(time_smooth, smooth_norm[,x], t10_vals[x])
    t90 <- cross_time(time_smooth, smooth_norm[,x], t90_vals[x])
    ifelse(is.na(t10) || is.na(t90), NA_real_, t90 - t10)
    }) -> rise_time_10_90
  
  # Half-height width, defined as the time interval between the point at which
  # the signal reaches half its maximum value and the (eventual) drop-off.
  t50_vals <- 0.50 * max_norm
  ROIs |> sapply(\(x){
    t50 <- cross_time(time_smooth, smooth_norm[,x], t50_vals[x])
    t_end <- ifelse(is.na(collapse_times),
                    time_vec[length(time_vec)], collapse_times)
    ifelse(is.na(t50), NA_real_, t_end - t50)
  }) -> half_height_width
  
  # store Descriptive Stats
  summary_tbl <- tibble(cell = ROIs)
  summary_tbl$F0 <- raw_traces |> sapply(baseline, b_type, b_param)
  summary_tbl$collapse_idx <- collapse_idx
  summary_tbl$collapse_times <- collapse_times
  summary_tbl$AUCs <- AUCs
  summary_tbl$max_norm <- max_norm
  summary_tbl$onset_time <- onset_time
  summary_tbl$rise_time_10_90 <- rise_time_10_90
  summary_tbl$half_height_width <- half_height_width
  
  # --- Save per-Cell Report ---------------------------------------------------
  
  out_summary <- file.path(out_dir, paste0(exp_id, "_perCellReport"))
  
  # Save per-cell stats
  write.csv(summary_tbl,
            file = paste0(out_summary, ".csv"),
            row.names = FALSE)
  
  # Save RDS of results for programmatic use
  saveRDS(summary_tbl,
          file = paste0(out_summary, ".rds"))

  message("  >>> Per-cell stats written to:\n",
          "      ", out_summary, " (csv + rds)")
  
}

