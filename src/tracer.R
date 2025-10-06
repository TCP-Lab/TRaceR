#!/usr/bin/env Rscript

# --- Packages -----------------------------------------------------------------

library(tidyverse)
library(zoo)    # for rollapply
library(pracma) # for trapz

# Function loading
source("./src/util_functions.R")

# --- Input Parsing ------------------------------------------------------------


in_dir <- "."                 # directory with CSV files
file_pattern <- "\\.csv$"     # pattern to identify CSV files
out_dir <- "analysis_output"  # where to save plots/reports
baseline_frac <- 0.10         # fraction of initial samples to compute F0 (default first 10%)


dir.create(out_dir, showWarnings = FALSE) # Move to Bash wrapper






# --- Data Loading -------------------------------------------------------------

# list CSVs
files <- list.files(in_dir, pattern = file_pattern, full.names = TRUE)
if (length(files) == 0) {
  warning(paste("No CSV files found in", in_dir, "input directory."))
  next
} # For final implementation within a super FOR cycling over multiple input dirs



# To collect final stats
all_summary_list <- list()



for (fpath in files) {
  
  
  # DEBUG !!!!!
  fpath <- files[1]
  
  message("Processing: ", fpath)
  
  # Read entire CSV skipping the first row, and using the second one for heading
  fpath |> read.csv(header = TRUE,
                    skip = 1,
                    check.names = TRUE,
                    stringsAsFactors = FALSE) |>
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
  ROIs <- colnames(raw_traces)
  n_samples <- length(time_vec)
  n_cells <- length(ROIs)
  
  
  
  
  # --- Plot raw traces --------------------------------------------------------
  
  p_raw <- plot_traces(raw_traces, time_vec,
                       title = paste0("Raw traces: ", basename(fpath)),
                       axis_labels = c("Time (s)", "Raw fluorescence (a.u.)"))
  p_raw
  

  # --- Normalize traces -------------------------------------------------------
  
  
  raw_traces |> lapply(normalize, b_type = "abs", b_param = 50) |>
    as.data.frame() -> norm_traces
  
  raw_traces |> lapply(normalize, b_type = "rel", b_param = baseline_frac) |>
    as.data.frame() -> norm_traces

  p_norm <- plot_traces(norm_traces, time_vec,
                        title = paste0("Normalized traces (dF/F): ", basename(fpath)),
                        axis_labels = c("Time (s)", "dF/F"))
  p_norm
  
  
  # --- Collapse Detection -----------------------------------------------------
  
  
  norm_traces |> sapply(extract) -> collapse_idx
  
  collapse_idx |> sapply(\(x)ifelse(is.na(x), NA_real_, time_vec[x])) -> collapse_times
  
  # --- Descriptive Stats ------------------------------------------------------
  
  # Compute the Area Under the Curve (AUC) up to collapse
  ROIs |> sapply(\(x)auc(norm_traces[[x]], collapse_idx[x], time_vec)) -> AUCs
  
  # Find maximum of the normalized signal, after gentle denoising
  # (remove single-point anomalies and mean-smooth)
  norm_traces |>
    sapply(\(x) x |>
             rollapply(width = 5, FUN = median, align = "center") |>
             rollapply(width = 5, FUN = mean, align = "center") |>
             max()) -> max_norm

  # Rising Time???

  # store
  summary_tbl <- tibble(cell = ROIs)
  summary_tbl$F0 <- raw_traces |> sapply(baseline, type = "abs", param = 50)
  summary_tbl$collapse_idx <- collapse_idx
  summary_tbl$collapse_times <- collapse_times
  summary_tbl$AUCs <- AUCs
  summary_tbl$max_norm <- max_norm


  
  

  
  # --- Save per-cell summary --------------------------------------------------
  
  out_summary_csv <- file.path(out_dir, paste0(basename(fpath), "_per_cell_stats.csv"))
  write.csv(summary_tbl, out_summary_csv, row.names = FALSE)
  message("Per-cell stats written to: ", out_summary_csv)
  
  
  # Save histograms for each stat
  stat_names <- c("F0", "collapse_times", "AUCs", "max_norm")
  hist_plots <- list()
  for (s in stat_names) {
    df_plot <- summary_tbl |> select(cell, !!sym(s)) |> rename(value = !!sym(s))
    p_hist <- ggplot(df_plot, aes(x = value)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.9) +
      theme_minimal() +
      labs(title = paste0(basename(fpath), " — histogram: ", s),
           x = s, y = "Count")
    #fname <- file.path(out_dir, paste0(basename(fpath), "_hist_", s, ".png"))
    #ggsave(filename = fname, plot = p_hist, width = 7, height = 5)
    hist_plots[[s]] <- p_hist
  }
  
  hist_plots[[4]]
  
  
  # --- Save per-experiment summary -------------------------------------------- 
  
  # Compute summary statistics across cells (mean, sd, sem, median, n) for each stat
  summary_tbl |>
    summarize_at(vars(F0, collapse_times, AUCs, max_norm),
                 list(mean = ~mean(., na.rm = TRUE),
                      sd   = ~sd(., na.rm = TRUE),
                      sem  = ~ (sd(., na.rm = TRUE) / sqrt(sum(!is.na(.)))),
                      median = ~median(., na.rm = TRUE),
                      n = ~sum(!is.na(.)))) |>
    pivot_longer(everything(), names_to = "stat_metric", values_to = "value") -> summary_stats
  
  
  
  out_summary_stats_csv <- file.path(out_dir, paste0(basename(fpath), "_summary_stats.csv"))
  write.csv(summary_stats, out_summary_stats_csv, row.names = FALSE)
  message("Summary stats written to: ", out_summary_stats_csv)
  

  
  
  # # Save combined per-file RDS of results for programmatic use
  # saveRDS(list(per_cell = summary_tbl, summary_stats = summary_stats, norm_traces = long_norm),
  #         file = file.path(out_dir, paste0(basename(fpath), "_results.rds")))
  
  # Collect global list
  all_summary_list[[basename(fpath)]] <- list(per_cell = summary_tbl, summary_stats = summary_stats)
}



# # Optionally combine all files into a single table and write
# combined_per_cell <- bind_rows(lapply(all_summary_list, function(x) x$per_cell), .id = "source_file")
# write.csv(combined_per_cell, file.path(out_dir, "combined_per_cell_stats.csv"), row.names = FALSE)
# message("Combined per-cell stats saved to: ", file.path(out_dir, "combined_per_cell_stats.csv"))
# 
# 
