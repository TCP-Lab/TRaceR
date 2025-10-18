#!/usr/bin/env Rscript

# --- Packages -----------------------------------------------------------------

# library(tidyverse)

# library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(tibble)
#library(r4tcpl)

# Function loading
source("./src/util_functions.R")

# --- Input Parsing ------------------------------------------------------------

# Extract command-line arguments.
target_dir <- commandArgs(trailingOnly = TRUE)[1]

# Interactive debug (from the project root directory)
target_dir <- "./data/out/KU_lessUV_NAC"

# --- Data Loading -------------------------------------------------------------

# List CSVs
file_pattern <- "_perCellReport\\.rds$"
files <- list.files(target_dir, pattern = file_pattern,
                    full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
if (length(files) == 0) {
  warning(paste("No experiment reports found in", target_dir, "directory."))
  quit(status = 1)
}

# Initialize the master list
data_points <- list()
# Retrieve the number of considered features
# WARNING!! it is: ncol() - 1 + 1 == ncol()
#           -1 (because of ROI names)
#           +1 (because of 'collapse_rate' addition)
# ...modify accordingly.
files[1] |> readRDS() |> ncol() -> n_vars

# Find experimental conditions...
files |> dirname() |> basename() |> unique() -> conditions
# ...and cycle over them
for (condition in conditions) {
  
  message("Processing: ", condition)

  # Create an empty tibble for each condition
  matrix(nrow = n_vars, ncol = 0) |> as.data.frame() -> data_points[[condition]]
  
  # Find biological replicates...
  files |> grep(file.path(target_dir, condition, "/"), x=_, fixed = TRUE) -> idx
  # ...and cycle over them
  for (i in idx) {
    
    files[i] |> basename() |> sub(file_pattern, "", x=_) -> exp_id
    
    # --- Experiment Stats -----------------------------------------------------
    
    # Compute summary statistics across cells (mean, SD, SEM, ...) for each var
    readRDS(files[i]) -> summary_tbl
    data.frame(Mean   = summary_tbl |> mean_values(),
               SD     = summary_tbl |> sd_values(),
               SEM    = summary_tbl |> sem_values(),
               Median = summary_tbl |> median_values(),
               N      = summary_tbl |> n_values(),
               NAs    = summary_tbl |> na_values()) -> summary_stats
    
    # --- Save Experiment Report -------------------------------------------------
  
    # Save experiment report
    out_summary <- file.path(target_dir, condition,
                           paste0(exp_id, "_ExperimentReport"))
    write.csv(summary_stats,
              file = paste0(out_summary, ".csv"),
              row.names = TRUE)
    message("Experiment report written to: ", out_summary)
    
    # Store in the master list
    data_points[[condition]][[exp_id]] <- summary_stats[, "Mean"]
  }
  # Add row names to each dataframe
  rownames(data_points[[condition]]) <- rownames(summary_stats)
}











# Still to be reviewed... for histograms
if (FALSE) {
  # Save histograms for each stat
  stat_names <- c("F0", "collapse_times", "AUCs", "max_norm")
  hist_plots <- list()
  for (s in stat_names) {
    df_plot <- summary_tbl |> select(cell, !!sym(s)) |> rename(value = !!sym(s))
    p_hist <- ggplot(df_plot, aes(x = value)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.9) +
      theme_minimal() +
      labs(title = paste0(exp_id, " — histogram: ", s),
           x = s, y = "Count")
    #fname <- file.path(out_dir, paste0(exp_id, "_hist_", s, ".png"))
    #ggsave(filename = fname, plot = p_hist, width = 7, height = 5)
    hist_plots[[s]] <- p_hist
  }
}




