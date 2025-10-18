#!/usr/bin/env Rscript

# --- Packages -----------------------------------------------------------------

# library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
#library(tibble)
library(ggplot2)
library(patchwork)
#library(r4tcpl)

# Function loading
source("./src/util_functions.R")

# --- Input Parsing ------------------------------------------------------------

# Extract command-line arguments.
in_dir <- commandArgs(trailingOnly = TRUE)[1]
out_dir <- commandArgs(trailingOnly = TRUE)[2]

# # Interactive debug (from the project root directory)
# in_dir <- "./data/in/KU_lessUV_NAC"
# out_dir <- "./data/out/KU_lessUV_NAC"

# --- Data Loading -------------------------------------------------------------

# List CSVs
file_pattern <- "_perCellReport\\.rds$"
files <- list.files(out_dir, pattern = file_pattern,
                    full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
if (length(files) == 0) {
  warning(paste("No experiment reports found in", out_dir, "directory."))
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

  # Create an empty data frame for each condition
  matrix(nrow = n_vars, ncol = 0) |> as.data.frame() -> data_points[[condition]]
  
  # Find biological replicates...
  files |> dirname() |> grep(file.path(out_dir, condition),
                             x=_, fixed = TRUE) -> idx
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
    out_summary <- file.path(out_dir, condition,
                           paste0(exp_id, "_ExperimentReport"))
    write.csv(summary_stats,
              file = paste0(out_summary, ".csv"),
              row.names = TRUE)
    message("Experiment report written to: ", out_summary)
    
    # Store in the master list
    data_points[[condition]][[exp_id]] <- summary_stats[, "Mean"]
  }
  # Add row names to each dataframe
  features <- rownames(summary_stats)
  rownames(data_points[[condition]]) <- features
}

# --- Hypothesis Testing ------------------------------------------------------- 

# Load the contrasts of interest
list.files(in_dir,
           pattern = "^comp.*\\.txt$",
           full.names = TRUE,
           ignore.case = TRUE) -> comp_files

if (length(comp_files) == 0) {
  cat("WARNING: comparison definition unavailable... skip this batch!")
  quit(status = 1)
} else if (length(comp_files) > 1) {
  cat("WARNING: multiple comparison definitions available: taking one...")
}
comp_files[1] |> read.delim(header = FALSE,
                         sep = "\n",
                         blank.lines.skip = TRUE,
                         comment.char = "#") |> unlist() -> comparisons

# Make all the planned comparisons
for (comp in comparisons) {
  
  # Parse the contrast
  comp |> strsplit(" -- ") |> unlist() |> {\(x)x[1]}() -> cond
  comp |> strsplit(" -- ") |> unlist() |> {\(x)x[2]}() -> ref
  
  # Perform t-test, catch errors (e.g., zero variance, low sample size, ...)
  features |> sapply(function(feature) {
    tryCatch(
      t.test(x = data_points[[cond]][feature,],
             y = data_points[[ref]][feature,],
             paired = FALSE,
             var.equal = FALSE,
             alternative = "two.sided")$p.val,
      error = function(e) NA_real_
      )
    }) -> pvals
  
  
  
  features |> lapply(function(feature) {
    ttest_boxplot(data_points[[ref]][feature,],
                  data_points[[cond]][feature,],
                  pvals[feature],
                  group_labels = c(ref, cond),
                  ylab = feature,
                  show.p = TRUE)
    }) -> ggplots_list
  
  ggplots_list |> arrange_plots_grid(nrow = 2,
                                     ncol = 4,
                                     tag_levels = "A",
                                     title = comp) -> combined
  
  # Save the Plot
  r4tcpl::savePlots(
    \(){print(combined)},
    width_px = 2000,
    figure_Name = paste0(comp, "_Combined"),
    figure_Folder = out_dir)

  
}






