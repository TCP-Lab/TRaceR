#!/usr/bin/env Rscript

# --- Packages -----------------------------------------------------------------

#library(r4tcpl)
library(stats)
library(dplyr, warn.conflicts = FALSE) |> suppressWarnings()
library(tidyr, warn.conflicts = FALSE) |> suppressWarnings()
#library(tibble)
library(ggplot2, warn.conflicts = FALSE) |> suppressWarnings()
library(patchwork)

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

message(" \nOpening Batch: ", out_dir)

# List CSVs
file_pattern <- "_perCellReport\\.rds$"
files <- list.files(out_dir, pattern = file_pattern,
                    full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
if (length(files) == 0) {
  warning("No experiment reports found in ", out_dir, " directory.")
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
  
  message(" \nProcessing Condition: ", condition)

  # Create an empty data frame for each condition to store in the master list
  matrix(nrow = n_vars, ncol = 0) |> as.data.frame() -> data_points[[condition]]
  
  # Find biological replicates...
  files |> dirname() |> grep(paste0("^", file.path(out_dir, condition), "$"),
                             x=_, fixed = FALSE) -> idx
  message("  > ", length(idx), " experiments found.")
  # ...and cycle over them
  for (i in idx) {
    
    files[i] |> basename() |> sub(file_pattern, "", x=_) -> exp_id
    message("   >> Processing: ", exp_id)
    
    # --- Experiment Stats -----------------------------------------------------
    
    # Compute summary statistics across cells (mean, SD, SEM, ...) for each var
    readRDS(files[i]) -> summary_tbl
    data.frame(Mean   = summary_tbl |> get_stat(mean),
               SD     = summary_tbl |> get_stat(sd),
               SEM    = summary_tbl |> get_stat(sem),
               Median = summary_tbl |> get_stat(median),
               IQR    = summary_tbl |> get_stat(IQR),
               N      = summary_tbl |> get_stat(get_size),
               NAs    = summary_tbl |> get_stat(get_nas)) -> summary_stats
    
    # Save experiment report
    out_summary <- file.path(out_dir, condition,
                           paste0(exp_id, "_ExperimentReport"))
    write.csv(summary_stats,
              file = paste0(out_summary, ".csv"),
              row.names = TRUE)
    message("    >>> Experiment report written to:\n",
            "      ", out_summary, " (csv)")
    
    # Store in the master list
    data_points[[condition]][[exp_id]] <- summary_stats[, "Mean"]
  }
  # Add row names to each dataframe
  features <- rownames(summary_stats)
  rownames(data_points[[condition]]) <- features
  
  # Export data points for possible third-party replotting
  con <- file(file.path(out_dir,
                        paste0("Batch_(",basename(out_dir),")_MeanPoints.csv")),
              open = "wt")
  for (i in seq_along(data_points)) {
    writeLines(paste0("### Condition ", names(data_points)[i]), con)
    write.table(data_points[[i]], con, sep = ",",
                row.names = TRUE, col.names = NA)
    # "By default there is no column name for a column of row names.
    # If col.names = NA and row.names = TRUE a blank column name is added..."
    writeLines("", con)
  }
  close(con)
}

# --- Hypothesis Testing ------------------------------------------------------- 

message(" \nHypothesis Testing ", out_dir)

# Load the contrasts of interest
list.files(in_dir,
           pattern = "^comp.*\\.txt$",
           full.names = TRUE,
           ignore.case = TRUE) -> comp_files

if (length(comp_files) == 0) {
  warning("Comparison definition unavailable... skip this batch!")
  quit(status = 1)
} else if (length(comp_files) > 1) {
  warning("Multiple comparison definitions available: taking one...")
}
comp_files[1] |> read.delim(header = FALSE,
                            sep = "\n",
                            blank.lines.skip = TRUE,
                            comment.char = "#") |> unlist() -> comparisons
message(length(comparisons), " contrasts to run")

# Select the Features Of Interest (fois)
fois <- c("F0", "max_norm", "collapse_rate", "collapse_times",
          "AUCs", "onset_time", "rise_time_10_90", "half_height_width")

# Make all the planned comparisons
for (comp in comparisons) {
  
  message("  > Contrast: ", comp)
  
  # Parse the contrast
  comp |> strsplit(" -- ") |> unlist() |> {\(x)x[1]}() -> cond
  comp |> strsplit(" -- ") |> unlist() |> {\(x)x[2]}() -> ref
  
  # Check group name consistency
  if(any(! c(cond, ref) %in% names(data_points))) {
    message(paste0("WARNING: Inconsistent group names in comparison (",
                   comp, "). Contrast skipped!"))
    next
  }
  
  # Perform t-test, catching errors (e.g., zero variance, low sample size, ...)
  message("   >> Doing the t-test")
  
  # arcsine square-root transformation to stabilize variance
  asin_sqrt <- function(x) {asin(sqrt(x/100))}
  
  fois |> sapply(function(feature) {
    x <- as.numeric(data_points[[ref]][feature,])
    y <- as.numeric(data_points[[cond]][feature,])
    if (feature == "collapse_rate") {
      x <- asin_sqrt(x)
      y <- asin_sqrt(y)
    }
    tryCatch(
      t.test(x = x,
             y = y,
             paired = FALSE,
             var.equal = FALSE,
             alternative = "two.sided")$p.val,
      error = function(e) NA_real_
      )
    }) -> pvals_t
  
  # Perform Wilcoxon Rank-Sum (aka Mann-Whitney U) test
  message("   >> Doing the U-test")
  fois |> sapply(function(feature) {
    tryCatch(
      wilcox.test(x = as.numeric(data_points[[ref]][feature,]),
                  y = as.numeric(data_points[[cond]][feature,]),
                  paired = FALSE,
                  alternative = "two.sided")$p.val,
      error = function(e) NA_real_
    )
  }) -> pvals_U
  
  # Draw a box plot per feature
  message("   >> Making box plots")
  fois |> lapply(function(feature) {
    ttest_boxplot(x = data_points[[ref]][feature,],
                  y = data_points[[cond]][feature,],
                  group_labels = c(ref, cond),
                  t.pval = pvals_t[feature],
                  U.pval = pvals_U[feature],
                  ylab = feature,
                  show.p = TRUE)
    }) -> ggplots_list
  
  # Merge all features into a per-contrast graphical report
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

