# Still to be reviewed...
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
  
  
  out_summary_stats_csv <- file.path(out_dir, paste0(exp_id, "_summary_stats.csv"))
  write.csv(summary_stats, out_summary_stats_csv, row.names = FALSE)
  message("Summary stats written to: ", out_summary_stats_csv)
  
  # # Save combined per-file RDS of results for programmatic use
  # saveRDS(list(per_cell = summary_tbl, summary_stats = summary_stats, norm_traces = long_norm),
  #         file = file.path(out_dir, paste0(exp_id, "_results.rds")))
  
  # Collect global list
  all_summary_list[[exp_id]] <- list(per_cell = summary_tbl, summary_stats = summary_stats)
}