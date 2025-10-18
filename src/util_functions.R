


# Safe numeric coercion: possible characters, strings, or empty cells (e.g.,
# such as the ones introduced by MetaFluor when pressing F5 for event marking)
# are converted to NAs.
as_num <- function(x) {
  as.numeric(as.character(x))
}

# Makes a plot of multiple calcium traces, provided that 'df' is a data frame
# containing ONLY fluorescence signals in all of its columns, while 't_vec' is
# the corresponding vector of time samples.
plot_traces <- function(df,
                        t_vec = seq(1, nrow(df), by=1),
                        title = "Traces",
                        axis_labels = c("x", "y"))
{
  # Prepare data frame in long format for plotting
  colnames(df) |>
    lapply(\(x)tibble(time = t_vec, cell = x, value = df[[x]])) |>
    bind_rows() -> long_df
  
  # Plot traces superimposed
  ggplot(long_df, aes(x = time, y = value, color = cell)) +
    geom_line(alpha = 0.6, linewidth = 1, show.legend = FALSE) +
    theme_minimal() +
    labs(title = title, x = axis_labels[1], y = axis_labels[2])
}

# Computes the baseline value of a single trace (vector 'sig')
baseline <- function(sig, type = "abs", param)
{
  if (type == "abs") {
    # Compute baseline over a fixed number of time samples
    n_baseline <- round(param, 0)
  } else if (type == "rel") {
    # Compute baseline (F0) as median of first baseline_frac portion
    n_baseline <- ceiling(length(sig) * param)
  }
  # Compute the baseline
  F0 <- median(sig[1:n_baseline], na.rm = TRUE)
  if (is.na(F0) || F0 == 0) F0 <- mean(sig[1:n_baseline], na.rm=TRUE) # fallback
  
  return(F0)
}

# Normalizes a single trace (vector 'sig')
normalize <- function(sig, b_type = "abs", b_param = 10)
{
  F0 <- baseline(sig, b_type, b_param)
  (sig - F0) / F0
}







# Feature extraction

# Detects signal collapse in a single trace (vector 'sig') by rolling median.
#
# NOTE: at the edges of the signal, where there are not enough samples to fill
# the window, rollapply() doesn't make any computation so that the returned
# vector is a ('rolling_win' - 1)-trimmed version of the original one. On the
# contrary, rollapply(... , fill = NA) keeps signal's original length by filling
# with NAs the returned vector. Try, e.g.:
#    rnorm(10, 1, 1) |> rollapply(5, median, align = "left", fill = NA)
# align="right" => index j uses: (j - rolling_win + 1) ... (j)
# align="left"  => index j uses: (j) ... (j + rolling_win - 1)
# align="center"=> index j uses: (j-(rolling_win-1)/2) ... (j+(rolling_win-1)/2)
# In this last case (default) an odd-valued length for the rolling window is
# recommended for a more accurate and predictable behavior.


# 
# norm_traces |> sapply(\(x)(x |> rollapply(width = 5,
#                                           FUN = median,
#                                           align = "right",
#                                           fill = NA) |> na.trim() |> extract() |> min()))

# thr = 0.1 empirically detected... try:
# norm_traces |>
#   lapply(\(x)convolve(x, rev(c(1, 1, 1, -1, -1, -1)/6), type = "filter")) |>
#   as.data.frame() |> plot_traces()

# c(1, -1) is equivalent to diff()

extract <- function(sig, win = c(1, 1, 1, -1, -1, -1), thr = 0.1) {
  
  # Pre-filter by a rolling median to remove possible single-time-point spikes,
  # "dropping points", or brief transients, likely due to artifacts from camera
  # sensors or alike, in any case not representing true collapse events we are
  # interested in (robust). Median removes spikes, without altering or smoothing
  # long lasting signal components.
  # Also reduces the average SD of the derivative (diff) signal
  
  # compute rolling median of previous window (aligned to right so index j uses j-rolling_win+1 .. j)
  rolling_win <- 5 # odd number here, please!
  sig |> rollapply(width = rolling_win,
                   FUN = median,
                   align = "center") -> sig
  
  
  win <- win / sum(abs(win))  # normalize
  sig |> convolve(rev(win), type = "filter") |> 
    {\(x)c(min(x),which.min(x))}() -> biggest_drop
  
  # Both rollapply and convolve (type = "filter") returns a (win-1) shorted signal
  # We need to correct for them
  biggest_drop[2] <- biggest_drop[2] + (rolling_win-1)/2 + length(win)/2
  
  ifelse(biggest_drop[1] < -thr, biggest_drop[2], NA_integer_)
}




# AUC of a single normalized trace (vector 'sig') up to collapse (if found) or
# full duration otherwise
auc <- function(sig, collapse_idx = length(sig) + 1, time_vec)
{
  end_idx <- ifelse(is.na(collapse_idx), length(sig), collapse_idx - 1)
  if (end_idx > 2) {
    trapz(time_vec[1:end_idx], sig[1:end_idx])
  } else {
    NA_real_
  }
}

# Find the first time crossing a target value (upward) by linear interpolation
cross_time <- function(x, y, target)
{
  inds <- which(!is.na(y) & y >= target)
  if (length(inds) == 0) return(NA_real_)
  idx <- inds[1]
  if (idx == 1) return(x[1])
  # linear interpolation between idx-1 and idx
  t0 <- x[idx - 1]; t1 <- x[idx]
  y0 <- y[idx - 1]; y1 <- y[idx]
  if (is.na(y0) || is.na(y1) || y1 == y0) return(t1)
  t_cross <- t0 + (target - y0) * (t1 - t0)/(y1 - y0)
  return(t_cross)
}

# Perform a gentle denoising of a 1D signal: remove single-point anomalies by
# rolling median and smooth by rolling mean.
# REMEMBER that the returned numeric vector is always 2*(win_width-1) samples
# shorter that the original one !!
denoise <- function(sig, win_width = 5) {
  sig |>
    rollapply(width = win_width, FUN = median, align = "center") |>
    rollapply(width = win_width, FUN = mean, align = "center")
}

# --- per-experiment stats -----------------------------------------------------

mean_values <- function(df) {
  # Compute the collapsing rate as percentage of collapsing cells per experiment
  df$collapse_idx |> is.na() |> {\(x)!x}() |> sum() -> not_NAs
  rate <- (not_NAs / nrow(df)) * 100
  names(rate) <- "collapse_rate"
  # Compute the mean values of all the vars
  df |> sapply(is.numeric) -> idx
  df[,idx] |> sapply(mean, na.rm = TRUE) |> append(rate, after = 3)
}

sd_values <- function(df) {
  rate <- NA
  names(rate) <- "collapse_rate"
  # Compute the SDs of all the vars
  df |> sapply(is.numeric) -> idx
  df[,idx] |> sapply(sd, na.rm = TRUE) |> append(rate, after = 3)
}

sem_values <- function(df) {
  rate <- NA
  names(rate) <- "collapse_rate"
  # Compute the SEMs of all the vars
  df |> sapply(is.numeric) -> idx
  df[,idx] |> sapply(\(x){
    sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  }) |> append(rate, after = 3)
}

median_values <- function(df) {
  rate <- NA
  names(rate) <- "collapse_rate"
  # Compute the SDs of all the vars
  df |> sapply(is.numeric) -> idx
  df[,idx] |> sapply(median, na.rm = TRUE) |> append(rate, after = 3)
}

n_values <- function(df) {
  rate <- NA
  names(rate) <- "collapse_rate"
  # Compute the SDs of all the vars
  df |> sapply(is.numeric) -> idx
  df[,idx] |> sapply(\(x)sum(!is.na(x))) |> append(rate, after = 3)
}

na_values <- function(df) {
  rate <- NA
  names(rate) <- "collapse_rate"
  # Compute the SDs of all the vars
  df |> sapply(is.numeric) -> idx
  df[,idx] |> sapply(\(x)sum(is.na(x))) |> append(rate, after = 3)
}

# --- Plot Comparisons ---------------------------------------------------------

# Plot boxplots with jittered points of two samples, also  displaying
# significance as the appropriate number of asterisks above the two boxes
ttest_boxplot <- function(x,
                          y,
                          pval,
                          group_labels = c("Group1", "Group2"),
                          xlab = "",
                          ylab = "",
                          box.width = 0.3,
                          jitter.width = 0.15,
                          point.size = 1.8,
                          palette = c("#1b9e77", "#d95f02"),
                          show.p = FALSE,
                          star.size = 6,
                          p_digits = 3)
{
  # Input sanitation
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(group_labels) < 2) group_labels <- c("Group1", "Group2")
  if (length(x) == 0 || length(y) == 0) stop("Both x and y must contain at least one (non-NA) numeric value.")
  # Remove NAs
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) == 0 || length(y) == 0) stop("After removing NAs, one group is empty.")

  
  # Map p-value to asterisks (common convention)
  p_to_stars <- function(p) {
    if (is.na(p)) return(NA_character_)
    if (p <= 0.001) return("***")
    if (p <= 0.01)  return("**")
    if (p <= 0.05)  return("*")
    return("")  # no asterisk if not significant at 0.05
  }
  stars <- p_to_stars(pval)
  
  # Build data.frame for plotting
  df <- data.frame(
    value = c(x, y),
    group = factor(rep(group_labels, times = c(length(x), length(y))), levels = group_labels)
  )
  
  # Compute annotation position (one line spanning group 1 and 2)
  y_max <- max(df$value, na.rm = TRUE)
  y_min <- min(df$value, na.rm = TRUE)
  range_val <- y_max - y_min
  if (range_val == 0) {
    # all values equal; choose a small offset
    h <- abs(y_max) * 0.1
    if (h == 0) h <- 0.5
  } else {
    h <- range_val * 0.15
  }
  y_line <- y_max + h
  y_text <- y_line + h * 0.25
  
  # Build the plot
  plt <- ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(width = box.width, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = jitter.width, size = point.size, alpha = 0.7, aes(color = group)) +
    scale_fill_manual(values = palette, guide = "none") +
    scale_color_manual(values = palette, guide = "none") +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 14)
  
  # Add significance line and stars (only if we have a star string to show; empty string means non-significant)
  if (!is.na(stars) && stars != "") {
    plt <- plt +
      geom_segment(aes(x = 1, xend = 2, y = y_line, yend = y_line), inherit.aes = FALSE, size = 0.6) +
      geom_segment(aes(x = 1, xend = 1, y = y_line, yend = y_line - h * 0.08), inherit.aes = FALSE, size = 0.6) +
      geom_segment(aes(x = 2, xend = 2, y = y_line, yend = y_line - h * 0.08), inherit.aes = FALSE, size = 0.6) +
      annotate("text", x = 1.5, y = y_text, label = stars, size = star.size, fontface = "bold")
  } else if (!is.na(pval) && show.p) {
    # Optionally show p-value even if not significant
    p_label <- paste0("p = ", formatC(pval, digits = p_digits, format = "g"))
    plt <- plt +
      geom_segment(aes(x = 1, xend = 2, y = y_line, yend = y_line), inherit.aes = FALSE, size = 0.4, linetype = "dashed") +
      annotate("text", x = 1.5, y = y_text, label = p_label, size = 6)
  } else if (is.na(pval) && show.p) {
    plt <- plt + annotate("text", x = 1.5, y = y_text, label = "p = NA", size = 6)
  }
  
  # If user wants the p-value numeric on the plot in addition to stars, add it
  if (show.p && !is.na(pval) && stars != "") {
    p_label <- paste0("p = ", signif(pval, digits = p_digits))
    plt <- plt + annotate("text", x = 1.5, y = y_text - h * 0.45, label = p_label, size = 6)
  }
  
  return(plt)
}

# Arrange a list of ggplot objects on a single canvas in an ncol x nrow grid
arrange_plots_grid <- function(plots,
                               nrow = 2,
                               ncol = 4,
                               tag_levels = NULL,   # e.g. "A" (letters) or "1" (numbers) or NULL (no tags)
                               title = NULL,
                               blank_plot = NULL)
{
  # accept a single ggplot or a list
  if (inherits(plots, "gg") || inherits(plots, "ggplot")) {
    plots <- list(plots)
  } else if (!is.list(plots)) {
    stop("'plots' must be a list of ggplot objects or a single ggplot object.")
  }
  
  total_cells <- ncol * nrow
  n_plots <- length(plots)
  
  if (n_plots > total_cells) {
    warning(sprintf("You provided %d plots but the grid has %d cells: only the first %d plots will be used.",
                    n_plots, total_cells, total_cells))
    plots <- plots[1:total_cells]
  }
  
  # create a blank ggplot if needed
  if (is.null(blank_plot)) {
    blank_plot <- ggplot() + theme_void()
  }
  
  if (n_plots < total_cells) {
    plots <- c(plots, rep(list(blank_plot), total_cells - n_plots))
  }
  
  # Build annotation arguments for patchwork::plot_annotation
  ann_args <- list()
  if (!is.null(tag_levels)) ann_args$tag_levels <- tag_levels
  if (!is.null(title)) ann_args$title <- title
  
  # Arrange with patchwork
  combined <- wrap_plots(plots, ncol = ncol, nrow = nrow)
  if (length(ann_args) > 0) {
    combined <- combined + do.call(plot_annotation, ann_args)
  }
  
  #invisible(combined)
  return(combined)
}
