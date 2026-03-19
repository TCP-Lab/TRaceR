
# --- General ------------------------------------------------------------------

# Robust wrapper to make a function completely silent while still returning the
# result. This allow suppressing messages, warnings, printed output, and package
# startup messages... but not ERRORS, which will still stop execution, as it is
# usually desirable.
silent <- function(expr) {
  invisible(capture.output(
    suppressMessages(
      suppressWarnings(
        result <- eval.parent(substitute(expr))
      )
    )
  ))
  result
}

# Safe numeric coercion: possible characters, strings, or empty cells (e.g.,
# such as the ones introduced by MetaFluor when pressing F5 for event marking)
# are converted to NAs.
as_num <- function(x)
{
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

# Find the n-th top-most max within a vector of numbers
nth.max <- function(x, n)
{
  # Double '-' because 'decreasing' options is not available for partial sorting
  -sort(-x, partial = n)[n]
}

# Find the position (index) of the n-th top-most max within a vector of numbers
which.nth.max <- function(x, n)
{
  which(x == nth.max(x,n))
}

# Pad a matrix with zeros to reach a desired dimension
pad_matrix <- function(M, nrow_target, ncol_target)
{
  out <- matrix(0, nrow_target, ncol_target)
  out[seq_len(nrow(M)), seq_len(ncol(M))] <- M
  out
}

# Find the n top-most relative maxima (i.e., peaks) within a signal (numeric
# vector), in terms of both their values (first output element) and
# positions/indexes (second output element).
nth.peaks <- function(x, n, min_h = 5)
{
  findpeaks(x,
            minpeakdistance = 10, # minimum number of samples between peaks
            minpeakheight = min_h, # ignore peaks below a threshold
            threshold = 0.5 # required drop from the peak to neighbors
  ) -> peak_list
  
  if(is.null(peak_list)) {
    peak_list <- matrix(0,1,4) # Standard 'findpeaks' output has 4 columns
  }
  if (nrow(peak_list) < n) {
    peak_list |> pad_matrix(n,4) -> peak_list
  }
  
  peak_list[order(peak_list[,1], decreasing = TRUE)[1:n], ]
}

# --- Feature extraction -------------------------------------------------------

# Detects signal collapse in a single trace (vector 'sig') by convolution of a
# downward step. By using such a step as a convolution window, the signal
# drop-off is mapped to a (negative) peak, while possible peaks in the original
# signal are mapped to biphasic spikes, just as in the case of a standard
# derivative operator... indeed:
#
#     convolve(x, rev(c(1, -1)), type = "filter")
#
# is exactly equivalent to
#
#     diff(x)
#
# however, a wider convolution window helps make the detection more robust.
# Based on this, after convolution, it will be sufficient to apply a negative
# threshold (thr) to identify the collapses, and a "proximity filter" to
# effectively distinguish such drop-off transients from spikes.
step_convolve <- function(sig,
                          median_width = 5, # odd number here, please!
                          step_win = c(1, 1, 1, -1, -1, -1))
{
  win <- step_win / sum(abs(step_win)) # Normalize convolution window
  
  # Pre-filter by a rolling median to remove possible single-time-point spikes,
  # "dropping points", or brief transients, likely due to artifacts from camera
  # sensors or alike, in any case not representing the true collapse events we
  # are interested in (robustness). Median removes spikes, without altering or
  # smoothing long-lasting signal components.
  sig |> rollapply(width = median_width,
                   FUN = median,
                   align = "center") |> convolve(rev(win),
                                                 type = "filter")
}

# NOTE: As an additional safeguard, the initial portion of the signal is "muted"
# and is not included in the detection of collapse events. This prevents the
# sporadic activity typically observed at the beginning of recordings (in the
# form of rapid, multiple spikes) from invalidating a subsequent drop-off which,
# although real, may produce a smaller deflection in the derivative plot.
extract <- function(sig,
                    median_width = 5, # odd number here, please!
                    step_win = c(1, 1, 1, -1, -1, -1),
                    protect = 0, # to exclude the initial part of the signal (up to index 'protect')
                    thr = 5)
{
  sig |> step_convolve(median_width, step_win) |>
    {\(x) c(rep(0,protect), x[(protect+1):length(x)])}() |> # cancel the protected part
    {\(x) c(min(x),which.min(x),nth.peaks(x,3,min_h=thr)[,2])}() -> biggest_drop
  
  # Both rollapply(...) and convolve(..., type = "filter") return a
  # (win-1)-trimmed signal, requiring correction.
  correct <- function(x){ceiling((x-1)/2)}
  biggest_drop[2] <- biggest_drop[2] + correct(median_width) + correct(length(step_win))
  
  # Exclude points above the negative threshold and those points that have a
  # top-most maximum within the range 'k*length(step_win)' ("proximity filter").
  k <- 3
  ifelse(biggest_drop[1] < -thr &&
           abs(biggest_drop[2]-biggest_drop[3]) > k*length(step_win) &&
           abs(biggest_drop[2]-biggest_drop[4]) > k*length(step_win) &&
           abs(biggest_drop[2]-biggest_drop[5]) > k*length(step_win),
         biggest_drop[2], NA_integer_)
}

# AUC of a single normalized trace (vector 'sig') up to collapse (if found) or
# full duration otherwise.
auc <- function(sig, collapse_idx = length(sig) + 1, time_vec)
{
  # Take only the positive parts of the signal (no negative AUCs)
  sig[sig < 0] <- 0
  # Set the endpoint
  end_idx <- ifelse(is.na(collapse_idx), length(sig), collapse_idx - 1)
  if (end_idx > 2) {
    trapz(time_vec[1:end_idx], sig[1:end_idx])
  } else {
    NA_real_
  }
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

# Find the first time-point that crosses a reference/target value (upward) using
# linear interpolation.
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

# --- per-Experiment stats -----------------------------------------------------

# Functions to be called by 'get_stat()'. Dots (...) are important for the
# 'na.rm = TRUE' argument used within 'sapply()'.
sem <- function(x, ...) {sd(x, ...) / sqrt(sum(!is.na(x)))}
get_size <- function(x, ...) {sum(!is.na(x))}
get_nas <- function(x, ...) {sum(is.na(x))}

# General function to compute statistics from the 'summary_tbl'
get_stat <- function(df, FUN, ...)
{
  FUN <- match.fun(FUN) # safer pattern
  
  if (identical(FUN, mean)) {
    # Compute the collapsing rate as percentage of drop-offs per experiment
    df$collapse_idx |> is.na() |> {\(x)!x}() |> sum() -> not_NAs
    rate <- (not_NAs / nrow(df)) * 100
  } else {
    rate <- NA
  }
  names(rate) <- "collapse_rate"
  
  # Compute the FUN statistics for all the vars
  df |> sapply(is.numeric) -> idx
  df[,idx] |> sapply(FUN, na.rm = TRUE) |> append(rate, after = 3)
}

# --- Plot Comparisons ---------------------------------------------------------

# Map p-values to asterisks (common convention)
p_to_stars <- function(p, symb = "*")
{
  if (is.na(p))   return(NA_character_)
  if (p <= 0.001) return(paste(rep(symb,3), collapse = ""))
  if (p <= 0.01)  return(paste(rep(symb,2), collapse = ""))
  if (p <= 0.05)  return(paste(rep(symb,1), collapse = ""))
  return("")  # no asterisk if not significant at 0.05
}

# Is the 'p_to_stars' output from a significant comparison
is.sig <- function(stars) { !is.na(stars) && stars != "" }

# Plot boxplots with jittered points of two samples, also  displaying
# significance as the appropriate number of asterisks above the two boxes
ttest_boxplot <- function(x,
                          y,
                          xlab = "",
                          ylab = "",
                          group_labels = c("Group1", "Group2"),
                          t.pval = NA_real_,
                          U.pval = NA_real_,
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
  # Remove NAs and NaNs
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) == 0 || length(y) == 0) {
    message("After removing NAs, one group is empty.")
    return(ggplot() + theme_void() +
             annotate("text", x = 0, y = 0, label = "No data available"))
  }
  
  # p-values and stars
  pvals <- c(t = unname(t.pval),
             U = unname(U.pval))
  stars <- c(t = p_to_stars(t.pval),
             U = p_to_stars(U.pval, "+"))
  
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
  y_t.stars <- y_line + h * 0.50
  y_U.stars <- y_line + h * 0.25
  y_t.text <- y_line - h * 0.25
  y_U.text <- y_line - h * 0.75
  
  # Build the plot
  plt <- ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(width = box.width, linewidth = 0.7, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = jitter.width, height = 0, size = point.size, alpha = 0.7, aes(color = group)) +
    scale_fill_manual(values = palette, guide = "none") +
    scale_color_manual(values = palette, guide = "none") +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 14)
  
  # Add significance line and stars (only if we have a star string to show)
  if (is.sig(stars["t"]) || is.sig(stars["U"])) {
    plt <- plt +
      geom_segment(aes(x = 1, xend = 2, y = y_line, yend = y_line), inherit.aes = FALSE, size = 1) +
      geom_segment(aes(x = 1, xend = 1, y = y_line, yend = y_line - h * 0.08), inherit.aes = FALSE, size = 1) +
      geom_segment(aes(x = 2, xend = 2, y = y_line, yend = y_line - h * 0.08), inherit.aes = FALSE, size = 1) +
      annotate("text", x = 1.5, y = y_t.stars, label = stars["t"], size = star.size, fontface = "bold") +
      annotate("text", x = 1.5, y = y_U.stars, label = stars["U"], size = star.size, fontface = "bold")
  }
  
  # Optionally show actual p-values (even if not significant)
  if (show.p) {
    for (test in names(pvals)) {
      if (!is.na(pvals[test])) {
        p_label <- paste0("p[", test, "] = ", formatC(pvals[test], digits = p_digits, format = "g"))
        if (stars[test] == "") {
          plt <- plt +
            geom_segment(aes(x = 1, xend = 2, y = y_line, yend = y_line), inherit.aes = FALSE, size = 0.4, linetype = "dashed")
        }
        plt <- plt +
          annotate("text", x = 1.5, y = ifelse(test == "t", y_t.text, y_U.text), label = p_label, size = 4)
      } else {
        plt <- plt +
          annotate("text", x = 1.5, y = ifelse(test == "t", y_t.text, y_U.text), label = "p = NA", size = 4)
      }
    }
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
