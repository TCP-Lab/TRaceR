


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











