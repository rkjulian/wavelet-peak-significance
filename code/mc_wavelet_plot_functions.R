# Shared plotting functions for RvT4 wavelet analysis scripts
#
# Sourced by all RvT4 plotting and analysis scripts.

# Function to create formatted labels for intensity in plots
inten_label <- function(break_value) {
  expressions <- vector("list", length(break_value))

  for (i in seq_along(break_value)) {
    if (break_value[i] == 0 || is.na(break_value[i])) {
      expressions[[i]] <- parse(text = "0")[[1]]
    } else {
      mantissa <- break_value[i] / 10^floor(log10(abs(break_value[i])))
      exponent <- floor(log10(abs(break_value[i])))
      text <- sprintf("%.2f*x*10^%d", mantissa, as.integer(exponent))
      expressions[[i]] <- parse(text = text)[[1]]
    }
  }

  as.expression(expressions)
}

# Function to sanitize sample names for filenames
sanitize_filename <- function(name) {
  # Replace problematic characters with underscores
  name <- gsub("[/\\:*?\"<>|]", "_", name)
  # Replace spaces with underscores
  name <- gsub(" ", "_", name)
  # Remove leading/trailing underscores
  name <- gsub("^_+|_+$", "", name)
  return(name)
}

# Function to print summary plots for chemical noise characterization
#
# Parameters:
#   instrument     - instrument identifier (e.g. "6500", "7500")
#   peak_data      - named list of peak data frames (one per sample)
#   csv_prefix     - filename prefix for CSV output (e.g. "RvT4_peak_data_")
#   compound_name  - compound label for titles (e.g. "Quant" or "Qual")
#   precursor_mz   - precursor m/z value for titles (e.g. "361.1")
#   product_mz     - product m/z value for titles (e.g. "211.1" or "193.1")
print_summary_plots <- function(instrument, peak_data, csv_prefix,
                                compound_name, precursor_mz, product_mz) {
  # Combine all samples' data into a single data frame
  peak_data_combined <- bind_rows(peak_data, .id = "Sample")

  # Save complete data to CSV
  peak_data_combined |>
    select(-prominence) |>
    write_csv(file.path("results", "tables", paste0(csv_prefix, instrument, ".csv")))

  # Filter out non-finite values for plotting
  peak_data_valid <- peak_data_combined |>
    filter(is.finite(time), is.finite(power))

  # Build dynamic title components
  time_title <- bquote(
    "Time Distribution - RvT4" ~ .(compound_name) ~ "(m/z" ~ .(precursor_mz) %->% .(product_mz) ~ ") -" ~ .(instrument)
  )
  power_title <- bquote(
    "Peak Power Distribution - RvT4" ~ .(compound_name) ~ "(m/z" ~ .(precursor_mz) %->% .(product_mz) ~ ") -" ~ .(instrument)
  )

  # Create time distribution plot
  plot3 <- ggplot(peak_data_valid, aes(x = time)) +
    geom_histogram(bins = 25, fill = "blue", color = "black", alpha = 0.7) +
    labs(
      title = time_title,
      x = "Time (min)", y = "Frequency"
    ) +
    theme_minimal()

  # Create power distribution plot
  plot4 <- peak_data_valid |>
    mutate(power_db = 10 * log10(power + 1e-10)) |>
    filter(is.finite(power_db)) |>
    ggplot(aes(x = power_db)) +
    geom_histogram(aes(y = after_stat(density)),
      bins = 15,
      fill = "lightblue",
      color = "black",
      alpha = 0.5
    ) +
    geom_density(
      color = "darkred",
      linewidth = 1
    ) +
    labs(
      title = power_title,
      x = "Power (dB)", y = "Density"
    ) +
    theme_minimal()

  # Arrange plots vertically
  dist_plots <- plot3 / plot4
  print(dist_plots)
}

# Plot a 3-panel null wavelet analysis result for a single sample.
# Panel 1: Chromatogram with expected RT
# Panel 2: Wavelet power spectrum (asinh-scaled)
# Panel 3: Significance testing (-ln(p) with FWER correction)
#
# Caller is responsible for opening a PDF device and setting par() layout.
#
# Parameters:
#   retention_time - numeric vector of time points
#   sample_y       - numeric vector of intensity values
#   adjusted_power - matrix of wavelet power [n_times x n_scales]
#   scales         - numeric vector of wavelet scales
#   p_adjusted     - numeric vector of FWER-adjusted p-values per time point
#   expected_rt    - expected retention time for reference line
#   sample_name    - sample label for title
#   alpha          - significance level (default 0.05)
#   n_sim          - number of MC simulations (for uncertainty band)
#   compound_name  - compound label for title (e.g. "Quant" or "Qual")
#   precursor_mz   - precursor m/z for title (e.g. "361.1")
#   product_mz     - product m/z for title (e.g. "211.1" or "193.1")
#   plot_coi       - whether to plot cone of influence shading (default FALSE)
plot_null_wavelet_sample <- function(retention_time, sample_y, adjusted_power,
                                     scales, p_adjusted, expected_rt,
                                     sample_name, alpha, n_sim,
                                     compound_name, precursor_mz, product_mz,
                                     plot_coi = FALSE, analyte_name = "RvT4") {
  n_points <- length(retention_time)
  significant_times <- which(p_adjusted < alpha)

  # Plot 1: Chromatogram with expected RT
  plot(retention_time, sample_y,
    type = "l", col = "blue",
    main = bquote(.(analyte_name) ~ .(compound_name) ~ "(m/z" ~ .(precursor_mz) %->% .(product_mz) ~ "):" ~ .(sample_name)),
    xlab = "Time", ylab = "Amplitude",
    xaxs = "i"
  )
  abline(v = expected_rt, col = "black", lty = 2)

  # Plot 2: Wavelet power spectrum (asinh-scaled)
  power_transformed <- asinh(adjusted_power)

  image(
    x = retention_time,
    y = scales,
    z = power_transformed,
    col = viridis(256),
    xlab = "Time",
    ylab = "Scale",
    main = "Mexican Hat Wavelet Power Spectrum (asinh scaled)"
  )

  # Plot 3: Significance testing
  y_limits <- range(-log(p_adjusted[is.finite(-log(p_adjusted))]))
  y_range <- diff(y_limits)
  y_min <- y_limits[1] - 0.05 * y_range
  y_max <- y_limits[2] + 0.05 * y_range

  plot(retention_time, -log(p_adjusted),
    type = "l",
    main = "Significance Testing (Holm FWER-adjusted minimum p-values)",
    xlab = "Time", ylab = "-ln(p)",
    xaxs = "i", ylim = c(y_min, max(y_max, -log(alpha)))
  )

  # Cone of influence shading
  if (plot_coi) {
    dt <- mean(diff(retention_time))
    max_scale <- max(scales)
    coi_width_time <- 5 * max_scale * dt

    rect(retention_time[1], y_min, retention_time[1] + coi_width_time, y_max,
      col = rgb(0.5, 0.5, 0.5, 0.3), border = NA
    )
    rect(retention_time[n_points] - coi_width_time, y_min, retention_time[n_points], y_max,
      col = rgb(0.5, 0.5, 0.5, 0.3), border = NA
    )
  }

  # MC uncertainty band around significance threshold
  mc_se_p <- sqrt(alpha * (1 - alpha) / n_sim)
  mc_se_neglogp <- mc_se_p / alpha
  rect(retention_time[1], -log(alpha) - 2 * mc_se_neglogp,
    retention_time[n_points], -log(alpha) + 2 * mc_se_neglogp,
    col = rgb(1, 0, 0, 0.15), border = NA
  )

  abline(h = -log(alpha), col = "red", lty = 2)
  points(retention_time[significant_times],
    -log(p_adjusted[significant_times]),
    col = "red", pch = 20
  )
  abline(v = expected_rt, col = "black", lty = 2)
}
