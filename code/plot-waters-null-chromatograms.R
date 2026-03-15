library(tidyverse)
library(truncnorm)
library(viridis)
library(grid)
library(gridExtra)
library(jsonlite)

# Clear the environment
rm(list = ls())

# Set seed for reproducibility
set.seed(42)

# Function to create Gaussians for chemical noise
generate_gaussian <- function(mu, amplitude, scale, time) {
  # Step size in the time series
  delta_time <- mean(diff(time)) # Average spacing between time points

  # Sigma in data points and time units
  sigma_points <- scale # Gaussian width in data points
  sigma_time <- sigma_points * delta_time # Gaussian width in time units

  # Adjust mu to the closest time point in the time series
  mu_closest <- time[which.min(abs(time - mu))]

  # Generate the Gaussian
  gaussian <- amplitude * exp(-((time - mu_closest)^2) / (2 * sigma_time^2))
  return(gaussian)
}

# Sampling function for extracting values from power density
sample_from_kde <- function(density, n) {
  sample(density$x, size = n, prob = density$y, replace = TRUE)
}

# Generate a null hypothesis chromatogram
generate_null_chromatogram <- function(retention_time, peak_data, power_density,
                                       lambda_events, num_new_events,
                                       blank_noise_mean, blank_noise_std) {
  n_points <- length(retention_time)

  # Draw number of existing chemical noise events from Poisson
  # Cap at pool size to prevent sample() error (cap essentially never triggers)
  num_existing <- min(rpois(1, lambda_events), nrow(peak_data))

  # Sample existing events
  if (num_existing > 0) {
    sampled_times <- sample(peak_data$time, num_existing)
    sampled_scales <- sample(peak_data$scale, num_existing)
    sampled_powers <- sample_from_kde(power_density, num_existing)
  } else {
    sampled_times <- numeric(0)
    sampled_scales <- numeric(0)
    sampled_powers <- numeric(0)
  }

  # Generate new events
  new_times <- runif(num_new_events, min(retention_time), max(retention_time))
  new_scales <- sample(peak_data$scale, num_new_events, replace = TRUE)
  new_powers <- sample_from_kde(power_density, num_new_events)

  # Initialize arrays for the chromatogram and components
  chemical_noise <- numeric(length(retention_time))
  sampled_component <- numeric(length(retention_time))
  new_component <- numeric(length(retention_time))

  # Generate sampled peaks (first color)
  for (i in seq_along(sampled_times)) {
    # Amplitude is sqrt of wavelet power from empirical data
    amplitude <- sqrt(sampled_powers[i])

    gaussian <- generate_gaussian(
      sampled_times[i],
      amplitude,
      sampled_scales[i],
      retention_time
    )
    sampled_component <- sampled_component + gaussian
    chemical_noise <- chemical_noise + gaussian
  }

  # Generate new peaks (second color)
  for (i in seq_along(new_times)) {
    # Amplitude is sqrt of wavelet power from empirical data
    amplitude <- sqrt(new_powers[i])

    gaussian <- generate_gaussian(
      new_times[i],
      amplitude,
      new_scales[i],
      retention_time
    )
    new_component <- new_component + gaussian
    chemical_noise <- chemical_noise + gaussian
  }

  # Generate blank noise (electronic noise) from empirical blank distribution
  electronic_noise <- rtruncnorm(n_points,
    a = 0, b = Inf,
    mean = blank_noise_mean,
    sd = blank_noise_std
  )

  # Noisy signal = chemical peaks + blank noise
  total_signal <- chemical_noise + electronic_noise

  # Return all components
  return(list(
    time = retention_time,
    total = total_signal,
    sampled = sampled_component,
    new = new_component,
    electronic = electronic_noise,
    sampled_times = sampled_times,
    sampled_scales = sampled_scales,
    sampled_powers = sampled_powers,
    new_times = new_times,
    new_scales = new_scales,
    new_powers = new_powers
  ))
}

# Plot a noise-free chromatogram
plot_noise_free_chromatogram <- function(chromatogram, title = "Noise-Free Chromatogram") {
  # Create a data frame for plotting the noise-free chromatogram
  # This is just the chemical noise (sampled + new components) without electronic noise
  df <- tibble(
    Time = chromatogram$time,
    Total = chromatogram$sampled + chromatogram$new # Chemical noise only (no electronic noise)
  )

  # Create data frames for the peak locations
  sampled_peaks_df <- tibble(
    Time = chromatogram$sampled_times,
    Type = "Sampled"
  )

  new_peaks_df <- tibble(
    Time = chromatogram$new_times,
    Type = "New"
  )

  # Create the plot with vertical lines for peaks
  p <- ggplot() +
    # Black line for noise-free chromatogram
    geom_line(data = df, aes(x = Time, y = Total), color = "black", linewidth = 0.5) +
    # Blue vertical dashed lines for sampled peaks
    geom_vline(
      data = sampled_peaks_df, aes(xintercept = Time),
      color = "blue", linetype = "dashed", alpha = 0.7
    ) +
    # Red vertical dashed lines for new peaks
    geom_vline(
      data = new_peaks_df, aes(xintercept = Time),
      color = "red", linetype = "dashed", alpha = 0.7
    ) +
    # Add labels
    labs(
      title = title,
      subtitle = paste(
        "Noise-Free Version | Blue: Sampled peaks (", nrow(sampled_peaks_df), ")",
        " | Red: New peaks (", nrow(new_peaks_df), ")",
        sep = ""
      ),
      x = "Retention Time",
      y = "Intensity"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    )

  return(p)
}

# Plot a noisy chromatogram
plot_noisy_chromatogram <- function(chromatogram, title = "Noisy Chromatogram") {
  # Create a data frame for plotting the total chromatogram (with noise)
  df <- tibble(
    Time = chromatogram$time,
    Total = chromatogram$total # Total signal including electronic noise
  )

  # Create data frames for the peak locations
  sampled_peaks_df <- tibble(
    Time = chromatogram$sampled_times,
    Type = "Sampled"
  )

  new_peaks_df <- tibble(
    Time = chromatogram$new_times,
    Type = "New"
  )

  # Create the plot with vertical lines for peaks
  p <- ggplot() +
    # Black line for total chromatogram
    geom_line(data = df, aes(x = Time, y = Total), color = "black", linewidth = 0.5) +
    # Blue vertical dashed lines for sampled peaks
    geom_vline(
      data = sampled_peaks_df, aes(xintercept = Time),
      color = "blue", linetype = "dashed", alpha = 0.7
    ) +
    # Red vertical dashed lines for new peaks
    geom_vline(
      data = new_peaks_df, aes(xintercept = Time),
      color = "red", linetype = "dashed", alpha = 0.7
    ) +
    # Add labels
    labs(
      title = title,
      subtitle = paste(
        "With Electronic Noise | Blue: Sampled peaks (", nrow(sampled_peaks_df), ")",
        " | Red: New peaks (", nrow(new_peaks_df), ")",
        sep = ""
      ),
      x = "Retention Time",
      y = "Intensity"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    )

  return(p)
}

# -----------------------------------------------------------------------------
# 1. Load Lorazepam quant results and chemical noise statistics (Waters only)
# -----------------------------------------------------------------------------

# Load pre-computed Lorazepam quantitative wavelet results from JSON
simulation_results <- fromJSON(
  file.path("results", "json", "Lorazepam_quant_results.json"),
  simplifyVector = FALSE
)

cat("Generating null hypothesis chromatograms for Waters instrument\n")
cat("Found", length(simulation_results), "samples\n")

# Use the first sample's time grid (all samples from same instrument share identical grid)
retention_time <- unlist(simulation_results[[1]]$time)

# Load chemical noise statistics for Lorazepam quant (Waters)
peak_data <- read_csv(
  file.path("results", "tables", "Lorazepam_quant_peak_data_waters.csv")
)
power_density <- density(peak_data$power, from = 0)
power_density$y <- power_density$y / sum(power_density$y)

# Compute lambda for Poisson event count from per-sample peak counts
events_per_sample <- peak_data |>
  group_by(Sample) |>
  summarise(n_events = n(), .groups = "drop")
lambda_events <- mean(events_per_sample$n_events)

cat("  Waters: lambda_events =", round(lambda_events, 2),
    "(from", nrow(events_per_sample), "samples)\n")

# Electronic noise from Waters reagent blank (ROS estimates from pooled blanks)
noise_pooled <- read_csv(file.path("results", "tables", "waters-noise-pooled.csv"))
lorazepam_quant_noise <- noise_pooled |>
  filter(Compound == "Lorazepam", Chrom == "Quant")
blank_noise_mean <- lorazepam_quant_noise$ROS_Mean
blank_noise_std <- lorazepam_quant_noise$ROS_SD
num_new_events <- 1

# -----------------------------------------------------------------------------
# 2. Generate null hypothesis chromatograms
# -----------------------------------------------------------------------------

set.seed(123)
chromatograms <- list()
for (i in 1:10) {
  chromatograms[[i]] <- generate_null_chromatogram(
    retention_time = retention_time,
    peak_data = peak_data,
    power_density = power_density,
    lambda_events = lambda_events,
    num_new_events = num_new_events,
    blank_noise_mean = blank_noise_mean,
    blank_noise_std = blank_noise_std
  )
}

# -----------------------------------------------------------------------------
# 3. Create plots (5 pairs of noise-free and noisy)
# -----------------------------------------------------------------------------

noise_free_plots <- list()
noisy_plots <- list()

for (i in 1:5) {
  noise_free_plots[[i]] <- plot_noise_free_chromatogram(chromatograms[[i]],
    title = paste("Chromatogram", i)
  )
  noisy_plots[[i]] <- plot_noisy_chromatogram(chromatograms[[i]],
    title = paste("Chromatogram", i)
  )
}

# Arrange plots in pairs (noise-free in column 1, noisy in column 2)
paired_plots <- list()
for (i in 1:5) {
  paired_plots[[2 * i - 1]] <- noise_free_plots[[i]]
  paired_plots[[2 * i]] <- noisy_plots[[i]]
}

# Create title for the Lorazepam dilution series (Waters instrument)
title_expr <- bquote(
  "Lorazepam Quant (m/z 321.2" %->% " 275): Null Hypothesis Chromatograms (Waters)"
)

# Create a PDF with 5 rows and 2 columns with overall title
pdf(file.path("figures", "main_figures", "waters-null-chromatograms.pdf"),
    width = 11, height = 8.5)
grid.arrange(
  grobs = paired_plots,
  ncol = 2,
  nrow = 5,
  top = textGrob(title_expr, gp = gpar(fontsize = 14, fontface = "bold"))
)
dev.off()

# Print a message to confirm completion
cat(
  "\nGenerated 10 representative null hypothesis chromatograms for",
  "Waters instrument and saved to PDF file.\n"
)
cat("File created: figures/main_figures/waters-null-chromatograms.pdf\n")
