# Waters Chemical Noise Comparison Plots
# Creates a single 2x3 figure: time histograms (left) and power distributions
# (right) for Lorazepam, Clonidine, and Gabapentin (all Waters, Quant transitions)

# Load required packages
library(jsonlite)
library(tidyverse)
library(patchwork)

# Clear environment
rm(list = ls())

# Read JSON files with simplifyVector = FALSE to preserve list structure
lorazepam_data <- fromJSON(file.path("results", "json",
                                     "Lorazepam_quant_chemical_noise.json"),
                           simplifyVector = FALSE)
clonidine_data <- fromJSON(file.path("results", "json",
                                     "Clonidine_quant_chemical_noise.json"),
                           simplifyVector = FALSE)
gabapentin_data <- fromJSON(file.path("results", "json",
                                      "Gabapentin_chemical_noise.json"),
                            simplifyVector = FALSE)

# Extract chemical noise peaks from all entries (single instrument per file)
extract_peaks <- function(json_data) {
  map_dfr(json_data, function(result) {
    peaks <- result$detected_peaks
    if (is.null(peaks) || length(peaks) == 0) {
      return(tibble(time = numeric(0), power = numeric(0)))
    }
    bind_rows(lapply(peaks, function(p) {
      if (!is.null(p$time) && !is.null(p$power)) {
        tibble(time = p$time, power = p$power)
      } else {
        tibble(time = numeric(0), power = numeric(0))
      }
    }))
  })
}

# Extract data for each compound
lorazepam_peaks  <- extract_peaks(lorazepam_data)
clonidine_peaks  <- extract_peaks(clonidine_data)
gabapentin_peaks <- extract_peaks(gabapentin_data)

# Create time histogram function with fixed bins
create_time_histogram <- function(data, title) {
  if (nrow(data) == 0) {
    return(ggplot() +
           ggtitle(title) +
           theme_minimal() +
           annotate("text", x = 0.5, y = 0.5, label = "No data"))
  }

  ggplot(data, aes(x = time)) +
    geom_histogram(bins = 25, fill = "blue", color = "black", alpha = 0.7) +
    labs(x = "Time (min)", y = "Frequency") +
    ggtitle(title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    )
}

# Create power distribution function with fixed bins
create_power_distribution <- function(data, title) {
  if (nrow(data) == 0) {
    return(ggplot() +
           ggtitle(title) +
           theme_minimal() +
           annotate("text", x = 0.5, y = 0.5, label = "No data"))
  }

  # Convert to dB scale
  data_db <- data |>
    mutate(power_db = 10 * log10(power + 1e-10)) |>
    filter(is.finite(power_db))

  ggplot(data_db, aes(x = power_db)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 15,
                   fill = "lightblue",
                   color = "black",
                   alpha = 0.5) +
    geom_density(color = "darkred", linewidth = 1) +
    labs(x = "Power (dB)", y = "Density") +
    ggtitle(title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    )
}

# Create time histograms (left column: a, c, e)
time_a <- create_time_histogram(lorazepam_peaks,
  expression("(a) Lorazepam Quant (m/z 321.2" %->% "275)"))
time_c <- create_time_histogram(clonidine_peaks,
  expression("(c) Clonidine Quant (m/z 232.2" %->% "215.1)"))
time_e <- create_time_histogram(gabapentin_peaks,
  expression("(e) Gabapentin Quant (m/z 172.3" %->% "154.3)"))

# Create power distributions (right column: b, d, f)
power_b <- create_power_distribution(lorazepam_peaks,
  expression("(b) Lorazepam Quant (m/z 321.2" %->% "275)"))
power_d <- create_power_distribution(clonidine_peaks,
  expression("(d) Clonidine Quant (m/z 232.2" %->% "215.1)"))
power_f <- create_power_distribution(gabapentin_peaks,
  expression("(f) Gabapentin Quant (m/z 172.3" %->% "154.3)"))

# Combine into 2x3 layout: time histograms left, power distributions right
combined_plot <- (time_a | power_b) / (time_c | power_d) / (time_e | power_f)

# Save combined figure
pdf(file.path("figures", "main_figures",
              "waters-chemical-noise-combined.pdf"),
    width = 11, height = 10)
print(combined_plot)
dev.off()

# Save standalone time histograms (1 column, 3 rows)
time_plot <- time_a / time_c / time_e
pdf(file.path("figures", "main_figures",
              "waters-chemical-noise-time-histograms.pdf"),
    width = 6, height = 10)
print(time_plot)
dev.off()

# Save standalone power distributions (1 column, 3 rows)
power_plot <- power_b / power_d / power_f
pdf(file.path("figures", "main_figures",
              "waters-chemical-noise-power-distributions.pdf"),
    width = 6, height = 10)
print(power_plot)
dev.off()

cat("Plots saved to:\n")
cat("  - figures/main_figures/waters-chemical-noise-combined.pdf\n")
cat("  - figures/main_figures/waters-chemical-noise-time-histograms.pdf\n")
cat("  - figures/main_figures/waters-chemical-noise-power-distributions.pdf\n")
