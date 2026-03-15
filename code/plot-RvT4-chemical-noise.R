# RvT4 Chemical Noise Comparison Plots
# Creates 1x2 panel plots comparing chemical noise across SCIEX instruments
# Layout: 6500 top, 7500 bottom

# Load required packages
library(jsonlite)
library(tidyverse)
library(patchwork)

# Clear environment
rm(list = ls())

# Read JSON file with simplifyVector = FALSE to preserve list structure
rvt4_data <- fromJSON(file.path("results", "json",
                                "RvT4_chemical_noise.json"),
                      simplifyVector = FALSE)

# Extract chemical noise peaks for a specific instrument
extract_peaks <- function(json_data, instrument_filter) {
  all_peaks <- map_dfr(json_data, function(result) {
    if (result$instrument == instrument_filter) {
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
    } else {
      tibble(time = numeric(0), power = numeric(0))
    }
  })

  return(all_peaks)
}

# Extract data for each instrument
rvt4_6500 <- extract_peaks(rvt4_data, "6500")
rvt4_7500 <- extract_peaks(rvt4_data, "7500")

# Create time histogram function with fixed bins
create_time_histogram <- function(data, title) {
  if (nrow(data) == 0) {
    return(ggplot() +
           ggtitle(title) +
           theme_minimal() +
           annotate("text", x = 0.5, y = 0.5, label = "No data"))
  }

  ggplot(data, aes(x = time)) +
    geom_histogram(bins = 25, fill = "blue",
                   color = "black", alpha = 0.7) +
    labs(x = "Time (min)", y = "Frequency") +
    ggtitle(title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10,
                                face = "bold"),
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
      plot.title = element_text(hjust = 0.5, size = 10,
                                face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    )
}

# Create time histograms for both instruments
time_a <- create_time_histogram(rvt4_6500,
  expression("(a) RvT4 Quant (m/z 361.1" %->% "211.1) 6500"))
time_b <- create_time_histogram(rvt4_7500,
  expression("(b) RvT4 Quant (m/z 361.1" %->% "211.1) 7500"))

# Create power distributions for both instruments
power_a <- create_power_distribution(rvt4_6500,
  expression("(a) RvT4 Quant (m/z 361.1" %->% "211.1) 6500"))
power_b <- create_power_distribution(rvt4_7500,
  expression("(b) RvT4 Quant (m/z 361.1" %->% "211.1) 7500"))

# Combine into 1-column, 2-row layouts
time_plot <- time_a / time_b
power_plot <- power_a / power_b

# Save time histogram plot
pdf(file.path("figures", "main_figures",
              "RvT4-chemical-noise-time-histograms.pdf"),
    width = 6, height = 8)
print(time_plot)
dev.off()

# Save power distribution plot
pdf(file.path("figures", "main_figures",
              "RvT4-chemical-noise-power-distributions.pdf"),
    width = 6, height = 8)
print(power_plot)
dev.off()

cat("Plots saved to:\n")
cat("  - figures/main_figures/RvT4-chemical-noise-time-histograms.pdf\n")
cat("  - figures/main_figures/RvT4-chemical-noise-power-distributions.pdf\n")
