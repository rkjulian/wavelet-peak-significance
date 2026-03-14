# Load required packages
library(MSnbase) # For reading mzML files
library(tidyverse) # Data manipulation and visualization
library(MassSpecWavelet) # Wavelet analysis for mass spectrometry
library(patchwork) # Combining ggplot2 plots
library(ggh4x) # Enhanced ggplot2 facets
library(forcats) # Working with factors
library(jsonlite) # For JSON export


# Clear the environment
rm(list = ls())

# Source shared functions
source(file.path("code", "mc_wavelet_functions.R"))
source(file.path("code", "mc_wavelet_plot_functions.R"))

# Set seed for reproducibility
set.seed(42)

# Define instrument-specific data
instrument_data <- list(
  "6500" = list(
    file_list = c(
      "MW-RvT4ex3-plasma-MS4-MF-std mix.mzML",
      "MW-RvT4ex3-plasma-MS4-MF-1-plasma WTchow.mzML",
      "MW-RvT4ex3-plasma-MS4-MF-2-plasma WTchow.mzML",
      "MW-RvT4ex3-plasma-MS4-MF-3-plasma WTchow.mzML",
      "MW-RvT4ex3-plasma-MS4-MF-4-plasma WTchow.mzML",
      "MW-RvT4ex3-plasma-MS4-MF-5-plasma WTWD.mzML",
      "MW-RvT4ex3-plasma-MS4-MF-6-plasma WTWD.mzML",
      "MW-RvT4ex3-plasma-MS4-MF-7-plasma WTWD.mzML",
      "MW-RvT4ex3-plasma-MS4-MF-8-plasma WTWD.mzML"
    ),
    sample_list = c(
      "STD",
      "Chow 1", "Chow 2", "Chow 3", "Chow 4",
      "WD 1", "WD 2", "WD 3", "WD 4"
    ),
    quant_index = 117,
    qual_index = 116,
    is_index = 33,
    is_rt = 0
  ),
  "7500" = list(
    file_list = c(
      "JD_ACP_ApoE_RvT_052022-Std mix.mzML",
      "JD_ACP_ApoE_RvT_052022-Sample 1.mzML",
      "JD_ACP_ApoE_RvT_052022-Sample 2.mzML",
      "JD_ACP_ApoE_RvT_052022-Sample 3.mzML",
      "JD_ACP_ApoE_RvT_052022-Sample 4.mzML",
      "JD_ACP_ApoE_RvT_052022-Sample 5.mzML",
      "JD_ACP_ApoE_RvT_052022-Sample 6.mzML"
    ),
    sample_list = c(
      "STD",
      "Chow 1", "Chow 2", "Chow 3",
      "WD 1", "WD 2", "WD 3"
    ),
    quant_index = c(99, 34, 34, 34, 34, 34, 34),
    qual_index = c(98, 33, 33, 33, 33, 33, 34),
    is_index = c(32, 17, 17, 17, 17, 17, 17),
    is_rt = 0
  )
)

# Cross-scale consistency threshold for chemical noise peak detection.
# A peak must be detected at >= this many wavelet scales to be considered a
# real chromatographic feature rather than electronic noise. Electronic noise
# produces incoherent fluctuations that appear at only 1-2 scales, while real
# peaks produce coherent wavelet responses across multiple scales.
# With 31 scales (1.5 to 3.0 by 0.05), 5 scales = ~16% of the scale range.
min_scales_detected <- 5

# Create combined data frame
combined_data <- map_dfr(names(instrument_data), ~ create_instrument_data(
  .x,
  instrument_data[[.x]]$file_list,
  instrument_data[[.x]]$sample_list,
  instrument_data[[.x]]$quant_index,
  instrument_data[[.x]]$qual_index,
  instrument_data[[.x]]$is_index,
  instrument_data[[.x]]$is_rt
))

# Process data for both instruments
data_dir <- "data/S-BSST880"
chromatogram_data <- empty_chrom()

for (i in 1:nrow(combined_data)) {
  row <- combined_data[i, ]
  results <- process_srm_data(row, data_dir)

  # Add the instrument column to the results
  results$chromatogram_data$instrument <- row$instrument

  chromatogram_data <- bind_rows(chromatogram_data, results$chromatogram_data)
  combined_data$is_rt[i] <- results$is_rt
}

# Calculate quant offsets
quant_offset_6500 <- calculate_quant_offset(combined_data, "6500", chromatogram_data, "Qual", "STD")
quant_offset_7500 <- calculate_quant_offset(combined_data, "7500", chromatogram_data, "Qual", "STD")

combined_data <- combined_data |>
  mutate(quant_rt = case_when(
    instrument == "6500" ~ is_rt + quant_offset_6500,
    instrument == "7500" ~ is_rt + quant_offset_7500,
    TRUE ~ NA_real_
  )) |>
  filter(sample != "STD")

# Open PDF file
pdf(file.path("figures", "supplemental_figures", "RvT4_Chemical_Noise_Analysis.pdf"),
  width = 8.5, height = 11, onefile = TRUE
)

# Initialize lists to store peak times and powers across all samples
current_instrument <- NULL
all_peak_data <- list()
scale_list <- seq(1.5, 3, 0.05) # Starting at 1.5 to skip noise scales

cat("Cross-scale consistency filter: peaks must be detected at >= ",
    min_scales_detected, " out of ", length(scale_list), " scales\n", sep = "")

# Initialize list to store results for JSON export
json_results <- list()

# Iterate over all samples in combined_data
for (i in 1:nrow(combined_data)) {
  row <- combined_data[i, ]

  # Set initial instrument if needed
  if (is.null(current_instrument)) {
    current_instrument <- row$instrument
  }

  # Check for instrument change
  if (row$instrument != current_instrument && i > 1) {
    # Print summary plots for previous instrument
    print_summary_plots(current_instrument, all_peak_data,
      csv_prefix = "RvT4_peak_data_",
      compound_name = "Quant", precursor_mz = "361.1", product_mz = "211.1")

    # Reset for new instrument
    current_instrument <- row$instrument
    all_peak_data <- list()
  }

  # Process a sample (single row in combined_data)
  sample_name <- row$sample

  # Filter the current sample's chromatogram data
  current_traces <- chromatogram_data |>
    filter(Sample == sample_name, instrument == row$instrument) |>
    mutate(Compound = fct_relevel(Compound, "IS", "Quant", "Qual"))

  # Extract the Quant trace for wavelet analysis
  y <- current_traces |>
    filter(Compound == "Quant") |>
    pull(Intensity)

  t <- current_traces |>
    filter(Compound == "Quant") |>
    pull(RT)

  # Perform Continuous Wavelet Transform (CWT)
  wCoefs <- cwt(y, scales = scale_list, wavelet = "mexh")

  # Truncate negative coefficients
  wavelet_result_truncated <- wCoefs
  wavelet_result_truncated[wavelet_result_truncated < 0] <- 0

  # Compute the wavelet power spectrum
  wavelet_power <- wavelet_result_truncated^2
  df <- as_tibble(wavelet_power)
  df$t <- t

  # Find peaks using the previously defined logic
  all_peaks <- list()
  for (scale_name in names(df)[-ncol(df)]) { # Exclude the 't' column

    scale_data <- df |>
      select(t, all_of(scale_name))

    # Find all peaks above a minimal floor (cross-scale filter applied later).
    # Using 1 (not 0) avoids zero-weight issues in merge_nearby_peaks.
    peaks <- find_peaks(scale_data[[scale_name]], scale_data$t,
      scale = as.numeric(scale_name),
      min_prominence_absolute = 1
    )

    # Store peaks for this scale
    if (nrow(peaks) > 0) {
      all_peaks[[scale_name]] <- data.frame(
        Scale = as.numeric(scale_name),
        PeakTime = peaks$time,
        PeakPower = peaks$power
      )
    }
  }

  # Combine peaks from all scales for this sample
  all_peaks_df <- bind_rows(all_peaks)
  proximity_threshold <- 0.03

  # Cross-scale consistency filtering: group peaks by time proximity
  # across scales and require detection at multiple scales. This
  # discriminates real chromatographic features (coherent across scales)
  # from electronic noise artifacts (incoherent, 1-2 scales).
  empty_peaks <- data.frame(
    index = integer(0), time = numeric(0), power = numeric(0),
    prominence = numeric(0), scale = numeric(0)
  )

  if (nrow(all_peaks_df) > 0 && "PeakTime" %in% names(all_peaks_df)) {
    # Remove any rows with NA/NaN times (from zero-prominence merging artifacts)
    all_peaks_df <- all_peaks_df[!is.na(all_peaks_df$PeakTime), ]
  }

  if (nrow(all_peaks_df) > 0 && "PeakTime" %in% names(all_peaks_df)) {
    # Sort by time for sequential grouping
    all_peaks_df <- all_peaks_df[order(all_peaks_df$PeakTime), ]

    # Assign groups: consecutive peaks within tolerance belong to same feature
    cross_scale_tolerance <- 0.05
    group_ids <- integer(nrow(all_peaks_df))
    group_ids[1] <- 1
    gid <- 1
    if (nrow(all_peaks_df) > 1) {
      for (j in 2:nrow(all_peaks_df)) {
        if (all_peaks_df$PeakTime[j] - all_peaks_df$PeakTime[j - 1] > cross_scale_tolerance) {
          gid <- gid + 1
        }
        group_ids[j] <- gid
      }
    }
    all_peaks_df$group <- group_ids

    # Summarize each group: count unique scales, weighted mean time, max power
    group_summary <- all_peaks_df |>
      group_by(group) |>
      summarise(
        n_scales = n_distinct(Scale),
        time = weighted.mean(PeakTime, PeakPower),
        power = max(PeakPower),
        scale = Scale[which.max(PeakPower)],
        .groups = "drop"
      ) |>
      filter(n_scales >= min_scales_detected) |>
      filter(abs(time - row$quant_rt) > proximity_threshold)

    cat("  ", sample_name, ": ", nrow(all_peaks_df), " per-scale peaks -> ",
        nrow(group_summary), " cross-scale features (>= ", min_scales_detected,
        " scales)\n", sep = "")

    if (nrow(group_summary) > 0) {
      dominant_peaks_df <- data.frame(
        index = seq_len(nrow(group_summary)),
        time = group_summary$time,
        power = group_summary$power,
        prominence = group_summary$power,
        scale = group_summary$scale
      )
    } else {
      dominant_peaks_df <- empty_peaks
    }
  } else {
    cat("  ", sample_name, ": 0 per-scale peaks\n", sep = "")
    dominant_peaks_df <- empty_peaks
  }

  # Step 3: Store Peak Data for the Current Sample
  all_peak_data[[sample_name]] <- dominant_peaks_df

  # Plot the current sample traces, power, and selected chemical noise peaks

  df_long <- df |>
    pivot_longer(cols = -t, names_to = "row", values_to = "power") |>
    mutate(scale = as.numeric(row))

  quant_offset <- if (row$instrument == "6500") quant_offset_6500 else quant_offset_7500

  # Create the plot
  plot1 <- current_traces |>
    filter(Compound %in% c("IS", "Quant")) |>
    ggplot(aes(x = RT, y = Intensity, color = Compound)) +
    geom_line() +
    labs(x = "RT", y = "Intensity") +
    ggtitle(label = bquote("RvT4 Quant (m/z 361.1" %->% "211.1): " ~ .(sample_name) ~ " (" * .(row$instrument) * ")")) +
    theme(strip.text.x = element_text(size = 12, face = "bold"))

  # Add vlines only if values are valid (not NA and finite)
  if (!is.na(row$is_rt) && is.finite(row$is_rt)) {
    plot1 <- plot1 + geom_vline(xintercept = row$is_rt, color = "black")
  }
  if (!is.na(row$quant_rt) && is.finite(row$quant_rt)) {
    plot1 <- plot1 + geom_vline(xintercept = row$quant_rt, color = "black", linetype = "dashed")
  }
  if (nrow(dominant_peaks_df) > 0 && any(!is.na(dominant_peaks_df$time))) {
    valid_times <- dominant_peaks_df$time[!is.na(dominant_peaks_df$time) & is.finite(dominant_peaks_df$time)]
    if (length(valid_times) > 0) {
      plot1 <- plot1 + geom_vline(xintercept = valid_times, color = "red", linetype = "dashed")
    }
  }

  # Define custom scales for each facet
  scales_y <- list(
    "IS" = scale_y_continuous(labels = inten_label),
    "Quant" = scale_y_continuous(labels = inten_label)
  )

  # Apply facets with custom scales
  plot1 <- plot1 +
    ggh4x::facet_wrap2(
      ~ factor(ifelse(Compound == "IS", "IS", "Quant"),
        levels = c("IS", "Quant")
      ),
      nrow = 2,
      scales = "free_y",
      axes = "y"
    ) +
    ggh4x::facetted_pos_scales(y = scales_y) +
    scale_color_discrete(breaks = c("IS", "Quant"))

  plot2 <- df_long |>
    ggplot(aes(x = t, y = power, color = scale, group = row)) +
    geom_line() +
    labs(x = "Retention Time (min)", y = "Wavelet Power") +
    ggtitle("Quant") +
    scale_color_gradient(name = "Scale", low = "darkblue", high = "darkred") +
    scale_y_continuous(labels = inten_label) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12),
      legend.box.just = "right",
      legend.position = "right"
    )

  if (!is.na(row$quant_rt) && is.finite(row$quant_rt)) {
    plot2 <- plot2 + geom_vline(xintercept = row$quant_rt, color = "black", linetype = "dashed")
    # Add gray dashed lines for proximity threshold bounds
    plot2 <- plot2 + geom_vline(xintercept = row$quant_rt - proximity_threshold, color = "gray", linetype = "dashed")
    plot2 <- plot2 + geom_vline(xintercept = row$quant_rt + proximity_threshold, color = "gray", linetype = "dashed")
  }
  if (nrow(dominant_peaks_df) > 0 && any(!is.na(dominant_peaks_df$time))) {
    valid_times <- dominant_peaks_df$time[!is.na(dominant_peaks_df$time) & is.finite(dominant_peaks_df$time)]
    if (length(valid_times) > 0) {
      plot2 <- plot2 + geom_vline(xintercept = valid_times, color = "red", linetype = "dashed")
    }
  }

  # Calculate x-axis range and breaks for alignment
  x_range <- range(current_traces$RT, na.rm = TRUE)
  x_breaks <- seq(floor(x_range[1] * 10) / 10,
    ceiling(x_range[2] * 10) / 10,
    by = 0.1
  )

  # Apply matching x-axis settings to both plots
  plot1 <- plot1 +
    scale_x_continuous(breaks = x_breaks, limits = c(min(x_breaks), max(x_breaks)))

  plot2 <- plot2 +
    scale_x_continuous(breaks = x_breaks, limits = c(min(x_breaks), max(x_breaks)))

  # Print the combined plot for current sample
  combined_plot <- plot1 / plot2
  print(combined_plot)

  # Store data for JSON export
  # Align IS and Quant intensities onto the Quant time axis (t)
  rt_IS_raw <- current_traces |>
    filter(Compound == "IS") |>
    pull(RT)

  is_raw <- current_traces |>
    filter(Compound == "IS") |>
    pull(Intensity)

  rt_Quant_raw <- current_traces |>
    filter(Compound == "Quant") |>
    pull(RT)

  quant_raw <- current_traces |>
    filter(Compound == "Quant") |>
    pull(Intensity)

  intensity_IS <- approx(rt_IS_raw, is_raw, xout = t, rule = 2)$y
  intensity_Quant <- approx(rt_Quant_raw, quant_raw, xout = t, rule = 2)$y

  json_results[[i]] <- list(
    instrument      = row$instrument,
    sample          = row$sample,
    is_rt           = row$is_rt,
    quant_rt        = row$quant_rt,
    quant_offset    = quant_offset,
    # Single canonical time axis for wavelet and chromatogram
    wavelet_time    = t,
    intensity_IS    = intensity_IS,
    intensity_Quant = intensity_Quant,
    wavelet_power   = wavelet_power,
    scales          = scale_list,
    detected_peaks  = dominant_peaks_df |> select(time, power, scale)
  )
}

# After the loop ends, print final instrument's summary plots
print_summary_plots(current_instrument, all_peak_data,
  csv_prefix = "RvT4_peak_data_",
  compound_name = "Quant", precursor_mz = "361.1", product_mz = "211.1")

dev.off()

# Export results to JSON
json_text <- toJSON(json_results, pretty = TRUE, auto_unbox = TRUE)
write(json_text, file = file.path("results", "json", "RvT4_chemical_noise.json"))
