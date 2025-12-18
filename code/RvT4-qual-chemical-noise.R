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

# Set seed for reproducibility
set.seed(42)

# NOTE: This script analyzes the Qual (qualifier) transition for chemical noise
# while using the Quant transition as the retention time reference standard

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

# Helper function to create a data frame for instrument data
create_instrument_data <- function(instrument, file_list, sample_list, quant_index, qual_index, is_index, is_rt) {
  tibble(
    instrument = instrument,
    filename = file_list,
    sample = sample_list,
    quant_index = as.integer(quant_index),
    qual_index = as.integer(qual_index),
    is_index = as.integer(is_index),
    is_rt = as.numeric(is_rt)
  )
}

# Helper function to create an empty chromatogram data frame
empty_chrom <- function() {
  tibble(
    RT = numeric(),
    Intensity = numeric(),
    Compound = character(),
    Sample = character(),
    instrument = character()
  )
}

# Helper function to process SRM data for a single row
process_srm_data <- function(row, data_dir) {
  srm_filename <- file.path(data_dir, "mzML", row$filename)
  srm <- readSRMData(srm_filename)

  compounds <- c("Quant", "Qual", "IS")
  indices <- c(row$quant_index, row$qual_index, row$is_index)

  chromatogram_data <- map2_dfr(compounds, indices, ~ {
    tibble(
      RT = rtime(srm[.y]),
      Intensity = intensity(srm[.y]),
      Compound = .x,
      Sample = row$sample
    )
  })

  is_rt <- chromatogram_data |>
    filter(Compound == "IS") |>
    slice(which.max(Intensity)) |>
    pull(RT)

  list(chromatogram_data = chromatogram_data, is_rt = is_rt)
}

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

# Calculate quant offsets (using Quant as RT reference standard)
calculate_quant_offset <- function(data, instrument) {
  quant_rt <- chromatogram_data |>
    filter(instrument == !!instrument, Sample == "STD", Compound == "Quant") |>
    slice(which.max(Intensity)) |>
    pull(RT)

  is_rt <- data |>
    filter(instrument == !!instrument, sample == "STD") |>
    pull(is_rt)

  quant_rt - is_rt
}

quant_offset_6500 <- calculate_quant_offset(combined_data, "6500")
quant_offset_7500 <- calculate_quant_offset(combined_data, "7500")

combined_data <- combined_data |>
  mutate(
    quant_rt = case_when(
      instrument == "6500" ~ is_rt + quant_offset_6500,
      instrument == "7500" ~ is_rt + quant_offset_7500,
      TRUE ~ NA_real_
    ),
    qual_expected_rt = case_when(
      instrument == "7500" ~ 15.37,
      TRUE ~ NA_real_
    )
  ) |>
  filter(sample != "STD")

# Function to print summary plots
print_summary_plots <- function(instrument, peak_data) {
  # Combine all samples' data into a single data frame
  peak_data_combined <- bind_rows(peak_data, .id = "Sample")

  # Save complete data to CSV
  peak_data_combined |>
    select(-prominence) |>
    write_csv(file.path("results", "tables", paste0("RvT4_qual_peak_data_", instrument, ".csv")))

  # Filter out non-finite values for plotting
  peak_data_valid <- peak_data_combined |>
    filter(is.finite(time), is.finite(power))

  # Create time distribution plot
  plot3 <- ggplot(peak_data_valid, aes(x = time)) +
    geom_histogram(bins = 25, fill = "blue", color = "black", alpha = 0.7) +
    labs(
      title = bquote("Time Distribution - RvT4 Qual (m/z 361.1" %->% "193.1) - " ~ .(instrument)),
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
      title = bquote("Peak Power Distribution - RvT4 Qual (m/z 361.1" %->% "193.1) - " ~ .(instrument)),
      x = "Power (dB)", y = "Density"
    ) +
    theme_minimal()

  # Arrange plots vertically
  dist_plots <- plot3 / plot4
  print(dist_plots)
}

find_peaks <- function(data, t, scale, min_power = 1e-5, min_prominence_absolute = 1e4) {
  window_size <- 3
  smoothed_data <- stats::filter(data, rep(1 / window_size, window_size), sides = 2)
  smoothed_data[is.na(smoothed_data)] <- 0

  peaks <- numeric()
  prominences <- numeric()

  for (i in (window_size + 1):(length(smoothed_data) - window_size)) {
    if (smoothed_data[i] >= 0.95 * max(smoothed_data[(i - window_size):(i + window_size)]) &&
      smoothed_data[i] > min_power) {
      # Calculate prominence
      left_bound <- i
      while (left_bound > 1) {
        if (smoothed_data[left_bound - 1] > smoothed_data[i]) break
        left_bound <- left_bound - 1
      }

      right_bound <- i
      while (right_bound < length(smoothed_data)) {
        if (smoothed_data[right_bound + 1] > smoothed_data[i]) break
        right_bound <- right_bound + 1
      }

      left_min <- min(smoothed_data[left_bound:i])
      right_min <- min(smoothed_data[i:right_bound])
      reference_height <- max(left_min, right_min)
      prominence <- smoothed_data[i] - reference_height

      if (prominence >= min_prominence_absolute) {
        peaks <- c(peaks, i)
        prominences <- c(prominences, prominence)
      }
    }
  }

  peak_data <- data.frame(
    index = peaks,
    time = t[peaks],
    power = data[peaks],
    prominence = prominences,
    scale = scale
  )

  # Add within-scale merging
  if (nrow(peak_data) > 0) {
    peak_data <- merge_nearby_peaks(peak_data, min_separation = 0.04) # Tighter threshold for same scale
  }

  return(peak_data)
}

merge_nearby_peaks <- function(peak_data, min_separation = 0.05) {
  # Check if input is valid
  if (!is.data.frame(peak_data) || nrow(peak_data) == 0) {
    return(peak_data)
  }

  if (nrow(peak_data) <= 1) {
    return(peak_data)
  }

  # Check if required columns exist
  required_cols <- c("index", "time", "power", "prominence", "scale")
  if (!all(required_cols %in% names(peak_data))) {
    stop(
      "Missing required columns in peak_data. Need: ",
      paste(required_cols, collapse = ", ")
    )
  }

  # Sort by prominence if it exists
  if (length(peak_data$prominence) > 0) {
    peak_data <- peak_data[order(-peak_data$prominence), ]
  }

  merged <- list()
  used_indices <- numeric()

  for (i in seq_len(nrow(peak_data))) {
    if (length(peak_data$index[i]) == 0 || is.na(peak_data$index[i])) next
    if (peak_data$index[i] %in% used_indices) next

    nearby <- which(abs(peak_data$time - peak_data$time[i]) < min_separation)
    nearby <- nearby[!(peak_data$index[nearby] %in% used_indices)]

    if (length(nearby) > 0) {
      weights <- peak_data$prominence[nearby]
      merged[[length(merged) + 1]] <- data.frame(
        index = peak_data$index[i],
        time = weighted.mean(peak_data$time[nearby], weights),
        power = max(peak_data$power[nearby]),
        prominence = max(peak_data$prominence[nearby]),
        scale = peak_data$scale[which.max(peak_data$prominence[nearby])]
      )
      used_indices <- c(used_indices, peak_data$index[nearby])
    }
  }

  # Check if any peaks were merged
  if (length(merged) == 0) {
    return(peak_data[0, ]) # Return empty data frame with same structure
  }

  result <- dplyr::bind_rows(merged)
  return(result)
}

# Open PDF file
pdf(file.path("figures", "supplemental_figures", "RvT4_Qual_Chemical_Noise_Analysis.pdf"),
  width = 8.5, height = 11, onefile = TRUE
)

# Initialize lists to store peak times and powers across all samples
current_instrument <- NULL
all_peak_data <- list()
scale_list <- seq(1.5, 3, 0.05) # Starting at 1.5 to skip noise scales

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
    print_summary_plots(current_instrument, all_peak_data)

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

  # Extract the Qual trace for wavelet analysis
  y <- current_traces |>
    filter(Compound == "Qual") |>
    pull(Intensity)

  t <- current_traces |>
    filter(Compound == "Qual") |>
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

    # Find peaks
    peaks <- find_peaks(scale_data[[scale_name]], scale_data$t,
      scale = as.numeric(scale_name)
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

  # Step 1: Remove Peaks Near `quant_rt` (using Quant as RT reference)
  proximity_threshold <- 0.03
  filtered_peaks <- all_peaks_df |>
    filter(abs(PeakTime - row$quant_rt) > proximity_threshold)

  # Step 1b: For 7500 instrument, also remove peaks near the expected Qual compound at 15.38 min
  if (row$instrument == "7500" && !is.na(row$qual_expected_rt)) {
    filtered_peaks <- filtered_peaks |>
      filter(abs(PeakTime - row$qual_expected_rt) > proximity_threshold)
  }

  # Step 2: Merge nearby peaks across all scales
  if (nrow(filtered_peaks) > 0) {
    # Transform filtered_peaks to match required column names
    peaks_for_merging <- data.frame(
      index = 1:nrow(filtered_peaks),
      time = filtered_peaks$PeakTime,
      power = filtered_peaks$PeakPower,
      prominence = filtered_peaks$PeakPower, # Using power as prominence
      scale = filtered_peaks$Scale
    )
    dominant_peaks_df <- merge_nearby_peaks(peaks_for_merging)
  } else {
    dominant_peaks_df <- data.frame(
      index = integer(0),
      time = numeric(0),
      power = numeric(0),
      prominence = numeric(0),
      scale = numeric(0)
    )
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
    filter(Compound %in% c("IS", "Qual")) |>
    ggplot(aes(x = RT, y = Intensity, color = Compound)) +
    geom_line() +
    labs(x = "RT", y = "Intensity") +
    ggtitle(label = bquote("RvT4 Qual (m/z 361.1" %->% "193.1): " ~ .(sample_name) ~ " (" * .(row$instrument) * ")")) +
    theme(strip.text.x = element_text(size = 12, face = "bold"))

  # Add vlines only if values are valid (not NA and finite)
  if (!is.na(row$is_rt) && is.finite(row$is_rt)) {
    plot1 <- plot1 + geom_vline(xintercept = row$is_rt, color = "black")
  }
  if (!is.na(row$quant_rt) && is.finite(row$quant_rt)) {
    plot1 <- plot1 + geom_vline(xintercept = row$quant_rt, color = "black", linetype = "dashed")
  }
  if (row$instrument == "7500" && !is.na(row$qual_expected_rt) && is.finite(row$qual_expected_rt)) {
    plot1 <- plot1 + geom_vline(xintercept = row$qual_expected_rt, color = "blue", linetype = "dashed")
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
    "Qual" = scale_y_continuous(labels = inten_label)
  )

  # Apply facets with custom scales
  plot1 <- plot1 +
    ggh4x::facet_wrap2(
      ~ factor(ifelse(Compound == "IS", "IS", "Qual"),
        levels = c("IS", "Qual")
      ),
      nrow = 2,
      scales = "free_y",
      axes = "y"
    ) +
    ggh4x::facetted_pos_scales(y = scales_y) +
    scale_color_discrete(breaks = c("IS", "Qual"))

  plot2 <- df_long |>
    ggplot(aes(x = t, y = power, color = scale, group = row)) +
    geom_line() +
    labs(x = "Retention Time (min)", y = "Wavelet Power") +
    ggtitle("Qual") +
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
  if (row$instrument == "7500" && !is.na(row$qual_expected_rt) && is.finite(row$qual_expected_rt)) {
    plot2 <- plot2 + geom_vline(xintercept = row$qual_expected_rt, color = "blue", linetype = "dashed")
    # Add light blue dashed lines for proximity threshold bounds
    plot2 <- plot2 + geom_vline(xintercept = row$qual_expected_rt - proximity_threshold, color = "lightblue", linetype = "dashed")
    plot2 <- plot2 + geom_vline(xintercept = row$qual_expected_rt + proximity_threshold, color = "lightblue", linetype = "dashed")
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
  # Align IS, Quant, and Qual intensities onto the Qual time axis (t)
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

  rt_Qual_raw <- current_traces |>
    filter(Compound == "Qual") |>
    pull(RT)

  qual_raw <- current_traces |>
    filter(Compound == "Qual") |>
    pull(Intensity)

  intensity_IS <- approx(rt_IS_raw, is_raw, xout = t, rule = 2)$y
  intensity_Quant <- approx(rt_Quant_raw, quant_raw, xout = t, rule = 2)$y
  intensity_Qual <- approx(rt_Qual_raw, qual_raw, xout = t, rule = 2)$y

  json_results[[i]] <- list(
    instrument       = row$instrument,
    sample           = row$sample,
    is_rt            = row$is_rt,
    quant_rt         = row$quant_rt,
    qual_expected_rt = if (row$instrument == "7500") row$qual_expected_rt else NULL,
    quant_offset     = quant_offset,
    # Single canonical time axis for wavelet and chromatogram (Qual RT)
    wavelet_time     = t,
    intensity_IS     = intensity_IS,
    intensity_Quant  = intensity_Quant,
    intensity_Qual   = intensity_Qual,
    wavelet_power    = wavelet_power,
    scales           = scale_list,
    detected_peaks   = dominant_peaks_df |> select(time, power, scale)
  )
}

# After the loop ends, print final instrument's summary plots
print_summary_plots(current_instrument, all_peak_data)

dev.off()

# Export results to JSON
json_text <- toJSON(json_results, pretty = TRUE, auto_unbox = TRUE)
write(json_text, file = file.path("results", "json", "RvT4_qual_chemical_noise.json"))
