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
    "waters" = list(
        file_list = c(
            "20Feb24 Blank Plasma 1.mzML",
            "20Feb24 Blank Plasma 2.mzML",
            "20Feb24 Blank Plasma 3.mzML",
            "20Feb24 Spiked Plasma 0-25X 1.mzML",
            "20Feb24 Spiked Plasma 0-25X 2.mzML",
            "20Feb24 Spiked Plasma 0-25X 3.mzML",
            "20Feb24 Spiked Plasma 0-5X 1.mzML",
            "20Feb24 Spiked Plasma 0-5X 2.mzML",
            "20Feb24 Spiked Plasma 0-5X 3.mzML",
            "20Feb24 Spiked Plasma 1X 1.mzML",
            "20Feb24 Spiked Plasma 1X 2.mzML",
            "20Feb24 Spiked Plasma 1X 3.mzML",
            "20Feb24 Spiked Plasma 2X 1.mzML",
            "20Feb24 Spiked Plasma 2X 2.mzML",
            "20Feb24 Spiked Plasma 2X 3.mzML",
            "20Feb24 Spiked Plasma 4X 1.mzML",
            "20Feb24 Spiked Plasma 4X 2.mzML",
            "20Feb24 Spiked Plasma 4X 3.mzML",
            "20Feb24 Spiked Plasma 8X 1.mzML",
            "20Feb24 Spiked Plasma 8X 2.mzML",
            "20Feb24 Spiked Plasma 8X 3.mzML",
            "20Feb24 Spiked Plasma 16X 1.mzML",
            "20Feb24 Spiked Plasma 16X 2.mzML",
            "20Feb24 Spiked Plasma 16X 3.mzML",
            "20Feb24 Spiked Plasma 32X 1.mzML",
            "20Feb24 Spiked Plasma 32X 2.mzML",
            "20Feb24 Spiked Plasma 32X 3.mzML",
            "20Feb24 Spiked Plasma 64X 1.mzML",
            "20Feb24 Spiked Plasma 64X 2.mzML",
            "20Feb24 Spiked Plasma 64X 3.mzML"
        ),
        sample_list = c(
            "Blank 1", "Blank 2", "Blank 3",
            "0_25X 1", "0_25X 2", "0_25X 3",
            "0_5X 1", "0_5X 2", "0_5X 3",
            "1X 1", "1X 2", "1X 3",
            "2X 1", "2X 2", "2X 3",
            "4X 1", "4X 2", "4X 3",
            "8X 1", "8X 2", "8X 3",
            "16X 1", "16X 2", "16X 3",
            "32X 1", "32X 2", "32X 3",
            "64X 1", "64X 2", "64X 3"
        ),
        quant_index = 8,
        qual_index = 7,
        is_index = 9,
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

# Process data
data_dir <- file.path("data", "Waters_Tox")
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
quant_offset_waters <- calculate_quant_offset(combined_data, "waters", chromatogram_data, "Quant", "64X 1")

combined_data <- combined_data |>
    mutate(quant_rt = is_rt + quant_offset_waters)

# Open PDF file
pdf(file.path("figures", "supplemental_figures", "Clonidine_Quant_Chemical_Noise_Analysis.pdf"),
    width = 8.5, height = 11, onefile = TRUE
)

# Initialize lists to store peak times and powers across all samples
all_peak_data <- list()
scale_list <- seq(1.5, 3, 0.05) # Matches narrow chemical noise and real peak widths

cat("Cross-scale consistency filter: peaks must be detected at >= ",
    min_scales_detected, " out of ", length(scale_list), " scales\n", sep = "")

# Initialize list to store results for JSON export
json_results <- list()

# Iterate over all samples in combined_data
for (i in 1:nrow(combined_data)) {
    row <- combined_data[i, ]

    # Process a sample (single row in combined_data)
    sample_name <- row$sample

    # Filter the current sample's chromatogram data
    current_traces <- chromatogram_data |>
        filter(Sample == sample_name, instrument == row$instrument, RT >= 0.6, RT <= 3.0) |>
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
        # Time-binned cross-scale counting: snap each detected peak to its
        # nearest data time point, then count unique scales per bin. This
        # evaluates each location independently, avoiding the chaining
        # problem where dense noise peaks bridge gaps between real features
        # in sequential grouping.
        all_peaks_df$bin_idx <- sapply(all_peaks_df$PeakTime,
                                       function(pt) which.min(abs(t - pt)))

        bin_summary <- all_peaks_df |>
            group_by(bin_idx) |>
            summarise(
                n_scales = n_distinct(Scale),
                max_power = max(PeakPower),
                best_scale = Scale[which.max(PeakPower)],
                .groups = "drop"
            )

        # Build per-time-point vectors
        scale_counts <- integer(length(t))
        max_power_at_t <- numeric(length(t))
        best_scale_at_t <- numeric(length(t))
        scale_counts[bin_summary$bin_idx] <- bin_summary$n_scales
        max_power_at_t[bin_summary$bin_idx] <- bin_summary$max_power
        best_scale_at_t[bin_summary$bin_idx] <- bin_summary$best_scale

        # Filter out electronic noise spikes BEFORE local max detection.
        # Spikes (narrow in raw trace) produce coherent CWT power across
        # all scales but are not chromatographic features. Zeroing their
        # power first prevents them from suppressing nearby real peaks
        # in the local max competition.
        min_peak_width <- 3  # minimum half-height width in data points
        for (b in seq_len(nrow(bin_summary))) {
            k <- bin_summary$bin_idx[b]
            peak_val <- y[k]
            half_height <- peak_val / 2

            left <- k
            while (left > 1 && y[left - 1] > half_height) left <- left - 1

            right <- k
            while (right < length(y) && y[right + 1] > half_height) right <- right + 1

            if ((right - left + 1) < min_peak_width) {
                max_power_at_t[k] <- 0
                scale_counts[k] <- 0L
            }
        }

        # Zero power at all non-significant bins so they can't block
        # significant neighbors in the local max competition. A bin with
        # high wavelet power but low scale count (e.g., 2 scales) is not
        # a cross-scale-consistent feature and should not suppress one.
        significant <- which(scale_counts >= min_scales_detected)
        non_significant <- which(scale_counts > 0L & scale_counts < min_scales_detected)
        max_power_at_t[non_significant] <- 0

        cat("  ", sample_name, ": ", nrow(all_peaks_df), " per-scale peaks -> ",
            length(significant), " significant time points", sep = "")

        if (length(significant) > 0) {
            # Identify distinct features as local maxima in max_power_at_t
            # among significant time points. This avoids daisy-chain merging
            # where bridge points between real features (with >= 5 scales but
            # lower power) link distant peaks into one mega-group.
            half_window <- max(round(0.02 / mean(diff(t))), 3)

            feature_idx <- integer(0)
            for (k in significant) {
                nbr_start <- max(1, k - half_window)
                nbr_end <- min(length(t), k + half_window)
                if (max_power_at_t[k] >= max(max_power_at_t[nbr_start:nbr_end])) {
                    feature_idx <- c(feature_idx, k)
                }
            }

            group_summary <- tibble(
                time = t[feature_idx],
                power = max_power_at_t[feature_idx],
                n_scales = as.integer(scale_counts[feature_idx]),
                scale = best_scale_at_t[feature_idx]
            ) |>
                filter(abs(time - row$quant_rt) > proximity_threshold)

            cat(" -> ", nrow(group_summary), " cross-scale features (>= ",
                min_scales_detected, " scales)\n", sep = "")

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
            cat(" -> 0 cross-scale features\n", sep = "")
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

    quant_offset <- quant_offset_waters

    # Create facet column for the chromatogram data
    plot_traces <- current_traces |>
        filter(Compound %in% c("IS", "Quant")) |>
        mutate(panel = factor(ifelse(Compound == "IS", "IS", "Quant"),
                              levels = c("IS", "Quant")))

    # Create the plot
    plot1 <- plot_traces |>
        ggplot(aes(x = RT, y = Intensity, color = Compound)) +
        geom_line() +
        labs(x = "RT", y = "Intensity") +
        ggtitle(label = bquote("Clonidine Quant (m/z 232.2" %->% "215.1): " ~ .(sample_name) ~ " (" * .(row$instrument) * ")")) +
        theme(strip.text.x = element_text(size = 12, face = "bold"))

    # IS RT line appears in both facets
    if (!is.na(row$is_rt) && is.finite(row$is_rt)) {
        plot1 <- plot1 + geom_vline(xintercept = row$is_rt, color = "black")
    }
    # Quant RT and chemical noise lines only in the Quant facet
    quant_panel <- factor("Quant", levels = c("IS", "Quant"))
    if (!is.na(row$quant_rt) && is.finite(row$quant_rt)) {
        vline_df <- data.frame(xint = row$quant_rt, panel = quant_panel)
        plot1 <- plot1 + geom_vline(data = vline_df, aes(xintercept = xint),
                                     color = "black", linetype = "dashed")
    }
    if (nrow(dominant_peaks_df) > 0 && any(!is.na(dominant_peaks_df$time))) {
        valid_times <- dominant_peaks_df$time[!is.na(dominant_peaks_df$time) & is.finite(dominant_peaks_df$time)]
        if (length(valid_times) > 0) {
            vline_df <- data.frame(xint = valid_times, panel = quant_panel)
            plot1 <- plot1 + geom_vline(data = vline_df, aes(xintercept = xint),
                                         color = "red", linetype = "dashed")
        }
    }

    # Define custom scales for each facet
    scales_y <- list(
        "IS" = scale_y_continuous(labels = inten_label),
        "Quant" = scale_y_continuous(labels = inten_label)
    )

    # Apply facets with custom scales
    plot1 <- plot1 +
        ggh4x::facet_wrap2(~ panel, nrow = 2, scales = "free_y", axes = "y") +
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
    x_span <- diff(x_range)
    x_step <- if (x_span > 10) 2 else if (x_span > 5) 1 else 0.5
    x_breaks <- seq(floor(x_range[1] / x_step) * x_step,
                    ceiling(x_range[2] / x_step) * x_step,
                    by = x_step
    )
    x_limits <- c(min(x_breaks), max(x_breaks))

    # Apply matching x-axis settings to both plots
    plot1 <- plot1 +
        scale_x_continuous(breaks = x_breaks, limits = x_limits, expand = c(0, 0))

    plot2 <- plot2 +
        scale_x_continuous(breaks = x_breaks, limits = x_limits, expand = c(0, 0))

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

# Print summary plots
print_summary_plots("waters", all_peak_data,
                    csv_prefix = "Clonidine_quant_peak_data_",
                    compound_name = "Quant", precursor_mz = "232.2", product_mz = "215.1")

dev.off()

# Export results to JSON
json_text <- toJSON(json_results, pretty = TRUE, auto_unbox = TRUE)
write(json_text, file = file.path("results", "json", "Clonidine_quant_chemical_noise.json"))
