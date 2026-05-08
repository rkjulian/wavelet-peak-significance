# Ketamine Chemical Noise Supplemental Plots (plot-only)
#
# Reads pre-computed results from Ketamine_quant_chemical_noise.json and
# generates the multi-page supplemental PDF. Produces the same output as
# Ketamine-quant-chemical-noise.R without re-running the CWT analysis.

library(tidyverse)
library(viridis)
library(jsonlite)
library(patchwork)
library(ggh4x)
library(forcats)

# Clear the environment
rm(list = ls())

source(file.path("code", "mc_wavelet_plot_functions.R"))

# -----------------------------------------------------------------------------
# 1. Load pre-computed results from JSON
# -----------------------------------------------------------------------------

cat("Loading pre-computed results from Ketamine_quant_chemical_noise.json...\n")

json_results <- fromJSON(file.path("results", "json",
                                   "Ketamine_quant_chemical_noise.json"),
                         simplifyVector = FALSE)

cat("Loaded", length(json_results), "samples\n")

proximity_threshold <- 0.03

# Open multi-page PDF
pdf(file.path("figures", "supplemental_figures",
              "Ketamine_Quant_Chemical_Noise_Analysis.pdf"),
    width = 8.5, height = 11)

all_peak_data <- list()

# -----------------------------------------------------------------------------
# 2. Generate plots from saved data
# -----------------------------------------------------------------------------

for (i in seq_along(json_results)) {
    result <- json_results[[i]]

    # Extract data from JSON
    instrument_name <- result$instrument
    sample_name     <- result$sample
    t               <- unlist(result$wavelet_time)
    intensity_IS    <- unlist(result$intensity_IS)
    intensity_Quant <- unlist(result$intensity_Quant)
    is_rt           <- result$is_rt
    quant_rt        <- result$quant_rt
    wavelet_power   <- do.call(rbind, lapply(result$wavelet_power, unlist))
    scales          <- unlist(result$scales)
    detected_peaks  <- result$detected_peaks

    # Extract detected peak data
    if (is.null(detected_peaks) || length(detected_peaks) == 0) {
        peak_times <- numeric(0)
        peak_df <- tibble(time = numeric(0), power = numeric(0),
                          scale = numeric(0), prominence = numeric(0))
    } else if (is.data.frame(detected_peaks)) {
        peak_times <- detected_peaks$time
        peak_df <- as_tibble(detected_peaks)
        if (!"prominence" %in% names(peak_df)) {
            peak_df$prominence <- NA_real_
        }
    } else {
        peak_times <- unlist(lapply(detected_peaks, `[[`, "time"),
                            use.names = FALSE)
        peak_df <- bind_rows(lapply(detected_peaks, function(p) {
            tibble(
                time = p$time %||% NA_real_,
                power = p$power %||% NA_real_,
                scale = p$scale %||% NA_real_,
                prominence = p$prominence %||% NA_real_
            )
        }))
    }

    # Accumulate peak data for summary plots
    all_peak_data[[sample_name]] <- peak_df

    cat("Plotting sample:", sample_name, "(", instrument_name, ")\n")

    # Align lengths of time and intensity vectors
    n <- min(length(t), length(intensity_IS), length(intensity_Quant))
    t_use <- t[seq_len(n)]
    intensity_IS_use <- intensity_IS[seq_len(n)]
    intensity_Quant_use <- intensity_Quant[seq_len(n)]

    # Create chromatogram data frame with panel factor for faceting
    chromatogram_data <- tibble(
        RT = rep(t_use, 2),
        Intensity = c(intensity_IS_use, intensity_Quant_use),
        Compound = factor(rep(c("IS", "Quant"), each = length(t_use)),
                         levels = c("IS", "Quant")),
        panel = factor(rep(c("IS", "Quant"), each = length(t_use)),
                      levels = c("IS", "Quant"))
    )

    # ------------------------------------------------------------------
    # Plot 1: Chromatogram traces with detected peaks
    # ------------------------------------------------------------------

    plot1 <- chromatogram_data %>%
        ggplot(aes(x = RT, y = Intensity, color = Compound)) +
        geom_line() +
        labs(x = "RT", y = "Intensity") +
        ggtitle(label = bquote(
            "Ketamine Quant (m/z 238" %->% "125): " ~
            .(sample_name) ~ " (" * .(instrument_name) * ")")) +
        theme(strip.text.x = element_text(size = 12, face = "bold"))

    # IS RT line appears in both facets
    if (!is.na(is_rt) && is.finite(is_rt)) {
        plot1 <- plot1 + geom_vline(xintercept = is_rt, color = "black")
    }

    # Quant RT and chemical noise lines only in the Quant facet
    quant_panel <- factor("Quant", levels = c("IS", "Quant"))
    if (!is.na(quant_rt) && is.finite(quant_rt)) {
        vline_df <- data.frame(xint = quant_rt, panel = quant_panel)
        plot1 <- plot1 + geom_vline(data = vline_df, aes(xintercept = xint),
                                     color = "black", linetype = "dashed")
    }
    if (length(peak_times) > 0) {
        valid_times <- peak_times[!is.na(peak_times) & is.finite(peak_times)]
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

    # ------------------------------------------------------------------
    # Plot 2: Wavelet power spectrum
    # ------------------------------------------------------------------

    df <- as_tibble(wavelet_power, .name_repair = "unique")
    colnames(df) <- as.character(scales)
    df$t <- t

    df_long <- df %>%
        pivot_longer(cols = -t, names_to = "row", values_to = "power") %>%
        mutate(scale = as.numeric(row))

    plot2 <- df_long %>%
        ggplot(aes(x = t, y = power, color = scale, group = row)) +
        geom_line() +
        labs(x = "Retention Time (min)", y = "Wavelet Power") +
        ggtitle("Quant") +
        scale_color_gradient(name = "Scale",
                             low = "darkblue", high = "darkred") +
        scale_y_continuous(labels = inten_label) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 16),
              plot.subtitle = element_text(hjust = 0.5),
              axis.text = element_text(size = 11),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 11),
              legend.title = element_text(size = 12),
              legend.box.just = "right",
              legend.position = "right")

    if (!is.na(quant_rt) && is.finite(quant_rt)) {
        plot2 <- plot2 +
            geom_vline(xintercept = quant_rt, color = "black",
                       linetype = "dashed") +
            geom_vline(xintercept = quant_rt - proximity_threshold,
                       color = "gray", linetype = "dashed") +
            geom_vline(xintercept = quant_rt + proximity_threshold,
                       color = "gray", linetype = "dashed")
    }
    if (length(peak_times) > 0) {
        valid_times <- peak_times[!is.na(peak_times) & is.finite(peak_times)]
        if (length(valid_times) > 0) {
            plot2 <- plot2 +
                geom_vline(xintercept = valid_times, color = "red",
                           linetype = "dashed")
        }
    }

    # Calculate x-axis range and breaks for alignment
    x_range <- range(t, na.rm = TRUE)
    x_span <- diff(x_range)
    x_step <- if (x_span > 10) 2 else if (x_span > 5) 1 else 0.5
    x_breaks <- seq(floor(x_range[1] / x_step) * x_step,
                    ceiling(x_range[2] / x_step) * x_step,
                    by = x_step)
    x_limits <- c(min(x_breaks), max(x_breaks))

    # Apply matching x-axis settings to both plots
    plot1 <- plot1 +
        scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
        coord_cartesian(xlim = x_limits)
    plot2 <- plot2 +
        scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
        coord_cartesian(xlim = x_limits)

    # Print the combined plot for current sample
    combined_plot <- plot1 / plot2
    print(combined_plot)
}

# Print summary plots
print_summary_plots("sciex", all_peak_data,
    csv_prefix = "Ketamine_quant_peak_data_",
    compound_name = "Quant",
    precursor_mz = "238", product_mz = "125",
    analyte_name = "Ketamine")

dev.off()

cat("\nSupplemental chemical noise plots saved to:\n")
cat("  figures/supplemental_figures/Ketamine_Quant_Chemical_Noise_Analysis.pdf\n")
