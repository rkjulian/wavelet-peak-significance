library(MSnbase)
library(tidyverse)
library(MassSpecWavelet)
library(foreach)
library(doParallel)
library(truncnorm)
library(jsonlite)

# Clear the environment
rm(list = ls())

# Source shared functions
source(file.path("code", "mc_wavelet_functions.R"))

# Set seed for reproducibility
set.seed(42)

# -----------------------------------------------------------------------------
# 1. Read raw and meta data
# -----------------------------------------------------------------------------

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

# Calculate quant offset using highest concentration standard
quant_offset_waters <- calculate_quant_offset(
    combined_data, "waters", chromatogram_data, "Quant", "64X 1"
)

combined_data <- combined_data |>
    mutate(quant_rt = is_rt + quant_offset_waters)

# -----------------------------------------------------------------------------
# 2. Chemical noise parameters (constant for single instrument)
# -----------------------------------------------------------------------------

# Electronic noise from Waters reagent blank (ROS estimates from
# waters-noise-pooled.csv, Clonidine Quant row):
blank_noise_mean <- 1779.4147437
blank_noise_std <- 1304.5668835

# Load chemical noise statistics (from Clonidine-quant-chemical-noise.R)
peak_data <- read_csv(
    file.path("results", "tables", "Clonidine_quant_peak_data_waters.csv")
)

# Derive components needed for constructing null model
power_density <- density(peak_data$power, from = 0) # KDE, truncated at zero
power_density$y <- power_density$y / sum(power_density$y) # Renormalize

# Compute lambda for Poisson event count from per-sample peak counts
events_per_sample <- peak_data |>
    group_by(Sample) |>
    summarise(n_events = n(), .groups = "drop")
lambda_events <- mean(events_per_sample$n_events)

# Number of additional random chemical noise events per null simulation.
# Justified by: if a sample contained no analyte and a contaminant was
# detected at the expected RT, at least one new event must have occurred.
num_new_events <- 1

cat("waters: lambda_events =", round(lambda_events, 2),
    "(from", nrow(events_per_sample), "samples)\n")

# -----------------------------------------------------------------------------
# 3. Wavelet analysis and Monte Carlo simulation
# -----------------------------------------------------------------------------

# Initialize an empty list to store results for each file
simulation_results <- list()

# Loop through combined data
for (i in 1:nrow(combined_data)) {
    row <- combined_data[i, ]

    instrument_name <- row$instrument
    sample_name <- row$sample
    expected_rt <- row$quant_rt

    cat("Sample: ", sample_name, "\n")

    # Get the time axis (trimmed to match chemical noise characterization)
    retention_time <- chromatogram_data[
        chromatogram_data$instrument == instrument_name &
        chromatogram_data$Sample == sample_name, ] |>
        filter(Compound == "Quant", RT >= 0.6, RT <= 3.0) |>
        pull(RT)

    # Get the sample intensity values
    sample_y <- chromatogram_data[
        chromatogram_data$instrument == instrument_name &
        chromatogram_data$Sample == sample_name, ] |>
        filter(Compound == "Quant", RT >= 0.6, RT <= 3.0) |>
        pull(Intensity)

    # Save time series length
    n_points <- length(retention_time)

    # Parameters for CWT
    # MassSpecWavelet uses the mexh() function which is evaluated over [-8,8]
    # the maximum scale is set by the Nyquist criteria, or the limit of
    # this implementation (maximum scale is 64). The Nyquist criteria requires
    # two points per scale value so the minimum width is 16 and scale is
    # (n_points-1)/16
    scales <- seq(1, min((n_points - 1) / 16, 64), length.out = 100)
    wavelet_type <- "mexh" # Mexican hat wavelet

    # Compute CWT for the sample signal
    wavelet_result <- cwt(sample_y, scales = scales, wavelet = wavelet_type)

    # Truncate negative coefficients
    wavelet_result[wavelet_result < 0] <- 0

    # Compute the observed power in sample
    adjusted_power <- abs(wavelet_result)^2

    # Setup for parallel processing using fork via mclapply backend
    # registerDoParallel(cores=) uses fork-based mclapply internally,
    # avoiding PSOCK serialization and explicit cluster management
    n_cores <- detectCores() - 4
    registerDoParallel(cores = n_cores)

    # Define parameters
    n_sim <- 1000000 # Total simulations
    sims_per_worker <- ceiling(n_sim / n_cores) # Each worker accumulates internally

    # Initialize the count matrix
    counts <- matrix(0, nrow = nrow(adjusted_power), ncol = ncol(adjusted_power))

    cat("Running", n_sim, "simulations across", n_cores,
        "workers (", sims_per_worker, "per worker)\n")

    time_taken <- system.time({
        # Each worker runs sims_per_worker simulations and accumulates counts
        # internally, returning only the summed count matrix. This minimizes
        # inter-process communication from n_sim transfers to n_cores transfers.
        counts <- foreach(
            worker = seq_len(n_cores), .combine = "+"
        ) %dopar% {
            # Determine how many simulations this worker runs
            worker_n <- if (worker < n_cores) {
                sims_per_worker
            } else {
                n_sim - sims_per_worker * (n_cores - 1)
            }

            # Local accumulator for this worker
            local_counts <- matrix(0,
                nrow = nrow(adjusted_power),
                ncol = ncol(adjusted_power)
            )

            for (sim in seq_len(worker_n)) {
                # Draw number of existing chemical noise events from Poisson
                # Cap at pool size to prevent sample() error
                num_existing <- min(rpois(1, lambda_events), nrow(peak_data))

                if (num_existing > 0) {
                    sampled_times <- sample(peak_data$time, num_existing)
                    sampled_scales <- sample(peak_data$scale, num_existing)
                    sampled_powers <- sample_from_kde(power_density, num_existing)
                } else {
                    sampled_times <- numeric(0)
                    sampled_scales <- numeric(0)
                    sampled_powers <- numeric(0)
                }

                # Add new Gaussian event(s) at random positions
                new_times <- runif(
                    num_new_events,
                    min(retention_time),
                    max(retention_time)
                )
                new_scales <- sample(
                    peak_data$scale, num_new_events, replace = TRUE
                )
                new_powers <- sample_from_kde(power_density, num_new_events)

                # Combine all events
                gaussian_means <- c(sampled_times, new_times)
                gaussian_widths <- c(sampled_scales, new_scales)
                gaussian_amplitudes <- sqrt(c(sampled_powers, new_powers))

                chemical_noise <- numeric(length(retention_time))
                for (j in seq_along(gaussian_means)) {
                    chemical_noise <- chemical_noise +
                        generate_gaussian(
                            gaussian_means[j],
                            gaussian_amplitudes[j],
                            gaussian_widths[j],
                            retention_time
                        )
                }

                electronic_noise <- rtruncnorm(n_points,
                    a = 0, b = Inf,
                    mean = blank_noise_mean,
                    sd = blank_noise_std
                )

                # Total signal = chemical peaks + blank noise
                total_signal <- chemical_noise + electronic_noise

                noise_wavelet <- cwt(
                    total_signal, scales = scales, wavelet = wavelet_type
                )
                noise_wavelet[noise_wavelet < 0] <- 0
                simulated_power <- abs(noise_wavelet)^2

                local_counts <- local_counts + (simulated_power >= adjusted_power)
            }

            local_counts
        }
    })

    # Print time
    print(time_taken)

    # Clean up parallel backend
    stopImplicitCluster()

    # Compute p-values
    p_values <- counts / n_sim

    # p_values is [n_times x n_scales], with time points as rows and scales
    # as columns. To get minimum p-value across scales at each time point,
    # use apply() with MARGIN=1 (applies min to each row across its columns)
    p_values_min <- apply(p_values, 1, min) # One value per time point

    # Apply FWER correction to minimum p-values
    p_adjusted <- p.adjust(p_values_min, method = "holm") + .Machine$double.eps

    # Set significance level
    alpha <- 0.05

    # Save the current iteration's outputs into the simulation_results list
    simulation_results[[i]] <- list(
        instrument     = row$instrument,
        sample         = row$sample,
        expected_rt    = row$quant_rt,
        alpha          = alpha,
        n_sim          = n_sim,
        analyte_name   = "Clonidine",
        chrom_type     = "Quant",
        precursor_mz   = "232.2",
        product_mz     = "215.1",
        time           = retention_time,
        intensity      = sample_y,
        p_adjusted     = p_adjusted,
        adjusted_power = adjusted_power,
        scales         = scales
    )
}

# Convert the list to JSON
json_text <- toJSON(simulation_results, pretty = TRUE, auto_unbox = TRUE)

# Write the JSON text to a file
write(json_text, file = file.path("results", "json", "Clonidine_quant_results.json"))
