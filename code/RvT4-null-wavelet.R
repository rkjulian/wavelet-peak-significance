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
  "6500" = list(
    file_list = c(
      "MW-RvT4ex3-plasma-MS4-MF-blank.mzML",
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
      "6500 Blank", "6500 STD",
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
      "JD_ACP_ApoE_RvT_052022-blank.mzML",
      "JD_ACP_ApoE_RvT_052022-Std mix.mzML",
      "JD_ACP_ApoE_RvT_052022-Sample 1.mzML",
      "JD_ACP_ApoE_RvT_052022-Sample 2.mzML",
      "JD_ACP_ApoE_RvT_052022-Sample 3.mzML",
      "JD_ACP_ApoE_RvT_052022-Sample 4.mzML",
      "JD_ACP_ApoE_RvT_052022-Sample 5.mzML",
      "JD_ACP_ApoE_RvT_052022-Sample 6.mzML"
    ),
    sample_list = c(
      "7500 Blank", "7500 STD",
      "Chow 1", "Chow 2", "Chow 3",
      "WD 1", "WD 2", "WD 3"
    ),
    quant_index = c(34, 99, 34, 34, 34, 34, 34, 34),
    qual_index = c(33, 98, 33, 33, 33, 33, 33, 34),
    is_index = c(17, 32, 17, 17, 17, 17, 17, 17),
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
data_dir <- file.path("data", "S-BSST880")
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
quant_offset_6500 <- calculate_quant_offset(combined_data, "6500", chromatogram_data, "Qual", paste("6500", "STD"))
quant_offset_7500 <- calculate_quant_offset(combined_data, "7500", chromatogram_data, "Qual", paste("7500", "STD"))

combined_data <- combined_data |>
  mutate(quant_rt = case_when(
    instrument == "6500" ~ is_rt + quant_offset_6500,
    instrument == "7500" ~ is_rt + quant_offset_7500,
    TRUE ~ NA_real_
  ))

# Initialize an empty list to store results for each file
simulation_results <- list()

# Loop through combined data
for (i in 1:nrow(combined_data)) {
  row <- combined_data[i, ]

  instrument_name <- row$instrument
  sample_name <- row$sample
  expected_rt <- row$quant_rt

  if (instrument_name == "6500") {
    # Electronic noise from 6500 blank (external characterization):
    # Mean = 945.96, SD = 428.6769
    # Shapiro-Wilk p-value 0.065

    blank_noise_mean <- 945.96
    blank_noise_std <- 428.6769

    # Load chemical noise statistics (from RvT4-chemical-noise.R)
    peak_data <- read_csv(file.path("results", "tables", "RvT4_peak_data_6500.csv"))

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

    cat("  6500: lambda_events =", round(lambda_events, 2),
        "(from", nrow(events_per_sample), "samples)\n")

  } else if (instrument_name == "7500") {
    # Electronic noise from 7500 blank (external characterization):
    # Mean = 3567.7, SD = 429.622
    # Shapiro-Wilk p-value 0.16

    blank_noise_mean <- 3567.7
    blank_noise_std <- 429.622

    # Load chemical noise statistics (from RvT4-chemical-noise.R)
    peak_data <- read_csv(file.path("results", "tables", "RvT4_peak_data_7500.csv"))

    # Derive components needed for constructing null model
    power_density <- density(peak_data$power, from = 0) # KDE, truncated at zero
    power_density$y <- power_density$y / sum(power_density$y) # Renormalize

    # Compute lambda for Poisson event count from per-sample peak counts
    events_per_sample <- peak_data |>
      group_by(Sample) |>
      summarise(n_events = n(), .groups = "drop")
    lambda_events <- mean(events_per_sample$n_events)

    # Number of additional random chemical noise events per null simulation
    num_new_events <- 1

    cat("  7500: lambda_events =", round(lambda_events, 2),
        "(from", nrow(events_per_sample), "samples)\n")

  } else {
    print("Unknown Instrument. Stopping.")
    break
  }

  cat("Sample: ", sample_name, "\n")

  # get the time axis
  retention_time <- chromatogram_data[chromatogram_data$instrument == instrument_name &
    chromatogram_data$Sample == sample_name, ] |>
    filter(Compound == "Quant") |>
    pull(RT)

  # get the sample intensity values
  sample_y <- chromatogram_data[chromatogram_data$instrument == instrument_name &
    chromatogram_data$Sample == sample_name, ] |>
    filter(Compound == "Quant") |>
    pull(Intensity)

  # save time series length
  n_points <- length(retention_time)

  # --------------------------------------------------------------------------
  # 2. Wavelet Analysis and Monte Carlo Simulation
  # --------------------------------------------------------------------------

  # Parameters for CWT

  # MassSpecWavelet uses the mexh() function which is evaluated over [-8,8]
  # the maximum scale is set by the Nyquist criteria, or the limit of
  # this implementation (maximum scale is 64). The Nyquist criteria requires
  # two points per scale value so the minumum width is 16 and scale is
  # (n_points-1)/16

  scales <- seq(1, min((n_points - 1) / 16, 64), length.out = 100)
  wavelet_type <- "mexh" # Mexican hat wavelet

  # Compute CWT for the sample signal
  wavelet_result <- cwt(sample_y, scales = scales, wavelet = wavelet_type)

  # Truncate negative coefficients
  wavelet_result[wavelet_result < 0] <- 0

  # compute the observed power in sample
  adjusted_power <- abs(wavelet_result)^2

  # Setup for parallel processing using fork via mclapply backend
  # registerDoParallel(cores=) uses fork-based mclapply internally,
  # avoiding PSOCK serialization and explicit cluster management
  n_cores <- detectCores() - 4
  registerDoParallel(cores = n_cores)

  # Define parameters
  n_sim <- 100000 # Total simulations
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

  # p_values is [n_times x n_scales], with time points as rows and scales as columns
  # To get minimum p-value across scales at each time point, use apply() with MARGIN=1,
  # which applies min to each row (across its columns, i.e., across scales for that time)
  p_values_min <- apply(p_values, 1, min) # One value per time point

  # Apply FWER correction to minimum p-values
  p_adjusted <- p.adjust(p_values_min, method = "holm") + .Machine$double.eps

  # Set significance level
  alpha <- 0.05

  # Save the current iteration's outputs into the simulation_results list.
  simulation_results[[i]] <- list(
    instrument     = row$instrument,
    sample         = row$sample,
    expected_rt    = row$quant_rt,
    alpha          = alpha,
    n_sim          = n_sim,
    analyte_name   = "RvT4",
    precursor_mz   = "361.1",
    product_mz     = "211.1",
    time           = retention_time,
    intensity      = sample_y,
    p_adjusted     = p_adjusted,
    adjusted_power = adjusted_power,
    scales         = scales
  )
}

# Convert the list to JSON.
# The 'auto_unbox = TRUE' option makes sure that single values remain as scalars in JSON.
json_text <- toJSON(simulation_results, pretty = TRUE, auto_unbox = TRUE)

# Write the JSON text to a file
write(json_text, file = file.path("results", "json", "RvT4_results.json"))
