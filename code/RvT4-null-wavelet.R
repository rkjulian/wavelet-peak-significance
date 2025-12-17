library(MSnbase)
library(tidyverse)
library(MassSpecWavelet)
library(viridis)
library(foreach)
library(doParallel)
library(truncnorm)
library(jsonlite)

# Clear the environment
rm(list = ls())

# Set seed for reproducibility
set.seed(42)

PLOT_COI <- FALSE

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

# -----------------------------------------------------------------------------
# 1. Read raw and meta data
# -----------------------------------------------------------------------------

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
data_dir <- "large-data/S-BSST880"
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
calculate_quant_offset <- function(data, instrument) {
  quant_rt <- chromatogram_data |>
    filter(instrument == !!instrument, Sample == paste(instrument, "STD"), Compound == "Qual") |>
    slice(which.max(Intensity)) |>
    pull(RT)

  is_rt <- data |>
    filter(instrument == !!instrument, sample == paste(instrument, "STD")) |>
    pull(is_rt)

  quant_rt - is_rt
}

quant_offset_6500 <- calculate_quant_offset(combined_data, "6500")
quant_offset_7500 <- calculate_quant_offset(combined_data, "7500")

combined_data <- combined_data |>
  mutate(quant_rt = case_when(
    instrument == "6500" ~ is_rt + quant_offset_6500,
    instrument == "7500" ~ is_rt + quant_offset_7500,
    TRUE ~ NA_real_
  ))

# Sampling function for extracting values from power density
sample_from_kde <- function(density, n) {
  sample(density$x, size = n, prob = density$y, replace = TRUE)
}

# Function to create Guassians for chemical noise
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

# Setup the PDF file
pdf(file.path("analysis", "RvT4-quant-stats.pdf"),
  width = 11, height = 8.5, onefile = TRUE
)

par(mfrow = c(3, 1), mar = c(4, 4, 2, 2))

# Initialize an empty list to store results for each file
simulation_results <- list()

# Loop through combined data
for (i in 1:nrow(combined_data)) {
  row <- combined_data[i, ]

  instrument_name <- row$instrument
  sample_name <- row$sample
  expected_rt <- row$quant_rt

  if (instrument_name == "6500") {
    # Data from the 6500 blank sample:
    # Mean = 945.96
    # SD = 428.6769
    # Shapiro-Wilk p-value 0.065
    # Bootstrap SD bias: -5.230671
    # Blank number of events (mean) 5 (sd ~ 1)

    blank_noise_mean <- 945.96
    blank_noise_std <- 428.6769

    num_existing_events <- 5 # Number of existing events to sample
    num_new_events <- 1 # Number of new chemical noise events to generate

    # Load chemical noise statistics (from RvT4-chemical-noise.R)
    peak_data <- read_csv(file.path("analysis", "RvT4_peak_data_6500.csv"))

    # Derive components needed for constructing null model

    power_density <- density(peak_data$power, from = 0) # KDE, truncated at zero
    power_density$y <- power_density$y / sum(power_density$y) # Renormalize
  } else if (instrument_name == "7500") {
    # Data from the 7500 blank sample
    # W = 0.97932, p-value = 0.1568
    # Mean = 3567.7
    # SD = 429.622
    # Shapiro-Wilk p-value 0.16
    #
    # Blank number of events (mean) 10 (max 13) (sd ~ 1)

    blank_noise_mean <- 3567.7
    blank_noise_std <- 429.622

    num_existing_events <- 10 # Number of existing events to sample
    num_new_events <- 1 # Number of new chemical noise events to generate

    # Load chemical noise statistics (from chemical-noise.R)
    peak_data <- read_csv(file.path("analysis", "RvT4_peak_data_7500.csv"))

    # Derive components needed for constructing null model

    power_density <- density(peak_data$power, from = 0) # KDE, truncated at zero
    power_density$y <- power_density$y / sum(power_density$y) # Renormalize
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

  # setup for parallel processing: a cluster of cores using sockets
  cl <- makeCluster(detectCores() - 1, type = "PSOCK")
  registerDoParallel(cl)

  # Load required packages on workers
  invisible(clusterEvalQ(cl, {
    library(MassSpecWavelet) # needed for cwt() inside cluster
    library(truncnorm) # needed for the truncated rnorm
  }))


  # Export data objects needed in the loop
  clusterExport(cl, c(
    "peak_data", # characteristics of chemical noise
    "scales", # used in cwt()
    "wavelet_type", # used in cwt()
    "retention_time", # the time axis of the sample
    "adjusted_power", # observed power from sample
    "power_density", # power density distribution of chemical noise
    "sample_from_kde", # function to get a power value from density
    "generate_gaussian", # function to generate chemical noise peaks
    "num_existing_events", # number of chemical noise events per sample
    "num_new_events", # number of random chemical noise events
    "blank_noise_mean", # mean of normal distribution for blank noise
    "blank_noise_std" # std dev of blank noise
  ))

  # Define parameters
  n_sim <- 1000000 # Total simulations
  chunk_size <- 100000 # Number of simulations per chunk

  # Initialize the count matrix
  counts <- matrix(0, nrow = nrow(adjusted_power), ncol = ncol(adjusted_power))

  time_taken <- system.time(
    # Run in chunks
    for (start_sim in seq(1, n_sim, by = chunk_size)) {
      end_sim <- min(start_sim + chunk_size - 1, n_sim)
      cat("Processing simulations", start_sim, "to", end_sim, "\n")

      # Initialize chunk_counts to zero for each chunk
      chunk_counts <- matrix(0, nrow = nrow(adjusted_power), ncol = ncol(adjusted_power))

      # Run a chunk of simulations in parallel
      chunk_counts <- foreach(sim = start_sim:end_sim, .combine = "+") %dopar% {
        # Sequential loop over simulations (for debugging)
        # for (sim in start_sim:end_sim) {

        sampled_times <- sample(peak_data$time, num_existing_events)
        sampled_scales <- sample(peak_data$scale, num_existing_events)
        sampled_powers <- sample_from_kde(power_density, num_existing_events)

        # Add two new Gaussian events
        new_times <- runif(num_new_events, min(retention_time), max(retention_time)) # Uniformly sample times
        new_scales <- sample(peak_data$scale, num_new_events, replace = TRUE) # Sample scales from empirical distribution
        new_powers <- sample_from_kde(power_density, num_new_events) # Sample powers from KDE

        # Combine all events
        gaussian_means <- c(sampled_times, new_times)
        gaussian_widths <- c(sampled_scales, new_scales)
        gaussian_amplitudes <- sqrt(c(
          sampled_powers,
          new_powers
        )) # A = sqrt(P_wavelet)

        chemical_noise <- numeric(length(retention_time))
        for (i in seq_along(gaussian_means)) {
          chemical_noise <- chemical_noise +
            generate_gaussian(
              gaussian_means[i],
              gaussian_amplitudes[i],
              gaussian_widths[i],
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

        noise_wavelet <- cwt(total_signal, scales = scales, wavelet = wavelet_type)
        noise_wavelet[noise_wavelet < 0] <- 0
        simulated_power <- abs(noise_wavelet)^2

        # Compute the p-value: proportion of simulations exceeding observed power
        # count_exceeding <- simulated_power >= adjusted_power

        # Accumulate results
        # chunk_counts <- chunk_counts + count_exceeding


        return(simulated_power >= adjusted_power)
      }

      # Accumulate results from the chunk
      counts <- counts + chunk_counts
    }
  )

  # print time
  print(time_taken)

  # Clean up cluster
  stopCluster(cl)

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

  # Identify significant time points
  significant_times <- which(p_adjusted < alpha)

  # --------------------------------------------------------------------------
  # 3. Plots
  # --------------------------------------------------------------------------

  # Plot 1: Original signals with SNR
  plot(retention_time, sample_y,
    type = "l", col = "blue",
    main = paste0("Sample: ", sample_name),
    xlab = "Time", ylab = "Amplitude",
    xaxs = "i"
  )
  abline(v = expected_rt, col = "black", lty = 2)

  # Plot 2: Log-scaled Wavelet power spectrum
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

  # Plot 3: Minimum p-values with FWER correction
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

  # Add COI shading as full-height rectangles
  if (PLOT_COI) {
    # Calculate COI width for maximum scale
    dt <- mean(diff(retention_time))
    max_scale <- max(scales)
    coi_width_time <- 5 * max_scale * dt

    # Left edge
    rect(retention_time[1], y_min, retention_time[1] + coi_width_time, y_max,
      col = rgb(0.5, 0.5, 0.5, 0.3), border = NA
    )
    # Right edge
    rect(retention_time[n_points] - coi_width_time, y_min, retention_time[n_points], y_max,
      col = rgb(0.5, 0.5, 0.5, 0.3), border = NA
    )
  }

  abline(h = -log(alpha), col = "red", lty = 2)
  points(retention_time[significant_times],
    -log(p_adjusted[significant_times]),
    col = "red", pch = 20
  )
  abline(v = expected_rt, col = "black", lty = 2)

  # Save the current iteration's outputs into the simulation_results list.
  simulation_results[[i]] <- list(
    instrument     = row$instrument,
    sample         = row$sample,
    expected_rt    = row$quant_rt,
    alpha          = alpha,
    time           = retention_time,
    intensity      = sample_y,
    p_adjusted     = p_adjusted,
    adjusted_power = adjusted_power,
    scales         = scales
  )
}
dev.off()

# Convert the list to JSON.
# The 'auto_unbox = TRUE' option makes sure that single values remain as scalars in JSON.
json_text <- toJSON(simulation_results, pretty = TRUE, auto_unbox = TRUE)

# Write the JSON text to a file
write(json_text, file = "RvT4_results.json")
