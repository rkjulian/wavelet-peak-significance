# Shared analysis functions for RvT4 wavelet analysis scripts
#
# Sourced by RvT4 analysis scripts (null-wavelet and chemical-noise).

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

# Calculate RT offset between a reference compound peak and internal standard
#
# Parameters:
#   data              - combined_data tibble with instrument/sample/is_rt columns
#   instrument        - instrument identifier (e.g. "6500", "7500")
#   chromatogram_data - chromatogram tibble with RT/Intensity/Compound/Sample/instrument columns
#   compound          - compound to use as RT reference (e.g. "Qual" or "Quant")
#   std_name          - sample name for the standard (e.g. "STD" or "6500 STD")
calculate_quant_offset <- function(data, instrument, chromatogram_data, compound, std_name) {
  quant_rt <- chromatogram_data |>
    filter(instrument == !!instrument, Sample == std_name, Compound == compound) |>
    slice(which.max(Intensity)) |>
    pull(RT)

  is_rt <- data |>
    filter(instrument == !!instrument, sample == std_name) |>
    pull(is_rt)

  quant_rt - is_rt
}

# Sampling function for extracting values from kernel density estimate
sample_from_kde <- function(density, n) {
  sample(density$x, size = n, prob = density$y, replace = TRUE)
}

# Function to create Gaussians for chemical noise simulation
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

# Detect peaks in wavelet power data using prominence-based algorithm
#
# Parameters:
#   data                    - numeric vector of wavelet power values
#   t                       - time vector corresponding to data
#   scale                   - wavelet scale at which peaks are being detected
#   min_power               - minimum power threshold for peak candidates
#   min_prominence_absolute - minimum prominence for a peak to be retained
#   merge_min_separation    - minimum time separation for merging nearby peaks
find_peaks <- function(data, t, scale, min_power = 1e-5,
                       min_prominence_absolute = 1e4,
                       merge_min_separation = 0.03) {
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
    scale = rep(scale, length(peaks))
  )

  # Add within-scale merging
  if (nrow(peak_data) > 0) {
    peak_data <- merge_nearby_peaks(peak_data, min_separation = merge_min_separation)
  }

  return(peak_data)
}

# Merge nearby peaks within the same wavelet scale to avoid duplicates
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
