library(tidyverse)
library(viridis)
library(jsonlite)

# Clear the environment
rm(list = ls())

source(file.path("code", "mc_wavelet_plot_functions.R"))

PLOT_COI <- FALSE

# -----------------------------------------------------------------------------
# List of result JSON files to process (one per compound).
# Add new compounds here as they are analyzed.
# -----------------------------------------------------------------------------

result_files <- c(
    file.path("results", "json", "RvT4_results.json")
)

# -----------------------------------------------------------------------------
# Process each compound
# -----------------------------------------------------------------------------

for (json_file in result_files) {
    cat("Loading", json_file, "...\n")

    # Read JSON file
    simulation_results <- fromJSON(json_file, simplifyVector = FALSE)

    # Extract compound metadata from first sample entry
    analyte_name <- simulation_results[[1]]$analyte_name
    precursor_mz <- simulation_results[[1]]$precursor_mz
    product_mz   <- simulation_results[[1]]$product_mz

    cat("Processing", analyte_name, "-", length(simulation_results), "samples\n")

    # Derive output PDF name from analyte name
    pdf_name <- paste0(analyte_name, "-quant-plots.pdf")
    pdf(file.path("figures", "supplemental_figures", pdf_name),
        width = 11, height = 8.5, onefile = TRUE)

    par(mfrow = c(3, 1), mar = c(4, 4, 2, 2))

    # Set significance level
    alpha <- 0.05

    # Loop through each result
    for (i in seq_along(simulation_results)) {
        result <- simulation_results[[i]]

        # Extract data from JSON
        instrument_name <- result$instrument
        sample_name <- result$sample
        expected_rt <- result$expected_rt
        n_sim <- result$n_sim
        retention_time <- unlist(result$time)
        sample_y <- unlist(result$intensity)
        p_adjusted <- unlist(result$p_adjusted)
        p_adjusted[p_adjusted <= 0] <- .Machine$double.eps
        adjusted_power <- do.call(rbind, lapply(result$adjusted_power, unlist))
        scales <- unlist(result$scales)

        n_points <- length(retention_time)

        # Identify significant time points
        significant_times <- which(p_adjusted < alpha)

        cat("Plotting sample:", sample_name, "\n")

        plot_null_wavelet_sample(retention_time, sample_y, adjusted_power,
            scales, p_adjusted, expected_rt, sample_name, alpha, n_sim,
            compound_name = "Quant",
            precursor_mz = precursor_mz,
            product_mz = product_mz,
            plot_coi = PLOT_COI,
            analyte_name = analyte_name)
    }

    dev.off()

    cat("Saved to figures/supplemental_figures/", pdf_name, "\n", sep = "")
}
