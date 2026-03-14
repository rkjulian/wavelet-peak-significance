# waters-noise-pooled.R
# Compute descriptive statistics for Waters Tox reagent blank noise
# by pooling intensity values across all three blanks for each
# compound/chrom transition before computing statistics.
#
# Run from project root:
#   Rscript code/waters-noise-pooled.R

library(tidyverse)
library(MSnbase)
library(NADA)

# Clear environment
rm(list = ls())

# --- Configuration -----------------------------------------------------------

reagent_blanks <- c(
    "20Feb24 reagent blank 1.mzML",
    "20Feb24 reagent blank 2.mzML",
    "20Feb24 reagent blank 3.mzML"
)

# Transition table from docs/waters-tox-experiment.md
transitions <- tibble(
    Compound = c("Gabapentin", "Gabapentin",
                 "Clonidine",  "Clonidine",
                 "Lorazepam",  "Lorazepam"),
    Chrom    = c("Qual", "Quant",
                 "Qual", "Quant",
                 "Qual", "Quant"),
    Q1       = c(172.3, 172.3,
                 232.2, 232.2,
                 321.2, 321.2),
    Q3       = c(137.3, 154.3,
                 44,    215.1,
                 229.1, 275),
    index    = c(1L, 2L, 7L, 8L, 10L, 11L)
)

# Time window for noise extraction (minutes)
t_min <- 0.6
t_max <- 3.0

# --- Collect raw intensities from all blanks ---------------------------------

# One list entry per transition, each a numeric vector pooled across blanks
pooled <- vector("list", nrow(transitions))
for (r in seq_len(nrow(transitions))) {
    pooled[[r]] <- numeric(0)
}

for (b in seq_along(reagent_blanks)) {
    file_name <- reagent_blanks[b]
    srm_filename <- file.path("data", "Waters_Tox", "mzML", file_name)

    cat("Reading", file_name, "\n")
    srm <- readSRMData(srm_filename)

    for (r in seq_len(nrow(transitions))) {
        idx <- transitions$index[r]

        t_raw <- rtime(srm[[idx]])
        y_raw <- intensity(srm[[idx]])

        keep <- t_raw >= t_min & t_raw <= t_max
        pooled[[r]] <- c(pooled[[r]], y_raw[keep])
    }
}

# --- Compute statistics on pooled data ---------------------------------------

results_list <- list()

for (r in seq_len(nrow(transitions))) {
    y <- pooled[[r]]

    # NADA ROS (no transformation)
    censored <- y <= 0
    ros_fit  <- ros(y, censored,
                    forwardT = NULL, reverseT = NULL)

    results_list[[r]] <- tibble(
        Compound = transitions$Compound[r],
        Chrom    = transitions$Chrom[r],
        Q1       = transitions$Q1[r],
        Q3       = transitions$Q3[r],
        N        = length(y),
        Mean     = mean(y),
        SD       = sd(y),
        Median   = median(y),
        MAD      = mad(y),
        ROS_Mean = mean(ros_fit),
        ROS_SD   = sd(ros_fit),
        Censored = sum(censored) / length(y)
    )
}

# --- Save --------------------------------------------------------------------

summary_df <- bind_rows(results_list)

output_file <- file.path("results", "tables", "waters-noise-pooled.csv")
write_csv(summary_df, output_file)
cat("Saved", nrow(summary_df), "rows to", output_file, "\n")
