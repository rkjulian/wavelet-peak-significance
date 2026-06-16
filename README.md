# wavelet-peak-significance

Supporting code, data, and figures for:

> **Monte Carlo Wavelet Analysis for Objective Peak Detection in SRM LC-MS/MS Analysis**
>
> Randall K. Julian, Jr.<sup>1,2</sup>\*, Brian A. Rappold<sup>3,4</sup>, Fan Yi<sup>5</sup>, Stephen R. Master<sup>5,6,7</sup>
>
> 1. Indigo BioAutomation, Inc., Carmel, IN, USA
> 2. Department of Chemistry, Purdue University, West Lafayette, IN, USA
> 3. Laboratory Corporation of America Holdings, Research Triangle Park, NC, USA
> 4. School of Health Sciences, University of Iceland, Reykjavik, Iceland
> 5. Center for Diagnostic Innovation (CDI), Children's Hospital of Philadelphia
> 6. Department of Pathology and Laboratory Medicine, Children's Hospital of Philadelphia
> 7. Perelman School of Medicine, University of Pennsylvania
>
> \*Email: rkjulian@indigobio.com

## Overview

This repository contains all source data, analysis code, pre-computed results, and supplemental figures for the paper. The method introduces a wavelet-based Monte Carlo framework that empirically characterizes chemical noise in SRM LC-MS/MS chromatograms, builds a generative null model, and assigns family-wise error rate (FWER) controlled p-values for peak detection.

Three datasets are included:

- **Validation dilution series (Waters Xevo TQ-XS)**: Clonidine, Gabapentin, and Lorazepam spiked into human plasma at nine two-fold dilution levels (0.0195--5.00 ng/mL), each in triplicate
- **RvT4 lipid mediator reanalysis (SCIEX 6500 and 7500)**: Data from Walker et al. (2024), BioStudies accession S-BSST880
- **Ketamine external validation (SCIEX 6500)**: Clinical pain panel assay, 16 samples including a low calibration standard, matrix blank, and 14 biological samples (4 confirmed sub-LLOQ positives)

## Repository Structure

```
code/                        # R scripts for simulation, analysis, and plotting
data/
  S-BSST880/mzML/           # RvT4 mzML files (SCIEX 6500 and 7500)
  Waters_Tox/mzML/           # Dilution series mzML files (Waters Xevo TQ-XS)
  ketamine/mzML/             # Ketamine mzML files (SCIEX 6500)
results/
  json/                      # Pre-computed wavelet analysis results (JSON)
  tables/                    # Peak data and noise estimates (CSV)
figures/
  main_figures/              # Assembled figures for manuscript (see below)
  supplemental_figures/      # Per-sample analysis PDFs (see below)
```

## Figures

All figures are pre-generated PDFs organized into two directories. If you are looking for a specific plot, start here.

### Main Figures (`figures/main_figures/`)

Assembled multi-panel figures used in the manuscript text. Each file is a single-page PDF.

**Null model example chromatograms** -- Show how the generative null model produces synthetic chromatograms that preserve the chemical noise structure of the original data. Each panel overlays the observed chromatogram with several null-model realizations drawn from the fitted noise parameters.

| File | What it shows |
|------|---------------|
| `RvT4-null-chromatograms.pdf` | RvT4 observed vs. null-model chromatograms for the SCIEX 6500 and 7500 instruments. Demonstrates that the null model faithfully reproduces the baseline noise structure on both platforms. |
| `waters-null-chromatograms.pdf` | Lorazepam observed vs. null-model chromatograms from the Waters Xevo TQ-XS dilution series. Shows null model fidelity across concentration levels from blank through high-concentration samples. |

**Chemical noise time histograms** -- Histograms of peak retention times extracted from reagent blanks and low-concentration samples. Peaks that cluster at specific retention times indicate chemical noise (reproducible interferents) rather than random electronic noise.

| File | What it shows |
|------|---------------|
| `RvT4-chemical-noise-time-histograms.pdf` | Distribution of chemical noise peak locations across RvT4 blank and plasma samples on the SCIEX 6500 and 7500. Reveals instrument-specific interferent patterns. |
| `waters-chemical-noise-time-histograms.pdf` | Distribution of chemical noise peak locations for Clonidine, Gabapentin, and Lorazepam blanks on the Waters platform. Shows compound-specific interferent retention times. |

**Chemical noise power distributions** -- Kernel density estimates (KDEs) of wavelet power at scales corresponding to chemical noise peaks. These characterize the amplitude distribution of interferent signals that the null model must reproduce.

| File | What it shows |
|------|---------------|
| `RvT4-chemical-noise-power-distributions.pdf` | KDE of chemical noise peak power for RvT4 on SCIEX 6500 and 7500. Compares the power distributions between instruments. |
| `waters-chemical-noise-power-distributions.pdf` | KDE of chemical noise peak power for Clonidine, Gabapentin, and Lorazepam on the Waters platform. Compares power distributions across compounds. |

**Combined chemical noise comparison** -- Side-by-side comparison of time histograms and power distributions for the Waters dilution series compounds on a single page.

| File | What it shows |
|------|---------------|
| `waters-chemical-noise-combined.pdf` | Combined time histogram and power KDE panels for all three Waters compounds (Clonidine, Gabapentin, Lorazepam), enabling direct cross-compound comparison of chemical noise characteristics. |

### Supplemental Figures (`figures/supplemental_figures/`)

Multi-page PDFs with one page per sample. These provide the complete per-sample analysis for every sample in both datasets.

**Hypothesis testing (wavelet significance plots)** -- Each page is a four-panel figure for one sample showing: (1) the raw SRM chromatogram, (2) a CWT power scalogram (time vs. wavelet scale, colored by power), (3) the Monte Carlo null power envelope at the analyte scale, and (4) the FWER-adjusted p-value time series with a significance threshold line.

| File | Samples included |
|------|------------------|
| `Clonidine-quant-plots.pdf` | Clonidine dilution series: reagent blanks and nine concentration levels (0.0195X--64X), each in triplicate |
| `Gabapentin-quant-plots.pdf` | Gabapentin dilution series: reagent blanks and nine concentration levels (0.0195X--64X), each in triplicate |
| `Lorazepam-quant-plots.pdf` | Lorazepam dilution series: reagent blanks and nine concentration levels (0.0195X--64X), each in triplicate |
| `RvT4-quant-plots.pdf` | RvT4 on SCIEX 6500 and 7500: reagent blanks, standard mixes, and biological plasma samples (chow-diet and Western-diet mice) |
| `Ketamine-quant-plots.pdf` | Ketamine on SCIEX 6500: low calibration standard, matrix blank, and 14 biological samples (4 sub-LLOQ positives, 9 negatives, 1 chemical noise case) |

**Chemical noise characterization** -- Each page shows one sample's chemical noise analysis: (1) the chromatogram with detected chemical noise peaks marked, (2) a histogram of peak retention times, (3) a histogram of peak widths (wavelet scales), and (4) a histogram of peak powers. Together these panels define the empirical distributions used to parameterize the generative null model for that sample's transition.

| File | Samples included |
|------|------------------|
| `Clonidine_Quant_Chemical_Noise_Analysis.pdf` | All Clonidine reagent blanks and dilution levels |
| `Gabapentin_Quant_Chemical_Noise_Analysis.pdf` | All Gabapentin reagent blanks and dilution levels |
| `Lorazepam_Quant_Chemical_Noise_Analysis.pdf` | All Lorazepam reagent blanks and dilution levels |
| `RvT4_Chemical_Noise_Analysis.pdf` | All RvT4 blanks and biological samples on both SCIEX instruments |
| `Ketamine_Quant_Chemical_Noise_Analysis.pdf` | All 16 ketamine samples on SCIEX 6500 |

## Pre-Computed Results

Results in `results/json/` are the output of Monte Carlo simulations and serve as input to the plotting scripts. They contain CWT coefficients, p-values, peak detection calls, and chemical noise characterizations.

| File | Description |
|------|-------------|
| `{Compound}_quant_results.json` | Dilution series hypothesis testing (Clonidine, Gabapentin, Lorazepam) |
| `{Compound}_quant_chemical_noise.json` | Dilution series chemical noise models |
| `RvT4_results.json` | RvT4 hypothesis testing results |
| `RvT4_chemical_noise.json` | RvT4 chemical noise model |
| `Ketamine_quant_results.json` | Ketamine hypothesis testing results |
| `Ketamine_quant_chemical_noise.json` | Ketamine chemical noise model |

Peak data tables in `results/tables/` contain extracted peak metrics (retention time, area, width, p-value) for each sample:

| File | Description |
|------|-------------|
| `{Compound}_quant_peak_data_waters.csv` | Dilution series peaks (Clonidine, Gabapentin, Lorazepam) |
| `waters-noise-pooled.csv` | Pooled electronic noise estimates from reagent blanks |
| `RvT4_peak_data_6500.csv` | RvT4 peaks from SCIEX 6500 |
| `RvT4_peak_data_7500.csv` | RvT4 peaks from SCIEX 7500 |
| `Ketamine_quant_peak_data_sciex.csv` | Ketamine peaks from SCIEX 6500 |

## Running the Workflow

All R scripts assume the working directory is the repository root.

### Requirements

- R (>= 4.0)
- CRAN packages: `tidyverse`, `viridis`, `jsonlite`, `truncnorm`, `patchwork`, `ggh4x`, `forcats`, `foreach`, `doParallel`, `NADA`
- Bioconductor packages: `MSnbase` (mzML reading), `MassSpecWavelet` (CWT)

Install Bioconductor packages with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("MSnbase", "MassSpecWavelet"))
```

### Regenerate Figures from Pre-Computed Results

The plotting scripts read JSON files from `results/json/` and produce figure PDFs. No simulation is needed.

```bash
# Main figures -- null model chromatograms
Rscript code/plot-RvT4-null-chromatograms.R
Rscript code/plot-waters-null-chromatograms.R

# Main figures -- chemical noise comparison plots
Rscript code/plot-RvT4-chemical-noise.R
Rscript code/waters-chemical-noise-comparison-plots.R

# Supplemental figures -- hypothesis testing (per-sample)
Rscript code/waters-null-wavelet-plot.R
Rscript code/RvT4-null-wavelet-plot.R
Rscript code/Ketamine-null-wavelet-plot.R

# Supplemental figures -- chemical noise characterization (per-sample)
Rscript code/Clonidine-quant-chemical-noise.R
Rscript code/Gabapentin-quant-chemical-noise.R
Rscript code/Lorazepam-quant-chemical-noise.R
Rscript code/RvT4-chemical-noise.R
Rscript code/Ketamine-quant-chemical-noise.R

# Supplemental figures -- chemical noise (plot-only from JSON)
Rscript code/plot-waters-chemical-noise-supplemental.R
Rscript code/plot-RvT4-chemical-noise-supplemental.R
Rscript code/plot-Ketamine-chemical-noise-supplemental.R
```

### Re-Run Monte Carlo Simulations

Simulation scripts read mzML data, run Monte Carlo wavelet analysis, and write JSON results.

```bash
# Waters dilution series
Rscript code/Clonidine-quant-null-wavelet.R
Rscript code/Gabapentin-quant-null-wavelet.R
Rscript code/Lorazepam-quant-null-wavelet.R

# RvT4
Rscript code/RvT4-null-wavelet.R

# Ketamine
Rscript code/Ketamine-quant-null-wavelet.R
```

### Code Organization

| Script | Purpose |
|--------|---------|
| `mc_wavelet_functions.R` | Shared analysis functions (CWT, null model, Monte Carlo) |
| `mc_wavelet_plot_functions.R` | Shared plotting functions |
| `*-null-wavelet.R` | Monte Carlo simulation (mzML -> JSON) |
| `*-null-wavelet-plot.R` | Hypothesis testing plots (JSON -> supplemental PDF) |
| `*-chemical-noise.R` | Chemical noise characterization (mzML -> JSON + supplemental PDF) |
| `plot-*-chemical-noise-supplemental.R` | Chemical noise plots from pre-computed JSON (plot-only) |
| `plot-RvT4-null-chromatograms.R` | RvT4 null model chromatograms (main figure) |
| `plot-waters-null-chromatograms.R` | Waters null model chromatograms (main figure) |
| `plot-RvT4-chemical-noise.R` | RvT4 chemical noise comparison (main figures) |
| `waters-chemical-noise-comparison-plots.R` | Waters chemical noise comparison (main figures) |
| `plot-RvT4-chemical-noise-supplemental.R` | RvT4 chemical noise supplemental comparison |
| `plot-waters-chemical-noise-supplemental.R` | Waters chemical noise supplemental comparison |
| `waters-noise-pooled.R` | Pooled noise estimates across replicates |

## Data Sources

### Validation Dilution Series

Clonidine, Gabapentin, and Lorazepam reference standards spiked into human K2EDTA plasma. Serial two-fold dilutions from 5.00 ng/mL to 0.0195 ng/mL, prepared in triplicate. Analyzed on Waters Acquity I-Class Plus UPLC coupled with Xevo TQ Absolute XR. See manuscript Methods for full sample preparation and LC-MS/MS conditions.

### RvT4 Lipid Mediator

Raw data from Walker et al. (2024) obtained from BioStudies database accession [S-BSST880](https://www.ebi.ac.uk/biostudies/studies/S-BSST880). Plasma samples analyzed on SCIEX 6500 and 7500 instruments. SRM transition 361.1 -> 211.1 (quantifier).

> Walker, M.E., De Matteis, R., Perretti, M., Dalli, J., 2024. Resolvin T4 enhances macrophage cholesterol efflux to reduce vascular disease. Nat Commun 15, 975. https://doi.org/10.1038/s41467-024-44868-1

### Ketamine External Validation

Clinical pain panel assay acquired on a SCIEX 6500 instrument. SRM transition m/z 238 -> 125 (quantifier). Sixteen samples analyzed: one low calibration standard (STD), one matrix blank, and 14 biological samples. Four biological samples (Inj 8, 11, 13, 27) are confirmed ketamine-positive with signal intensity below the low calibration standard. See manuscript Methods for full assay details.

## Platform Compatibility

This code was developed and tested on macOS and Linux. It should also work on Windows, but Windows has not been extensively tested. If you encounter any platform-specific issues, please open an issue on this repository or contact rkjulian@indigobio.com.

The Monte Carlo simulation scripts use `foreach` with `doParallel` for parallel processing. On macOS and Linux, `registerDoParallel(cores=)` uses fork-based parallelism; on Windows it falls back to PSOCK clusters. The `.packages` and `.export` arguments in each `foreach()` call ensure that required packages and variables are available to worker processes on all platforms.

## License

Copyright 2025 The Authors. Licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE) for the full text.
