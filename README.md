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

Two datasets are included:

- **Validation dilution series (Waters Xevo TQ-XS)**: Clonidine, Gabapentin, and Lorazepam spiked into human plasma at nine two-fold dilution levels (0.0195--5.00 ng/mL), each in triplicate
- **RvT4 lipid mediator reanalysis (SCIEX 6500 and 7500)**: Data from Walker et al. (2024), BioStudies accession S-BSST880

## Repository Structure

```
code/                        # R scripts for simulation, analysis, and plotting
data/
  S-BSST880/mzML/           # RvT4 mzML files (SCIEX 6500 and 7500)
  Waters_Tox/mzML/           # Dilution series mzML files (Waters Xevo TQ-XS)
results/
  json/                      # Pre-computed wavelet analysis results (JSON)
  tables/                    # Peak data and noise estimates (CSV)
figures/
  supplemental_figures/      # Per-sample analysis PDFs (see below)
```

## Supplemental Figures

The `figures/supplemental_figures/` directory contains multi-page PDFs with per-sample chromatograms, wavelet power maps, and significance test results for every sample in both datasets. These are referenced from the paper and its supplemental information.

**Hypothesis testing (wavelet significance plots):**

| File | Contents |
|------|----------|
| `Clonidine-quant-plots.pdf` | Clonidine dilution series (blanks through 64X, triplicate) |
| `Gabapentin-quant-plots.pdf` | Gabapentin dilution series (blanks through 64X, triplicate) |
| `Lorazepam-quant-plots.pdf` | Lorazepam dilution series (blanks through 64X, triplicate) |
| `RvT4-quant-plots.pdf` | RvT4 on SCIEX 6500 and 7500 (blanks, standards, plasma samples) |

Each page shows: raw chromatogram, CWT power scalogram, Monte Carlo null power, and adjusted p-value time series for one sample.

**Chemical noise characterization:**

| File | Contents |
|------|----------|
| `Clonidine_Quant_Chemical_Noise_Analysis.pdf` | Chemical noise characterization for Clonidine |
| `Gabapentin_Quant_Chemical_Noise_Analysis.pdf` | Chemical noise characterization for Gabapentin |
| `Lorazepam_Quant_Chemical_Noise_Analysis.pdf` | Chemical noise characterization for Lorazepam |
| `RvT4_Chemical_Noise_Analysis.pdf` | Chemical noise peak distributions for RvT4 samples |

Each page shows: observed chemical noise peaks, empirical distributions of peak location, width (wavelet scale), and power used to build the null model.

## Pre-Computed Results

Results in `results/json/` are the output of Monte Carlo simulations and serve as input to the plotting scripts. They contain CWT coefficients, p-values, peak detection calls, and chemical noise characterizations.

| File | Description |
|------|-------------|
| `{Compound}_quant_results.json` | Dilution series hypothesis testing (Clonidine, Gabapentin, Lorazepam) |
| `{Compound}_quant_chemical_noise.json` | Dilution series chemical noise models |
| `RvT4_results.json` | RvT4 hypothesis testing results |
| `RvT4_chemical_noise.json` | RvT4 chemical noise model |

Peak data tables in `results/tables/` contain extracted peak metrics (retention time, area, width, p-value) for each sample:

| File | Description |
|------|-------------|
| `{Compound}_quant_peak_data_waters.csv` | Dilution series peaks (Clonidine, Gabapentin, Lorazepam) |
| `waters-noise-pooled.csv` | Pooled electronic noise estimates from reagent blanks |
| `RvT4_peak_data_6500.csv` | RvT4 peaks from SCIEX 6500 |
| `RvT4_peak_data_7500.csv` | RvT4 peaks from SCIEX 7500 |

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

The plotting scripts read JSON files from `results/json/` and produce the supplemental figure PDFs. No simulation is needed.

```bash
# Chemical noise characterization
Rscript code/Clonidine-quant-chemical-noise.R
Rscript code/Gabapentin-quant-chemical-noise.R
Rscript code/Lorazepam-quant-chemical-noise.R
Rscript code/RvT4-chemical-noise.R

# Waters dilution series hypothesis testing plots
Rscript code/waters-null-wavelet-plot.R

# RvT4 hypothesis testing plots
Rscript code/RvT4-null-wavelet-plot.R
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
```

### Code Organization

| Script | Purpose |
|--------|---------|
| `mc_wavelet_functions.R` | Shared analysis functions (CWT, null model, Monte Carlo) |
| `mc_wavelet_plot_functions.R` | Shared plotting functions |
| `*-null-wavelet.R` | Monte Carlo simulation (mzML -> JSON) |
| `*-null-wavelet-plot.R` | Hypothesis testing plots (JSON -> PDF) |
| `*-chemical-noise.R` | Chemical noise characterization (JSON -> PDF) |
| `waters-noise-pooled.R` | Pooled noise estimates across replicates |

## Data Sources

### Validation Dilution Series

Clonidine, Gabapentin, and Lorazepam reference standards spiked into human K2EDTA plasma. Serial two-fold dilutions from 5.00 ng/mL to 0.0195 ng/mL, prepared in triplicate. Analyzed on Waters Acquity I-Class Plus UPLC coupled with Xevo TQ Absolute XR. See manuscript Methods for full sample preparation and LC-MS/MS conditions.

### RvT4 Lipid Mediator

Raw data from Walker et al. (2024) obtained from BioStudies database accession [S-BSST880](https://www.ebi.ac.uk/biostudies/studies/S-BSST880). Plasma samples analyzed on SCIEX 6500 and 7500 instruments. SRM transition 361.1 -> 211.1 (quantifier).

> Walker, M.E., De Matteis, R., Perretti, M., Dalli, J., 2024. Resolvin T4 enhances macrophage cholesterol efflux to reduce vascular disease. Nat Commun 15, 975. https://doi.org/10.1038/s41467-024-44868-1

## License

See LICENSE file for details.
