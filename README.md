# wavelet-peak-significance

Supporting code and data for wavelet-peak-significance method described in "Monte Carlo Wavelet Analysis for Objective Peak Detection in LC-MS/MS"

## Abstract

Detection and quantification of low-level analytes in complex chromatographic-mass spectrometric data ultimately depend on objective criteria for deciding when an apparent peak is distinguishable from background. Conventional signal-to-noise and LOD/LOQ rules are typically justified under simple, constant-variance (white) noise models and do not control the probability of obtaining at least one spurious peak across a chromatogram, even though structured chemical noise and co-eluting interferences dominate at trace levels. We introduce a wavelet-based Monte Carlo framework for single reaction monitoring (SRM) LC-MS/MS traces that empirically characterizes chemical-noise peaks, builds a generative noise-only null model directly from these observations, and assigns per-time-point $p$-values with family-wise error control. The need for such objective, statistically grounded criteria has become particularly acute in studies of analytes near the limits of detection (LOD). As a motivating example, we apply our framework to lipid mediator research, illustrating how chemically realistic null models and chromatogram-level error control can alter calls near the LOD while preserving robust detection of established prostaglandin-positive controls. At an adjusted $p < 0.05$ (controlling the chance of at least one false-positive peak per chromatogram at about $5\%$), our criteria show strong concordance with the presence or absence of confirmatory product ions in this dataset, suggesting that our empirical noise model captures the main features of chemically structured background peaks.

## Authors

**Randall K. Julian, Jr.**<sup>1,2</sup>\*  
**Brian A. Rappold**<sup>3,4</sup>  
**Stephen R. Master**<sup>5,6</sup>

### Affiliations

1. Indigo BioAutomation, Inc., Carmel, IN, USA
2. Department of Chemistry, Purdue University, West Lafayette, IN, USA
3. Laboratory Corporation of America Holdings, Research Triangle Park, NC, USA
4. School of Health Sciences, University of Iceland, Reykjavik, Iceland
5. Department of Pathology and Laboratory Medicine, Children's Hospital of Philadelphia
6. Perelman School of Medicine, University of Pennsylvania

\*Email: rkjulian@indigobio.com
