# SEM Texture Analysis of Lignocellulosic Substrates

This repository provides an R pipeline for quantitative texture analysis of
scanning electron microscopy (SEM) images of lignocellulosic substrates.

The workflow extracts multiple regions of interest (ROIs) along a transect
and computes texture descriptors to compare **initial** and **residual**
substrates after fungal cultivation.

## Features
- Transect-based ROI sampling
- Automatic scale calibration
- Texture descriptors:
  - Laplacian variance
  - Otsu bright pixel percentage
  - GLCM contrast
  - GLCM homogeneity
  - GLCM energy (ASM)
  - GLCM entropy
- Paired initial vs residual comparison
- Publication-ready figures and tables

## Requirements
R ≥ 4.2  
Packages: EBImage, ggplot2, dplyr, tidyr, boot, patchwork

## Usage
1. Place SEM images in `data/SEM/`
2. Edit image paths in the script
3. Run `scripts/SEM_transect_ROI_pipeline.R`
4. Results are saved in `outputs/` and `figures/`

## License
Creative Commons Attribution–NonCommercial 4.0 International (CC BY-NC 4.0)
