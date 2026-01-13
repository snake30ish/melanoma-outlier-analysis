# Melanoma Survival Analysis & Outlier Detection

This repository contains the R code developed for the identification and analysis of biological and clinical outliers in cutaneous melanoma patients, using data from the TCGA-SKCM cohort.

## Overview

The analysis follows a multi-layer framework to detect heterogeneity at different levels:
1. **Biological Level**: Identification of functional pathway outliers using **Mahalanobis Distance** and **k-Nearest Neighbors (kNN)** with bootstrap-derived thresholds.
2. **Semi-parametric Level**: Detection of clinical outliers and influential observations in **Cox Proportional Hazards** models using Martingale residuals and DFBETA values.
3. **Stochastic Level**: Identification of stochastic outliers through **Probability Integral Transform (PIT)** residuals in parametric FHT/AFT (Log-Logistic) models.

## Repository Structure

* `analysis_script.R`: The main R script containing the full pipeline (preprocessing, modeling, and outlier detection).
* `pathway_outliers_18_clinical_profile.csv`: Detailed clinical profiles of the 18 identified biological outliers.
* `cox_outliers_profile.csv`: Data on patients identified as outliers or influential cases in the Cox model.
* `fht_outliers_profile.csv`: Patients identified as stochastic outliers (Early Deaths) via PIT residuals.
* clinical_plus_melanoma_pathway_scores_hospital_first.csv: the sample file of 460 patients with demographic and pathways scores
  

## Methodology

The pipeline integrates:
- **Functional Pathways**: Scores for 10 key oncogenic pathways (MAPK, Cell Cycle, Immune Checkpoint, etc.).
- **Survival Modeling**: Comparison between semi-parametric (Cox) and parametric (Log-Logistic) approaches.
- **Influence Analysis**: Leave-One-Out (LOO) diagnostics to ensure model robustness.

## Requirements

To run the analysis, the following R packages are required:
```R
install.packages(c("MASS", "FNN", "dplyr", "survival", "flexsurv", "survminer", "VennDiagram"))
