# qPCR Data Analysis and Normalization with R

## Overview

This project presents a modular R-based pipeline for analyzing and normalizing **quantitative PCR (qPCR)** data, tailored for **gene expression studies**—particularly those investigating **phosphorus metabolism** (e.g., comparisons of high vs. low phosphorus groups). The pipeline performs data cleaning, outlier detection, replicate merging, and housekeeping gene-based normalization to generate clean, analyzable results.

---

## Dataset and Tools

- **Input File**: `qPCR_Result.csv` – Exported qPCR data from LightCycler 480, including:
  - `Experiment.Name`, `SampleName`, `CP`, `Concentration`
  
- **Tools & Libraries**:
  - R (v4.0 or later)
  - Required packages:
    - `openxlsx`
    - `outliers`
    - `ggplot2`

---

## Goals

- Automate **qPCR data cleaning** and **replicate merging**.
- Detect and optionally remove **outliers** using Grubbs’ test.
- Normalize gene expression values using **housekeeping genes**.
- Export **publication-ready tables** and **visual summaries**.

---

## Repository Structure

```
PCR_Phosphorus_Low_High/
├── qPCR_Result.csv                  # Raw qPCR input data
├── PCR_Low.csv                      # Optional: log2-transformed group means
├── R/
│   ├── PCR_Normalization_Script.R   # Main R script for pipeline
│   ├── Data_Normalized.xlsx         # Final normalized output
```

---

## Reproducibility Instructions

1. **Set the working directory** in R:
```r
setwd("Z:/30Deser/user/hasan/Privat/PCR_Phosphorus_Low_High/R")
```

2. **Install required R packages**:
```r
install.packages("openxlsx")
install.packages("outliers")
install.packages("ggplot2")
```

3. **Customize script parameters**:
```r
datum <- c("2708021")
housekeeping_genes <- c("RPL32")
ErstesTranscript <- "CYP24A1"
norm_method <- "per_group"   # Options: "per_group" or "all"
normgroup <- "tissue"
outlier <- "n"               # Set to "y" to activate Grubbs’ test
```

4. **Run the main script**:
```r
source("PCR_Normalization_Script.R")
```

---

## Outputs

- `Data_Normalized.xlsx`: Final table of normalized expression values.
- `qPCR_VitD_per_group_tissue_RPL32.csv`: Flat format normalized data.
- `PCR_Low.csv`: Log₂-transformed mean values for groups (optional).
- `qPCR_outlier.txt`: Flags and details of detected outliers (if enabled).

---

## Key Features

- **Flexible Normalization**:
  - Normalize within groups (`"per_group"`) or across all samples (`"all"`).
  - Supports multiple housekeeping genes using the **geometric mean**.

- **Optional Outlier Filtering**:
  - Uses Grubbs’ test to detect technical replicate outliers.
  - Toggle on/off via `outlier` parameter.

- **Customizable Settings**:
  - Easily modify transcript targets, housekeeping genes, and group labels.

---

## Purpose

This qPCR normalization pipeline was developed as part of a gene expression study on **phosphorus metabolism** and is adaptable for other transcriptomic projects requiring robust and reproducible qPCR data handling.

---

## Author

**Maruf Hasan**  
Interests: Interests: Molecular aging | Bioinformatics | Translational health research
