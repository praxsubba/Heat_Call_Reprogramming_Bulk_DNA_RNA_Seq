# Bioinformatic Pipeline for Analyzing Development Reprogramming in Heat Call Exposed Embryos

## Epigenetic Reprogramming Analysis: Complete and Annotated Workflow

This repository provides a well-organized, thoroughly annotated, and reproducible R workflow for Reduced Representation Bisulfite Sequencing (RRBS), RNA-Sequencing of brain tissue in heat call exposed embryos. The code is structured for clarity and ready for submission to GitHub as a supplement for publication in _____.  

The pipeline covers the following key stages:

---

## Pipeline Overview

### 1. Setup and Data Preparation
- **Input:** Raw RRBS methylation files, sample metadata
- **Tasks:** Load packages, set working directory, load sample metadata, prepare methylation file lists

### 2. Data Import and Preprocessing
- **Tools:** `methylKit`, `tidyverse`, `dplyr`
- **Tasks:** Import methylation data, add sample-specific covariates (e.g., sex), filter samples based on availability

### 3. Quality Control and Filtering
- **Tasks:** Visualize methylation and coverage statistics for all samples, filter samples by coverage to ensure data quality

### 4. Differential Methylation Analysis
- **Tasks:** Merge samples, perform differential methylation analysis (including sex as a covariate), extract significantly differentially methylated CpGs, analyze treatment effects separately for males and females

### 5. Summary and Comparison
- **Tasks:** Summarize and compare differential methylation results across sexes, report overlap of differentially methylated CpGs

### 6. Gene Annotation
- **Tools:** `GenomicRanges`
- **Tasks:** Annotate differentially methylated CpGs with gene information using genomic coordinates, save annotated results for manual review

### 7. Interaction Analysis (Sex × Treatment)
- **Tasks:** Create interaction factors, re-analyze methylation data for interaction effects, perform pairwise comparisons between treatment groups, annotate and summarize results

### 8. Additional Analyses and Visualization
- **Tasks:** Load and analyze specific methylation files, generate publication-ready boxplots and interaction plots for selected regions

### 9. Gene Ontology (GO) Enrichment Analysis
- **Tools:** `goseq`, `geneLenDataBase`, `org.Gg.eg.db`
- **Tasks:** Prepare gene universe for GO analysis, perform enrichment analysis, visualize top enriched GO terms

---

## RNA-Sequencing Analysis (Co-Expression Network and Deconvolution)

### 10. Setup and Data Preparation
- **Input:** Raw RNA-seq count data
- **Tasks:** Load data, process gene length information, prepare metadata

### 11. Normalization and Quality Control
- **Tools:** DESeq2 VST normalization
- **Tasks:** Normalize counts, perform sample QC, filter outliers

### 12. Network Construction and Module Detection
- **Tools:** WGCNA
- **Tasks:** Construct co-expression network, detect gene modules

### 13. Module-Trait Relationships
- **Tasks:** Correlate modules with clinical/experimental traits

### 14. Gene Ontology Enrichment
- **Tasks:** Enrichment analysis for detected modules

### 15. Hub Gene Identification and Validation
- **Tasks:** Identify and validate hub genes within modules

### 16. Network Visualization (Green Module)
- **Tasks:** Visualize the co-expression network for selected module(s)

### 17. Cell-Type Enrichment Analysis
- **Tasks:** Assess module enrichment for specific cell types

### 18. Single-Cell RNA-seq Deconvolution Signature Preparation
- **Tasks:** Prepare signatures for deconvolution using single-cell data

### 19. Seurat Marker Gene Identification
- **Tasks:** Identify marker genes using Seurat (if applicable)

### 20. SLURM Batch Script for CIBERSORTx Deconvolution
- **Tasks:** Run CIBERSORTx deconvolution on HPC using provided SLURM script

---

**This pipeline is designed for robust, reproducible analysis of RRBS and RNA-seq data, integrating advanced normalization, network analysis, differential methylation, and cell-type deconvolution in a single workflow.**
