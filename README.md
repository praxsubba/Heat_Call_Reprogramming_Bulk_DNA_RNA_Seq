# Bioinformatic Pipeline for Analyzing Development Reprogramming in Heat Call Exposed Embryos

## Epigenetic Reprogramming Analysis: Complete and Annotated Workflow

This repository provides a well-organized, thoroughly annotated, and reproducible R workflow for Reduced Representation Bisulfite Sequencing, RNA-Sequencing of brain tissue in heat call exposed embryos. The code is structured for clarity and ready for submission to GitHub as a supplement for publication in _____  

The pipeline covers the following key stages:


1. **Setup and Data Preparation**
2. **Normalization and Quality Control**
3. **Network Construction and Module Detection**
4. **Module-Trait Relationships**
5. **Gene Ontology Enrichment**
6. **Hub Gene Identification and Validation**
7. **Network Visualization (Green Module)**
8. **Cell-Type Enrichment Analysis**
9. **Single-Cell RNA-seq Deconvolution Signature Preparation**
10. **Seurat Marker Gene Identification**
11. **SLURM Batch Script for CIBERSORTx Deconvolution**

---

## Pipeline Overview

### 1. Setup and Data Preparation
- **Input:** Raw RNA-seq count data
- **Tasks:** Load data, process gene length information, prepare metadata

### 2. Normalization and Quality Control
- **Tools:** DESeq2 VST normalization
- **Tasks:** Normalize counts, perform sample QC, filter outliers

### 3. Network Construction and Module Detection
- **Tools:** WGCNA
- **Tasks:** Construct co-expression network, detect gene modules

### 4. Module-Trait Relationships
- **Tasks:** Correlate modules with clinical/experimental traits

### 5. Gene Ontology Enrichment
- **Tasks:** Enrichment analysis for detected modules

### 6. Hub Gene Identification and Validation
- **Tasks:** Identify and validate hub genes within modules

### 7. Network Visualization (Green Module)
- **Tasks:** Visualize the co-expression network for selected module(s)

### 8. Cell-Type Enrichment Analysis
- **Tasks:** Assess module enrichment for specific cell types

### 9. Single-Cell RNA-seq Deconvolution Signature Preparation
- **Tasks:** Prepare signatures for deconvolution using single-cell data

### 10. Seurat Marker Gene Identification
- **Tasks:** Identify marker genes using Seurat (if applicable)

### 11. SLURM Batch Script for CIBERSORTx Deconvolution
- **Tasks:** Run CIBERSORTx deconvolution on HPC using provided SLURM script

---

**This pipeline is designed for robust, reproducible analysis of RNA-seq data, integrating advanced normalization, network analysis, and cell-type deconvolution in a single workflow.**

