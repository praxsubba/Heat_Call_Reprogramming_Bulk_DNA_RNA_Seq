# Bioinformatic Pipeline for Analyzing Development Reprogramming in Heat Call Exposed Embryos

## Epigenetic Reprogramming Analysis: Complete and Annotated Workflow

**This pipeline is designed for robust, reproducible analysis of RRBS and RNA-seq data, integrating advanced normalization, network analysis, differential methylation, differential gene expression, gene set enrichment, and cell-type deconvolution in a single workflow.**

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

## RNA-Seq Differential Expression Analysis (Supplemental R Pipeline)

This repository also contains a **fully annotated and reproducible R script** for RNA-seq differential expression and gene set enrichment analysis, including both global and filtered (e.g., hypothalamic) gene expression analyses. The pipeline is modular, avoids hard-coded file paths, and is structured for maximum transparency and reproducibility.


### Analysis Steps

1. **Load Required Libraries:** Includes all necessary R packages for data import, analysis, and visualization.
2. **Define File Paths and Sample Information:** Generic file paths for transcript-to-gene mapping, quantification files, sample metadata, and gene set files.
3. **Read and Prepare Data:** Transcript-to-gene mapping, quantification file import, and sample metadata processing.
4. **Import Quantification Data:** Uses `tximport` for efficient import of Salmon quantification files.
5. **Differential Expression Analysis:** Constructs DESeq2 objects, filters low-count genes, and runs differential expression.
6. **Results Extraction and Shrinkage:** Extracts and shrinks log2 fold changes for robust gene ranking.
7. **Gene Set Enrichment Analysis:** Loads gene sets (Hallmark, KEGG, GO), prepares indices, and runs enrichment tests.
8. **Hypothalamic Filtered Analysis:** Optionally filters genes to a hypothalamic gene list and repeats differential expression and enrichment.
9. **Optional PCA for Hypothalamic Analysis:** Performs variance stabilizing transformation and PCA for sample-level visualization.

---

