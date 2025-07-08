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

### RNA-Seq Differential Expression and Gene Set Enrichment Pipeline

### 21. Setup and Data Preparation
**Input:**  
Raw RNA-seq quantification files (Salmon output), sample metadata, transcript-to-gene mapping files, gene set files (Hallmark, KEGG, GO), hypothalamic gene list (optional)

**Tasks:**
- Load required R packages  
- Set up working directory (not hard-coded for reproducibility)  
- Load sample metadata and prepare sample information  
- Prepare transcript-to-gene mapping files  
- Prepare gene set files for enrichment analysis

### 22. Data Import and Preprocessing
**Tools:**  
`tximport`, `DESeq2`, `limma`, `qusage`, `dplyr`, `tibble`, `readr`, `ggplot2`, `apeglm`, `ggsci`

**Tasks:**
- Import RNA-seq quantification data using tximport  
- Combine sample metadata with quantification data  
- Prepare sample information for differential expression analysis  
- Optionally filter gene list for tissue-specific (e.g., hypothalamic) analysis

### 23. Differential Expression Analysis
**Tasks:**
- Construct DESeq2 object from imported data  
- Filter low-count genes  
- Perform differential expression analysis with DESeq2  
- Extract and shrink log2 fold changes for robust gene ranking  
- Save results for downstream analysis

### 24. Gene Set Enrichment Analysis
**Tasks:**
- Load gene set files (Hallmark, KEGG, GO)  
- Prepare indices for gene set testing  
- Run gene set enrichment tests (cameraPR)  
- Optionally repeat for filtered gene lists (e.g., hypothalamic genes)

### 25. Quality Control and Visualization
**Tasks:**
- Visualize sample relationships with PCA (optional)  
- Generate scree plots for variance explained by principal components (optional)  
- Save results and plots for publication

---
