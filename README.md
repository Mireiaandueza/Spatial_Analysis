# Spatial Trancriptomic Analysis of Visium 10X Genomics data

This is the GitHub repository associated with my master's thesis "Spatial characterization of the dynamic changes driven by an exon6 deletion variant in murine CTRB1" for the master in Bioinformatics and Computational Biology at the Autonomus University of Madrid.


<summary>Table of Contents</summary>
  <ol>
    <li><a href="#Installation">Installation</a></li>
    <li><a href="#Data availability">Scripts</a></li>
    <li><a href="#Scripts">Scripts</a></li>
    <li><a href="#Results">Results</a></li>
  </ol>



## Installation

- Clone the repository
- Create conda environment with necessary packages: this can be done by cloning the following conda environments:

```bash
conda env create -f env.yml
```

## Package versions

For this analysis R version 4.3.3 and 4.3.2 were used.


## Data availability

Data used for this project were generated in the Epithelial Carcinogenesis group at the CNIO. 
These were Visium 10X Genomics Spatial Transcriptomics data generated with pancreatic samples of a GEMM.
Data is not publicly available in this repository.

## Scripts

### 1. Quality control

For quality control of the files FASTQC tool was used with the following command:

```bash
fastqc ./data/*.fastq.gz -o quality/ 
```

### 2. Spaceranger software

Spaceranger software was run on each capture area (1_spaceranger.sh)

```bash
./1_spaceranger.sh
```

### 3. Downstream analysis with R
  - Create Seurat objects and add metadata from mannualy annotated CSVs
```{R}
2_add_metadata.R
```
  - Merge and clustering (Seurat, Harmony)
```{R}
3_merge.R # Pipeline
harmony.R # Harmony trial
```
  - Signature analysis  (VISION) and create figures
```{R}
4_vision.R
5_vision_figures.R
```    
  - Single cell integration (RCTD and FindAnchors)
```{R}
6_1_singlecell_integration.R
6_2_RCTD.R
6_3_RCTD_join.R
```  
  - Correlations and DE analysis
    

  - Spot deconvolution (BayesSpace)
```{R}
enhancing_BayesSpace.R
```  
  - Regulon Analysis (SCENIC)
```{R}
scenic.R
```
  - Other helper functions are also added


