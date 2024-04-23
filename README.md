# Spatial Analysis

This is the GitHub repository associated with my master's thesis "Spatial characterization of the dynamic changes driven by the PDAC risk variant  Ctrb1 Δ exon6  in mouse pancreas" for the master in Bioinformatics and Compuational Biology at the Autonomus University of Madrid.

This work was carried out at the Epithelial Carcinogenesis group at the CNIO.

<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#prerequisites">Prerequisites</a></li>
    <li><a href="#data availability">Scripts</a></li>
    <li><a href="#scripts">Scripts</a></li>
    <li><a href="#results">Results</a></li>
  </ol>
</details>


## Installation

- Clone the repository
- Create conda environment with necessary packages: this can be done by cloning the conda environment env.yml

```bash
conda env create -f env.yml
```

## Package versions

For this analysis R version xx was used.


## Data Availability

Data used for this project were generated in the Epithelial Carcinogenesis group at the CNIO. 
These were Visium 10X Genomics Spatial Transcriptomics data generated with panccreatic samples of GEMM mouse models.
Data is not publicly available in this repository.

## Scripts

1. Quality control
2. Spaceranger software
3. Downstream analysis with R
   - Merge and clustering (Seurat, Harmony)
   - Signature analysis (VISION)
   - Single cell integration (RCTD)
   - Correlations
   - Spot deconvolution (BayesSpace)
   - Regulon Analysis (SCENIC)

