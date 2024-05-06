# Spatial Trancriptomic Analysis of Visium 10X Genomics data

This is the GitHub repository associated with my master's thesis "Spatial characterization of the dynamic changes driven by an exon6 deletion variant in murine CTRB1" for the master in Bioinformatics and Computational Biology at the Autonomus University of Madrid.

This work was carried out at the Epithelial Carcinogenesis Group at the CNIO, under the supervision of Francisco X. Real, Jaime Martinez de Villarreal and Enrique Carrillo.


<summary>Table of Contents</summary>
  <ol>
    <li><a href="#Installation">Installation</a></li>
    <li><a href="#Data availability">Scripts</a></li>
    <li><a href="#Scripts">Scripts</a></li>
    <li><a href="#Results">Results</a></li>
  </ol>



## Installation

- Clone the repository
- Create conda environment with necessary packages: this can be done by cloning the conda environment env.yml

```bash
conda env create -f env.yml
```

## Package versions

For this analysis R version xx was used.


## Data availability

Data used for this project were generated in the Epithelial Carcinogenesis group at the CNIO. 
These were Visium 10X Genomics Spatial Transcriptomics data generated with panccreatic samples of a GEMM.
Data is not publicly available in this repository.

## Scripts

1. Quality control

For quality control of the files FASTQC tool was used with the following command:

```bash
fastqc ./data/*.fastq.gz -o quality/ 
```

2. Spaceranger software

Spaceranger software was run on each capture area (1_spaceranger.sh)

4. Downstream analysis with R
   - Merge and clustering (Seurat, Harmony)
   - Signature analysis (VISION)
   - Single cell integration (RCTD)
   - Correlations and DE analysis
   - Spot deconvolution (BayesSpace)
   - Regulon Analysis (SCENIC)

## Results
