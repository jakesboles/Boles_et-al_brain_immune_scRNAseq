# Boles_et-al_brain_immune_scRNAseq

<sup>**This repository is currently under construction**</sup>


### This repository contains code to reproduce analyses and plots from:
[**Microfluidics-free single-cell genomics reveals complex central-peripheral immune crosstalk in the mouse brain during peripheral inflammation** (*BioRxiv*)](https://www.biorxiv.org/content/10.1101/2023.10.05.561054v1)

Jake Sondag Boles<sup>1</sup>, Oihane Uriarte Huarte, & Malú Gámez Tansey, 2023

<sup><sup>1</sup> Analysis lead and contact (jake.boles@ufl.edu)</sup>

## Original single-cell RNA sequencing data:
Our data can be accessed via the NCBI GEO ([GSE245309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE245309)). In this GEO submission, we have included the raw FASTQ files and the three output files from the [PIPseeker pipeline](https://www.fluentbio.com/products/pipseeker-software-for-data-analysis/) (version 2). The script `00_prep_project_and_read_data.R` will load these three output files for each sample and assemble them into `R` objects. 

## Code:
This project was run in `R` (v4.2-4.3). This repository will load version-controlled packages using the `renv` package. Certain key analyses, including normalization and integration of data, were done through the University of Florida's HiPerGator 3.0 cluster with up to 24 CPUs and 185GB RAM. If you intend to run this analysis yourself, I would recommend loading this repository onto your institution's high-performance cluster.

To get started, open a new R session and run:
```
if (!require(usethis)) {
  install.packages("usethis")
}

library(usethis)

usethis::use_course(
  'jakesboles/Boles_et-al_brain_immune_scRNAseq',
  destdir = 'PATH/TO/DIRECTORY')
```

## Relevant external resources:
### PIPseq:
### PIPseeker:
### Seurat:
### The chooseR clustering procedure:
### The cross entropy test for differences in dimensionally-reduced cell embeddings:
### scCustomize: 
### Artifactual gene expression modules:
### Psuedo-bulked differential expression with DESeq2:
### Pseudo-bulked gene set enrichment analyses with clusterProfiler:
### Curated lists of genes that describe the heritability of Alzheimer's and Parkinson's diseases:
