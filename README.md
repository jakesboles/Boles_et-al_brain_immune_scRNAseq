# Boles_et-al_brain_immune_scRNAseq

**This repository is currently under construction**


### This repository contains code to reproduce analyses and plots from:
[**Microfluidics-free single-cell genomics reveals complex central-peripheral immune crosstalk in the mouse brain during peripheral inflammation** (*BioRxiv*)](https://www.biorxiv.org/content/10.1101/2023.10.05.561054v1)

Jake Sondag Boles<sup>1</sup>, Oihane Uriarte Huarte, & Malú Gámez Tansey, 2023

<sup><sup>1</sup> Analysis lead and contact (jake.boles@ufl.edu)</sup>

## Original single-cell RNA sequencing data:

## Code and relevant external resources:
This project was run in R v4.2-4.3. This repository will load version-controlled packages using the `renv` package. Certain key analyses, including normalization and integration of data, was done through the University of Florida's HiPerGator cluster. 

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
