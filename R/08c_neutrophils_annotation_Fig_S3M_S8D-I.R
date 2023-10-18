library(Seurat)
library(scCustomize)
library(patchwork)
library(dittoSeq)
library(tidyverse)
library(paletteer)

obj <- readRDS("data_objects/05_clean_neutrophils_sub.RDS")

#Need to run PrepSCTFindMarkers because the data were integrated in 05
obj <- PrepSCTFindMarkers(obj)