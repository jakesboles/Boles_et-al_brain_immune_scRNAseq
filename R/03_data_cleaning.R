library(Seurat)
library(scCustomize)
library(tidyverse)
library(ggplot2)
library(glmGamPoi)

#In the last script, you may notice lots of cells that are, for example, P2ry12+
#but are amidst a cluster of cells that are clearly T-cells...
#This script goes cluster by cluster and removes what is likely microglial contamination.
#After this cleaning, the object is reassembled for the annotation that went into the paper.

obj <- readRDS("data_objects/02_annotated.RDS")

#The gene module below is derived from https://www.nature.com/articles/s41593-022-01022-8
marsh_microglia_id <- list(c("P2ry12", "Fcrls", "Trem2", "Tmem119", "Cx3cr1", "Hexb",
                             "Tgfbr1","Sparc","P2ry13","Olfml3","Adgrg1","C1qa",
                             "C1qb","C1qc","Csf1r","Fcgr3","Ly86", "Laptm5"))
obj <- AddModuleScore(obj, marsh_microglia_id, 
                      name = "marsh_microglia_id")