library(Seurat)
library(scCustomize)
library(ggplot2)
library(tidyverse)
library(xlsx)

obj <- readRDS("data_objects/01_batch_integrated.RDS")

DefaultAssay(obj) <- "SCT"

dir.create("plots/02_cluster_id")
plots <- "plots/02_cluster_id"

#Cluster on integrated data with a low resolution for initial clustering
DefaultAssay(obj) <- "integrated"
obj <- FindClusters(obj, resolution = 0.2)

for (i in 0:18){
  markers <- FindConservedMarkers(obj,
                                  ident.1 = i,
                                  grouping.var = "Group",
                                  only.pos = T)
  write.xlsx(markers, file = paste0(plots_output, "batch_k5_conserved_markers.xlsx"),
             sheetName = paste0("cluster_", i), append = T)
}

obj <- PrepSCTFindMarkers(obj)
all_markers <- FindAllMarkers(obj)
saveRDS(all_markers, paste0(plots_output, "batchK5_all_markers.RDS"))
write.xlsx(all_markers, paste0(plots_output, "batchK5_all_markers.xlsx"))
saveRDS(obj, paste0(plots_output, "full_integrated_byBatch_k5_2_obj.RDS"))