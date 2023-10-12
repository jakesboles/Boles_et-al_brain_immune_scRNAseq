library(Seurat)
library(scCustomize)
library(ggplot2)
library(tidyverse)
library(xlsx)

obj <- readRDS("data_objects/01_batch_integrated.RDS")

DefaultAssay(obj) <- "SCT"

dir.create("plots/02_cluster_id")
plots <- "plots/02_cluster_id/"

stats <- "tabular_output/"

#Cluster on integrated data with a low resolution for initial clustering ----
DefaultAssay(obj) <- "integrated"
obj <- FindClusters(obj, resolution = 0.2)

DefaultAssay(obj) <- "SCT"

#Look at canonical cell markers to help with annotation ----
DimPlot_scCustom(obj)
ggsave(paste0(plots, "unlabeled_clusters.png"), units = "in", dpi = 600,
       height = 8, width = 8)

#Endothelial cells
FeaturePlot_scCustom(obj, 
                     features = c("Cldn5", "Ly6c1", "Podxl"),
                     label = T,
                     num_columns = 3,
                     slot = "data",
                     na_cutoff = 1) #12
ggsave(paste0(plots, "endothelial_markers.png"), units = "in",
       dpi = 600, height = 5, width = 15)

#IEGs to look at activated cells
FeaturePlot_scCustom(obj, 
                     features = c("Jun", "Fos"), 
                     label = T,
                     num_columns = 1,
                     slot = "data",
                     na_cutoff = 1) #LPS cells and 9
ggsave(paste0(plots, "jun-fos_umap.png"), units = "in",
       dpi = 600, height = 5, width = 10)

#NK cells and T-cells
FeaturePlot_scCustom(obj, 
                     features = c("Nkg7", "Cd3e", "Cd3g", "Trbc1"),
                     label = T,
                     num_columns = 2,
                     slot = "data",
                     na_cutoff = 1) #13
ggsave(paste0(plots, "t_nk.png"), units = "in",
       dpi = 600, height = 10, width = 10)

#Astrocytes 
FeaturePlot_scCustom(obj, 
                     features = c("Gja1", "Gfap", "Aqp4"), 
                     label = T,
                     num_columns = 3,
                     slot = "data",
                     na_cutoff = 1) #17
ggsave(paste0(plots, "astrocyte.png"), units = "in",
       dpi = 600, height = 5, width = 10)

#Microglia
FeaturePlot_scCustom(obj, 
                     features = c("P2ry12", "Hexb", "Tmem119", "Aif1"),
                     label = T,
                     num_columns = 2,
                     slot = "data",
                     na_cutoff = 1)
ggsave(paste0(plots, "microglia.png"), units = "in",
       dpi = 600, height = 10, width = 10)
FeaturePlot_scCustom(obj, 
                     features = c("Ngp", "S100a8", "S100a9", "Ly6g"),
                     label = T,
                     num_columns = 2,
                     slot = "data",
                     na_cutoff = 1) #11
ggsave(paste0(plots, "neutrophil.png"), units = "in",
       dpi = 600, height = 10, width = 10)

#Oligos
FeaturePlot_scCustom(obj, 
                     features = c("Olig1", "Olig2"), 
                     label = T,
                     num_columns = 2,
                     na_cutoff = 1,
                     slot = "data") #22
ggsave(paste0(plots, "oligodendrocytes.png"), units = "in",
       dpi = 600, height = 10, width = 10)

#Choroid plexus
FeaturePlot_scCustom(obj, 
                     features = c("Ttr", "Folr1", "Prlr", "Aqp1"), 
                     slot = "data",
                     label = T,
                     na_cutoff = 1,
                     num_columns = 2)
ggsave(paste0(plots_path, "choroid_markers.png"), dpi = 600,
       units = "in", height = 10, width = 10)

#Pericytes
FeaturePlot_scCustom(obj, 
                     features = c("Rgs5", "Acta2", "Pdgfrb", "Des"), 
                     slot = "data", 
                     label = T,
                     na_cutoff = 1,
                     num_columns = 2)
ggsave(paste0(plots_path, "pericyte_markers.png"), dpi = 600,
       units = "in", height = 10, width = 10)

#Meningeal/leptomeningeal cells
FeaturePlot_scCustom(obj, 
                     features = c("Dcn", "Igfbp2", "Foxc1", "Slc47a1"), 
                     slot = "data",
                     label = T,
                     na_cutoff = 1,
                     num_columns = 2)
ggsave(paste0(plots_path, "meningeal_markers.png"), dpi = 600,
       units = "in", height = 10, width = 10)

#Neurons
FeaturePlot_scCustom(obj, 
                     features = c("Slc17a7", "Gad2", "Tubb3"), 
                     label = T,
                     num_columns = 3,
                     slot = "data",
                     na_cutoff = 1)
ggsave(paste0(plots, "neurons.png"), units = "in",
       dpi = 600, height = 5, width = 15)

#Some monocyte/macrophage markers
FeaturePlot_scCustom(obj, 
                     features = c("Cd14", "Ccr2", "Cd163", "Mrc1", "Ly6c2", "Arg1"), 
                     label = T,
                     num_columns = 3,
                     slot = "data",
                     na_cutoff = 1)
ggsave(paste0(plots, "monocytes_macrophages.png"), units = "in",
       dpi = 600, height = 10, width = 15)

#CD45 and CD11b
FeaturePlot_scCustom(obj, 
                     features = c("Ptprc", "Itgam"), 
                     label = T,
                     num_columns = 2,
                     slot = "data",
                     na_cutoff = 1)
ggsave(paste0(plots, "cd45-cd11b.png"), units = "in",
       dpi = 600, height = 5, width = 10)

#compute conserved markers for each group ----
for (i in 0:18){
  markers <- FindConservedMarkers(obj,
                                  ident.1 = i,
                                  grouping.var = "Group",
                                  only.pos = T)
  write.xlsx(markers, file = paste0(stats, "cluster_id_1_conserved_markers.xlsx"),
             sheetName = paste0("cluster_", i), append = T)
}

#Find ALL markers to have a complementary source of markers ----
obj <- PrepSCTFindMarkers(obj)
all_markers <- FindAllMarkers(obj)

#Save objects and markers as needed ----
saveRDS(all_markers, "data_objects/cluster_id_1_all_markers.RDS")
write.xlsx(all_markers, paste0(stats, "cluster_id_1_all_markers.xlsx"))

#Add annotations ----
cluster_id <- c("Microglia1", "Microglia2", "Microglia3", "Microglia4", "Microglia5",
                "Monocyte/Macrophage1", "Microglia6", "Monocyte/Macrophage2", "Monocyte/Macrophage3", 
                "Neutrophil", "Monocyte/Macrophage4", "B-cell", "NK/T-cell",
                "Endothelial/Vascular1", "Choroid plexus", "Meninges/CNS", "Astrocyte", "Endothelial/Vascular2",
                "Microglia7")

names(cluster_id) <- levels(obj)
obj <- RenameIdents(obj, cluster_id)
DimPlot_scCustom(obj, 
                 label = T, 
                 label.box = F,
                 repel = T) +
  NoLegend()
ggsave(paste0(plots, "annotated_UMAP.png"), units = "in",
       dpi = 600, height = 10, width = 10)

saveRDS(obj, "data_objects/01_batch_integrated.RDS")