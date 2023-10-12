library(Seurat)
library(scCustomize)
library(dittoSeq)

dir.create("plots/04_cluster_id")
plots <- "plots/04_cluster_id/"

obj <- readRDS("data_objects/03_cleaned.RDS")

obj <- PrepSCTFindMarkers(obj)

DefaultAssay(obj) <- "SCT"

markers <- FindAllMarkers(obj)
write.csv(markers,
          file = "tabular_output/cleaned_data_markers.csv")

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

#Initial DimPlot for Supp. Fig. 3 ----
theme <- theme(axis.text = element_blank(),
               axis.ticks = element_blank(),
               axis.line = element_line(arrow = arrow(angle = 15, 
                                                      length = unit(0.5, "cm"), 
                                                      type = "closed")))

DefaultAssay(obj) <- "SCT"
Idents(obj) <- obj$seurat_clusters

DimPlot_scCustom(obj, label = T) + 
  NoLegend() +
  theme
ggsave(paste0(plots, "Fig_S3A.png"),
       units = "in", dpi = 600,
       height = 4, width = 4)

#Cell markers for Supp. Fig. 3 ----
markerplot <- function(feature){
  p <- FeaturePlot_scCustom(obj, features = feature,
                            pt.size = 0.1, raster = F, label = F,
                            na_cutoff = 1) +
    theme(plot.title = element_text(face = "bold.italic")) +
    theme
  
  return(p)
}

markersave <- function(gene){
  ggsave(filename = paste0(plots, "Fig_S3_", gene, ".png"),
         units = "in", dpi = 600,
         height = 3, width = 3.5)
}

markerplot("Ptprc")
markersave("ptprc")

markerplot("Itgam")
markersave("itgam")

markerplot("Cd3e")
markersave("cd3e")

markerplot("Nkg7")
markersave("nkg7")

markerplot("Ms4a1")
markersave("ms4a1")

markerplot("Ngp")
markersave("ngp")

markerplot("P2ry12")
markersave("p2ry12")

markerplot("Aif1")
markersave("aif1")

markerplot("Tmem119")
markersave("tmem119")

markerplot("Hexb")
markersave("hexb")

#Annotate ----
cluster_id <- c("Microglia0", "Microglia1", "Microglia2", "Microglia3", "Microglia4",
                "Microglia5", "Monocytes/Macrophages6", "Microglia7", "Microglia8",
                "Monocytes/Macrophages9", "Neutrophils10", "Monocytes/Macrophages11",
                "T/NK cells12", "Monocytes/Macrophages13", "Monocytes/Macrophages14",
                "Microglia15", "B-cells16", "Endothelial cells17", "Microglia18",
                "Choroid plexus19", "Astrocytes20", "Endothelial cells21", 
                "Leptomeningeal cells22", "Pericytes23", "Oligodendrocytes24")
names(cluster_id) <- levels(obj)

obj <- RenameIdents(obj, cluster_id)

obj$cluster_id <- Idents(obj)

obj@meta.data <- obj@meta.data %>%
  dplyr::mutate(cell_type = case_when(
    str_detect(cluster_id, regex("Microglia")) ~ "Microglia",
    str_detect(cluster_id, regex("Monocytes")) ~ "Monocytes/Macrophages",
    str_detect(cluster_id, regex("Neut")) ~ "Neutrophils",
    str_detect(cluster_id, regex("Endothelial")) ~ "Endothelial cells",
    str_detect(cluster_id, regex("Astro")) ~ "Astrocytes",
    str_detect(cluster_id, regex("Lepto")) ~ "Leptomeningeal cells",
    str_detect(cluster_id, regex("Pericytes")) ~ "Pericytes",
    str_detect(cluster_id, regex("T/NK")) ~ "T/NK cells",
    str_detect(cluster_id, regex("B-")) ~ "B-cells",
    str_detect(cluster_id, regex("Choroid")) ~ "Choroid plexus",
    str_detect(cluster_id, regex("Oligo")) ~ "Oligodendrocytes"))

Idents(obj) <- "cell_type"

obj@meta.data <- obj@meta.data %>%
  dplyr::mutate(sfig3_idents = if_else(cell_type == "Astrocytes" |
                                         cell_type == "Endothelial cells" | 
                                         cell_type == "Oligodendrocytes" |
                                         cell_type == "Choroid plexus" |
                                         cell_type == "Pericytes" | 
                                         cell_type == "Leptomeningeal cells",
                                       "CD45-", cell_type)) %>%
  mutate(sfig3_idents = if_else(cell_type == "Monocytes/Macrophages",
                                "Myeloid APCs", sfig3_idents)) %>%
  mutate(sfig3_idents = factor(sfig3_idents,
                               levels = sort(unique(sfig3_idents))))
Idents(obj) <- obj$sfig3_idents

pal <- dittoColors()

DimPlot_scCustom(obj, label = F,
                 colors_use = pal[19:24]) +
  theme
ggsave(paste0(plots, "Fig_S3H.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)