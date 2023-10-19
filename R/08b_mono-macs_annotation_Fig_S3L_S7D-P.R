library(Seurat)
library(scCustomize)
library(patchwork)
library(dittoSeq)
library(tidyverse)
library(paletteer)

obj <- readRDS("data_objects/05_clean_mono-macs_sub.RDS")

dir.create("plots/03_fine_cluster_id")
plots <- "plots/03_fine_cluster_id/"

#Need to run PrepSCTFindMarkers because the data were integrated in 05
obj <- PrepSCTFindMarkers(obj)

#Make unicolored plot for Fig S3 ----
pal <- dittoColors()
col <- pal[22]

theme <- theme(axis.text = element_blank(),
               axis.ticks = element_blank(),
               axis.line = element_line(arrow = arrow(angle = 15, 
                                                      length = unit(0.5, "cm"), 
                                                      type = "closed")))

obj@meta.data <- obj@meta.data %>%
  mutate(sfig3 = 1)
DimPlot_scCustom(obj,
                 label = F,
                 group.by = "sfig3",
                 colors_use = col) + 
  theme(
    legend.position = "none",
    plot.title = element_blank()
  ) +
  theme
ggsave(paste0(plots, "Fig_S3L.png"),
       units = "in", dpi = 600,
       height = 3, width = 3)

#Generate markers file ----
obj <- FindClusters(obj, 
                    resolution = 0.8) #resolution comes from chooseR

markers <- FindAllMarkers(obj, only.pos = F)

write.csv(markers,
          file = "tabular_data/mono-macs_markers.csv")

#Canonical cell markers FeaturePlots based on markers file  ----
theme <- theme(axis.text = element_blank(),
               axis.ticks = element_blank(),
               axis.line = element_line(arrow = arrow(angle = 15, 
                                                      length = unit(0.5, "cm"), 
                                                      type = "closed")))
fplot <- function(gene, filename){
  p <- FeaturePlot_scCustom(obj,
                            features = gene,
                            colors_use = viridis_light_high,
                            na_cutoff = 1) + 
    theme(
      plot.title = element_text(face = "bold.italic")
    ) +
    theme
  
  ggsave(p,
         filename = paste0(plots, filename, ".png"),
         units = "in", dpi = 600,
         height = 4, width = 4.5)
}

fplot("Mrc1", "Fig_S7D")
fplot("Cd163", "Fig_S7E")
fplot("Ms4a7", "Fig_S7F")
fplot("Itgax", "Fig_S7G")
fplot("Cd83", "Fig_S7H")
fplot("Cd86", "Fig_S7I")
fplot("Ccr2", "Fig_S7J")
fplot("Ly6c2", "Fig_S7K")
fplot("Arg1", "Fig_S7L")
fplot("Igkc", "Fig_S7M")
fplot("Vpreb3", "Fig_S7N")
fplot("Hba-a1", "Fig_S7O")

#Annotating and saving annotated object ----
levels(obj) <- as.character(0:17)

ids <- c("Monocytes", "Macrophages", "Macrophages", "Macrophages", "Monocytes",
         "Monocytes", "Macrophages", "Dendritic cells", "Macrophages", "Macrophages",
         "Monocytes", "Dendritic cells", "Dendritic cells", "Monocytes", "Monocytes",
         "Plasma cells", "Monocytes", "Erythrocytes")

levels(obj)
levels <- as.character(c(0:17))
names(ids) <- levels

obj <- RenameIdents(obj, ids)
levels(obj)

DimPlot_scCustom(obj, label = F,
                 pt.size = 2,
                 colors_use = ColorBlind_Pal()) +
  theme
ggsave(paste0(plots, "Fig_S7P.png"),
       units = "in", dpi = 600,
       height = 5, width = 6.5)

#Save object ----
saveRDS(obj, 
        file = "data_objects/06_annotated_mono-macs.RDS")