library(Seurat)
library(scCustomize)
library(patchwork)
library(dittoSeq)
library(tidyverse)
library(paletteer)

obj <- readRDS("data_objects/05_clean_bcells_sub.RDS")

dir.create("plots/03_fine_cluster_id")
plots <- "plots/03_fine_cluster_id/"

#Make unicolored plot for Fig S3 ----
pal <- dittoColors()
col <- pal[19]

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
                 colors_use = col,
                 pt.size = 2) + 
  theme(
    legend.position = "none",
    plot.title = element_blank()
  ) +
  theme
ggsave(paste0(plots, "Fig_S3I.png"),
       units = "in", dpi = 600,
       height = 3, width = 3)

#Generate markers file ----
obj <- FindClusters(obj, 
                    resolution = 0.6) #resolution comes from chooseR

markers <- FindAllMarkers(obj, only.pos = F)

write.csv(markers,
          file = "tabular_data/bcells_markers.csv")

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

fplot("Cd79a", "Fig_S4D")
fplot("Ms4a1", "Fig_S4E")
fplot("Ighm", "Fig_S4F")
fplot("Cd3e", "Fig_S4G")
fplot("Nkg7", "Fig_S4H")