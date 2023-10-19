library(Seurat)
library(scCustomize)
library(patchwork)
library(dittoSeq)
library(tidyverse)
library(paletteer)

obj <- readRDS("data_objects/05_clean_microglia_sub.RDS")

dir.create("plots/03_fine_cluster_id")
plots <- "plots/03_fine_cluster_id/"

#Need to run PrepSCTFindMarkers because the data were integrated in 05
obj <- PrepSCTFindMarkers(obj)

#Make unicolored plot for Fig S3 ----
pal <- dittoColors()
col <- pal[21]

theme <- theme(axis.text = element_blank(),
               axis.ticks = element_blank(),
               axis.line = element_line(arrow = arrow(angle = 15, 
                                                      length = unit(0.5, "cm"), 
                                                      type = "closed")))

DefaultAssay(obj) <- "SCT"

obj@meta.data <- obj@meta.data %>%
  mutate(sfig3 = 1)
DimPlot_scCustom(obj,
                 label = F,
                 group.by = "sfig3",
                 colors_use = col,
                 pt.size = 0.5) + 
  theme(
    legend.position = "none",
    plot.title = element_blank()
  ) +
  theme
ggsave(paste0(plots, "Fig_S3K.png"),
       units = "in", dpi = 600,
       height = 3, width = 3)

#Generate markers file ----
obj <- FindClusters(obj,
                    resolution = 0.4)

markers <- FindAllMarkers(obj, only.pos = F)

write.csv(markers,
          file = "tabular_data/microglia_markers.csv")

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

fplot("Aif1", "Fig_S6F")
fplot("C1qc", "Fig_S6G")
fplot("Hexb", "Fig_S6H")
fplot("P2ry12", "Fig_S6I")
fplot("Tmem119", "Fig_S6J")

#All clusters are microglia, so no additional annotation needed

#Save object ----
saveRDS(obj, 
        file = "data_objects/06_annotated_microglia.RDS")
