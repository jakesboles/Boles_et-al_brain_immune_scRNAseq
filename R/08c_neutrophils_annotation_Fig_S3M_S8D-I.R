library(Seurat)
library(scCustomize)
library(patchwork)
library(dittoSeq)
library(tidyverse)
library(paletteer)

obj <- readRDS("data_objects/05_clean_neutrophils_sub.RDS")

dir.create("plots/03_fine_cluster_id")
plots <- "plots/03_fine_cluster_id/"

#Need to run PrepSCTFindMarkers because the data were integrated in 05
obj <- PrepSCTFindMarkers(obj)

#Make unicolored plot for Fig S3 ----
pal <- dittoColors()
col <- pal[23]

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
ggsave(paste0(plots, "Fig_S3M.png"),
       units = "in", dpi = 600,
       height = 3, width = 3)

#Generate markers file ----
obj <- FindClusters(obj,
                    resolution = 0.6)

markers <- FindAllMarkers(obj, only.pos = F)

write.csv(markers,
          file = "tabular_data/neutrophils_markers.csv")

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

fplot("Mmp8", "Fig_S8D")
fplot("Ngp", "Fig_S8E")
fplot("S100a8", "Fig_S8F")
fplot("S100a9", "Fig_S8G")
fplot("P2ry12", "Fig_S8H")

#Annotate object ----
ids <- c(rep("Neutrophils", 8), "Doublets", "Neutrophils")
names(ids) <- as.character(c(0:9))

obj <- RenameIdents(obj, ids)
DimPlot_scCustom(obj, label = F,
                 pt.size = 1) +
  theme
ggsave(paste0(plots, "Fig_S8I.png"),
       units = "in", dpi = 600,
       height = 5, width = 6.5)

#Save object ----
saveRDS(obj, 
        file = "data_objects/06_annotated_neutrophils.RDS")