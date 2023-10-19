library(Seurat)
library(scCustomize)
library(patchwork)
library(dittoSeq)
library(tidyverse)
library(paletteer)

obj <- readRDS("data_objects/05_clean_cd45neg_sub.RDS")

dir.create("plots/03_fine_cluster_id")
plots <- "plots/03_fine_cluster_id/"

#Make unicolored plot for Fig S3 ----
pal <- dittoColors()
col <- pal[20]

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
ggsave(paste0(plots, "Fig_S3J.png"),
       units = "in", dpi = 600,
       height = 3, width = 3)

#Generate markers file ----
obj <- FindClusters(obj,
                    resolution = 0.6)

markers <- FindAllMarkers(obj, only.pos = F)

write.csv(markers,
          file = "tabular_data/nonimmune_markers.csv")

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

fplot("Acta2", "Fig_S5D")
fplot("Aqp4", "Fig_SFE")
fplot("Cldn5", "Fig_S5F")
fplot("Col1a2", "Fig_S5G")
fplot("Dcn", "Fig_S5H")
fplot("Folr1", "Fig_S5I")
fplot("Gja1", "Fig_S5J")
fplot("Olig1", "Fig_S5K")
fplot("S100b", "Fig_S5L")
fplot("Slc38a2", "Fig_S5M")

#Annotating and saving annotated object ----
levels(obj) <- as.character(0:16)

ids <- c("Choroid plexus cells", "Endothelial cells", "Astrocytes", "Endothelial cells",
         "Endothelial cells", "Pericytes", "Choroid plexus cells", "Astrocytes",
         "Leptomeningeal cells", "Pericytes", "Leptomeningeal cells", "Astrocytes",
         "Leptomeningeal cells", "Choroid plexus cells", "Oligodendrocytes",
         "Oligodendrocytes", "Astrocytes")

levels <- as.character(c(0:16))
names(ids) <- levels

obj <- RenameIdents(obj, ids)

DimPlot_scCustom(obj, label = F,
                 pt.size = 1,
                 colors_use = Dark2_Pal()) +
  theme
ggsave(paste0(plots, "Fig_S5N.png"),
       units = "in", dpi = 600,
       height = 5, width = 6.5)

#Save object ----
saveRDS(obj, 
        file = "data_objects/06_annotated_cd45neg.RDS")
