library(Seurat)
library(scCustomize)
library(patchwork)
library(dittoSeq)
library(tidyverse)
library(paletteer)

obj <- readRDS("data_objects/05_clean_tnk_sub.RDS")

dir.create("plots/03_fine_cluster_id")
plots <- "plots/03_fine_cluster_id/"

#Make unicolored plot for Fig S3 ----
pal <- dittoColors()
col <- pal[24]

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
ggsave(paste0(plots, "Fig_S3N.png"),
       units = "in", dpi = 600,
       height = 3, width = 3)

#Generate markers file ----
obj <- FindClusters(obj,
                    resolution = 0.6)

markers <- FindAllMarkers(obj, only.pos = F)

write.csv(markers,
          file = "tabular_data/t_nk_markers.csv")

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

fplot("Cd3e", "Fig_S9D")
fplot("Cd8b1", "Fig_S9E")
fplot("Gzma", "Fig_S9F")
fplot("Klrb1c", "Fig_S9G")
fplot("Trdc", "Fig_S9H")
fplot("Trgc1", "Fig_S9I")
fplot("Marco", "Fig_S9J")
fplot("Hba-a1", "Fig_S9K")
fplot("Cebpa", "Fig_S9L")
fplot("Mcpt4", "Fig_S9M")
fplot("Plac8", "Fig_S9N")
fplot("Calca", "Fig_S9O")

#Annotating and saving annotated object ----
levels(obj) <- as.character(0:11)

ids <- c("αβ T-cells", "αβ T-cells", "αβ T-cells", "Natural killer cells",
         "Natural killer cells", "Stem cells", "γδ T-cells", "Basophils", "Macrophages",
         "Stem cells", "Mast cells", "Erythrocytes")

levels(obj)
levels <- as.character(c(0:11))
names(ids) <- levels

obj <- RenameIdents(obj, ids)
levels(obj)

DimPlot_scCustom(obj, label = F,
                 pt.size = 2,
                 colors_use = ColorBlind_Pal()) +
  theme
ggsave(paste0(plots, "Fig_S9P.png"),
       units = "in", dpi = 600,
       height = 5, width = 6.5)

#Save object ----
saveRDS(obj, 
        file = "data_objects/06_annotated_tnk.RDS")