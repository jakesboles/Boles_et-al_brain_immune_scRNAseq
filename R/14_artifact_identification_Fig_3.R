#This script compares DMSO- and Inhibitors-treated samples, identifying "artifactual"
#populations in our dataset. Plots made here are found in Figure 3 of our pre-print. 
#For rigor's sake, all cell types are shown for certain plots, but only granulocytes, 
#macrophages, monocytes, and microglia were shown in the pre-print because these cell 
#types showed a significant difference in the cross entropy test from 13-12. 

library(Seurat)
library(scCustomize)
library(tidyverse)
library(ggplot2)
library(dittoSeq)

dir.create("plots/Fig_3")
plots <- "plots/Fig_3/"

#Aesthetic fixtures
theme <- theme(axis.text = element_blank(),
               axis.ticks = element_blank(),
               axis.line = element_line(arrow = arrow(angle = 15, 
                                                      length = unit(0.5, "cm"), type = "closed")))
pal <- DiscretePalette_scCustomize(num_colors = 20, palette = "varibow")
PalettePlot(pal)

pal_num <- c(2, 1, 5, 14, 9, 11, 12, 13, 19)

#Load objects and add module scores as needed
cells <- c("blymph", "cd45neg", "dendritic_cells", "granulocytes", 
           "macrophages", "microglia", "monocytes", "natural_killer", "tcells")
titles <- c("B-lymphocytes", "CD45-", "Dendritic cells", "Granulocytes",
            "Macrophages", "Microglia", "Monocytes", "Natural killer cells", "T-cells")

objs <- list()
for (i in seq_along(cells)){
  objs[[i]] <- readRDS(paste0("data_objects/08_", cells[i], ".RDS"))
}

names(objs) <- cells

marsh_allCNS_artifact <- list(c("Fos", "Junb", "Zfp36", "Jun", "Hspa1a", "Socs3",
                                "Rgs1", "Egr1", "Btg2", "Fosb", "Ier5", "Atf3", "Dusp1",
                                "Folr1", "Serpine1"))

marsh_myeloid_artifact <- list(c("Rgs1", "Nfkbiz", "Klf2", "Junb", "Dusp1", "Ccl3", 
                                 "Hspa1a", "Hsp90aa1", "Fos", "Hspa1b", "Jun", "Jund", 
                                 "Nfkbid", "Gem", "Ccl4", "Ier5", "Txnip", "Zfp36",
                                 "Egr1", "Atf3", "Rhob"))

for (i in seq_along(objs)){
  objs[[i]] <- AddModuleScore(objs[[i]], features = marsh_allCNS_artifact,
                              name = "marsh_allCNS_artifact")
  objs[[i]] <- AddModuleScore(objs[[i]], marsh_myeloid_artifact, 
                              name = "marsh_myeloid_artifact")
  
  objs[[i]]@meta.data <- objs[[i]]@meta.data %>%
    mutate(f3 = 1)
}

for (i in seq_along(objs)){
  #Show ex vivo treatment groups:
  p1 <- DimPlot_All_Samples(objs[[i]], 
                            meta_data_column = "Dissociation",
                            colors_use = pal[[pal_num[i]]]) +
    NoLegend() + 
    theme +
    ggtitle(titles[i])
  
  for (j in 1:2){
    p1[[j]] <- p1[[j]] + theme
  }
  
  ggsave(p,
         filename = paste0(plots, "by_group_", names(objs)[i], ".png"),
         units = "in", dpi = 600,
         height = 3, width = 6)
  
  #Artifact 1 from Marsh et al 2020 Nat Neurosci
  p2 <- FeaturePlot_scCustom(objs[[i]],
                            features = "marsh_allCNS_artifact1",
                            colors_use = viridis_inferno_light_high,
                            na_cutoff = NA) + 
    theme +
    ggtitle("Artifact module 1")
  
  ggsave(p2,
         filename = paste0(plots, "marsh_allCNS_artifact_", names(objs)[i], ".png"),
         units = "in", dpi = 600,
         height = 3, width = 3)
  
  #Artifact 2 from Marsh et al 2020 Nat Neurosci
  p <- FeaturePlot_scCustom(objs[[i]],
                            features = "marsh_myeloid_artifact1",
                            colors_use = viridis_inferno_light_high,
                            na_cutoff = NA) + 
    theme +
    ggtitle("Artifact module 2")
  
  ggsave(p,
         filename = paste0(plots, "marsh_myeloid_artifact_", names(objs)[i], ".png"),
         units = "in", dpi = 600,
         height = 3, width = 3)
}

microglia <- objs[[6]]
macrophages <- objs[[5]]
granulocytes <- objs[[4]]
monocytes <- objs[[7]]

#Annotate microglia to mark artifact clusters ----
microglia <- FindClusters(microglia,
                          resolution = 0.5) #resolution was selected by manual iteration

DimPlot_scCustom(microglia,
                 label = F) + 
  theme
ggsave(paste0(plots, "dimplot_microglia.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)

ids <- c(rep("Authentic", 4), rep("Artificial", 3), rep("Authentic", 11))
levels(microglia)
levels <- as.character(c(0:17))
names(ids) <- levels
microglia <- RenameIdents(microglia, ids)
levels(microglia)

DimPlot_scCustom(microglia, 
                 label = F,
                 colors_use = c("midnightblue", "deeppink")) + 
  theme
ggsave(paste0(plots, "dimplot_microglia_artifacts.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)

#Annotate macrophages to mark artifact clusters ----
macrophages <- FindClusters(macrophages,
                            resolution = 0.8) #resolution was selected by manual iteration

DimPlot_scCustom(macrophages,
                 label = F) + 
  theme
ggsave(paste0(plots, "dimplot_macrophages.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)
ids <- c(rep("Authentic", 6), rep("Artificial", 2), "Authentic", "Artificial",
         rep("Authentic", 6))
levels(macrophages)

levels <- as.character(c(0:15))
names(ids) <- levels
macrophages <- RenameIdents(macrophages, ids)
levels(macrophages)

DimPlot_scCustom(macrophages,
                 label = F,
                 colors_use = c("midnightblue", "deeppink")) + 
  theme
ggsave(paste0(plots, "dimplot_macrophages_artifacts.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)
#Annotate granulocytes to mark artifact clusters ----
granulocytes <- FindClusters(granulocytes,
                             resolution = 0.6) #resolution was selected by manual iteration

DimPlot_scCustom(granulocytes,
                 label = F,
                 pt.size = 3) + 
  theme
ggsave(paste0(plots, "dimplot_granulocytes.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)
DimPlot_scCustom(granulocytes,
                 group.by = "cell_type",
                 colors_use = c(pal[3], pal[10], pal[14]))

ids <- c("Authentic", "Artificial", rep("Authentic", 9))
levels <- as.character(c(0:10))
names(ids) <- levels

granulocytes <- RenameIdents(granulocytes, ids)
levels(granulocytes)

DimPlot_scCustom(granulocytes,
                 label = F,
                 colors_use = c("midnightblue", "deeppink"),
                 pt.size = 3) + 
  theme
ggsave(paste0(plots, "dimplot_granulocytes_artifacts.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)
#Annotate monocytes to mark artifact clusters ----
monocytes <- FindClusters(monocytes,
                          resolution = 1) #resolution was selected by manual iteration

DimPlot_scCustom(monocytes,
                 label = F,
                 pt.size = 2) + 
  theme
ggsave(paste0(plots, "dimplot_monocytes.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)

ids <- c(rep("Authentic", 13), "Artificial", rep("Authentic", 3))
levels <- as.character(c(0:16))
names(ids) <- levels

monocytes <- RenameIdents(monocytes, ids)
levels(monocytes)

DimPlot_scCustom(monocytes,
                 label = F,
                 colors_use = c("midnightblue", "deeppink"),
                 pt.size = 2) + 
  theme
ggsave(paste0(plots, "dimplot_monocytes_artifacts.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)

#Violin plots with module scores ----
pal <- DiscretePalette_scCustomize(num_colors = 20,
                                   palette = "polychrome")

#Macrophages
VlnPlot_scCustom(macrophages,
                 group.by = "integrated_snn_res.0.8",
                 features = "marsh_allCNS_artifact1",
                 #split.by = "Dissociation",
                 colors_use = pal[1:18],
                 raster = F,
                 pt.size = 0.001) + 
  ggtitle("Artifact module 1") +
  NoLegend() + 
  theme(
    axis.text.x = element_text(angle = 0)
  )
ggsave(paste0(plots, "marsh_allCNS_artifact_macrophages_vln.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

#Microglia
VlnPlot_scCustom(microglia,
                 group.by = "integrated_snn_res.0.5",
                 features = "marsh_allCNS_artifact1",
                 #split.by = "Dissociation",
                 raster = F,
                 pt.size = 0.001,
                 colors_use = pal[1:18]) + 
  ggtitle("Artifact module 1") +
  NoLegend() +
  theme(
    axis.text.x = element_text(angle = 0)
  )
ggsave(paste0(plots, "marsh_allCNS_artifact_microglia_vln.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

#Granulocytes
VlnPlot_scCustom(granulocytes,
                 group.by = "integrated_snn_res.0.6",
                 features = "marsh_allCNS_artifact1",
                 #split.by = "Dissociation",
                 colors_use = pal[1:18],
                 raster = F,
                 pt.size = 0.001) + 
  ggtitle("Artifact module 1") +
  NoLegend() + 
  theme(
    axis.text.x = element_text(angle = 0)
  )
ggsave(paste0(plots, "marsh_allCNS_artifact_granulocytes_vln.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

#Monocytes
VlnPlot_scCustom(monocytes,
                 group.by = "integrated_snn_res.1",
                 features = "marsh_allCNS_artifact1",
                 colors_use = pal[1:17],
                 raster = F,
                 pt.size = 0.001) + 
  ggtitle("Artifact module 1") + 
  NoLegend() + 
  theme(axis.text.x = element_text(angle = 0))
ggsave(paste0(plots, "marsh_allCNS_artifact_monocytes_vln.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

#Sample composition bar charts to show which clusters come from DMSO-treated cells ----
#Monocytes
monocytes$integrated_snn_res.1 <- factor(monocytes$integrated_snn_res.1,
                                         levels = as.character(c(0:16)))
dittoBarPlot(monocytes,
             var = "Dissociation",
             group.by = "integrated_snn_res.1",
             main = "% of clusters from each group",
             ylab = "% of cluster",
             xlab = "Identity",
             x.labels.rotate = F,
             retain.factor.levels = T,
             color.panel = c("deeppink", "midnightblue")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    axis.ticks.x = element_blank()
  )

ggsave(paste0(plots, "sample_comp_monocytes.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

#Microglia
microglia$integrated_snn_res.0.5 <- factor(microglia$integrated_snn_res.0.5,
                                           levels = as.character(c(0:17)))

dittoBarPlot(microglia, 
             var = "Dissociation",
             group.by = "integrated_snn_res.0.5",
             main = "% of clusters from each group",
             ylab = "% of cluster",
             xlab = "Identity",
             x.labels.rotate = F,
             retain.factor.levels = T,
             color.panel = c("deeppink", "midnightblue")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    axis.ticks.x = element_blank()
  )

ggsave(paste0(plots, "sample_comp_microglia.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

#Macrophages
dittoBarPlot(macrophages, 
             var = "Dissociation",
             group.by = "integrated_snn_res.0.8",
             main = "% of clusters from each group",
             ylab = "% of cluster",
             xlab = "Identity",
             x.labels.rotate = F,
             retain.factor.levels = T,
             color.panel = c("deeppink", "midnightblue")) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    axis.ticks.x = element_blank()
  )
ggsave(paste0(plots, "sample_comp_macrophages.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

#Granulocytes
dittoBarPlot(granulocytes, 
             var = "Dissociation",
             group.by = "integrated_snn_res.0.6",
             main = "% of clusters from each group",
             ylab = "% of cluster",
             xlab = "Identity",
             x.labels.rotate = F,
             retain.factor.levels = T,
             color.panel = c("deeppink", "midnightblue")) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    axis.ticks.x = element_blank()
  )
ggsave(paste0(plots, "sample_comp_granulocytes.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

#Mark and save objects for artifact removal ----
microglia[["artifact_status"]] <- Idents(object = microglia)
granulocytes[["artifact_status"]] <- Idents(object = granulocytes)
macrophages[["artifact_status"]] <- Idents(object = macrophages)
monocytes[["artifact_status"]] <- Idents(object = monocytes)

saveRDS(microglia,
        file = "data_objects/08_microglia_annotated.RDS")
saveRDS(macrophages,
        file = "data_objects/08_macrophages_annotated.RDS")
saveRDS(granulocytes,
        file = "data_objects/08_granulocytes_annotated.RDS")
saveRDS(monocytes,
        file = "data_objects/08_monocytes_annotated.RDS")