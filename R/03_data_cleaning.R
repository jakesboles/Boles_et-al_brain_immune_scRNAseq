library(Seurat)
library(scCustomize)
library(tidyverse)
library(ggplot2)
library(glmGamPoi)

#In the last script, you may notice lots of cells that are, for example, P2ry12+
#but are amidst a cluster of cells that are clearly T-cells...
#This script goes cluster by cluster and removes what is likely microglial contamination.
#After this cleaning, the object is reassembled for the annotation that went into the paper.

dir.create("plots/03_cleaning_microglial_contamination")
plots <- "plots/03_cleaning_microglial_contamination"

obj <- readRDS("data_objects/02_annotated.RDS")

#The gene module below is derived from https://www.nature.com/articles/s41593-022-01022-8
marsh_microglia_id <- list(c("P2ry12", "Fcrls", "Trem2", "Tmem119", "Cx3cr1", "Hexb",
                             "Tgfbr1","Sparc","P2ry13","Olfml3","Adgrg1","C1qa",
                             "C1qb","C1qc","Csf1r","Fcgr3","Ly86", "Laptm5"))
obj <- AddModuleScore(obj, marsh_microglia_id, 
                      name = "marsh_microglia_id")

#NK/T-cell ----
tnk <- subset(obj, idents = "NK/T-cell")

tnk <- SCTransform(tnk, vst.flavor = "v2",
                   vars.to.regress = c("nCount_RNA", "percent_mito"),
                   conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)

tnk <- FindClusters(tnk, resolution = 0.2)

FeaturePlot(tnk, 
            features = c("Aif1", "Tmem119", "Hexb", "P2ry12"), 
            slot = "counts")

DimPlot_scCustom(tnk, 
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                               num_colors = 8),
                 pt.size = 3)

VlnPlot_scCustom(tnk, 
                 features = c("Aif1", "Tmem119", "Hexb", "P2ry12"), 
                 slot = "counts",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 8))

FeaturePlot(tnk, 
            features = "marsh_microglia_id1", 
            cols = c("navy", "gold"))

VlnPlot_scCustom(tnk, 
                 features = "marsh_microglia_id1",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 8))

tnk.ave <- AverageExpression(tnk, assays = "SCT", return.seurat = T)
DoHeatmap(tnk.ave, features = unlist(TopFeatures(tnk[["pca"]], balanced = TRUE)), size = 3,
          draw.lines = FALSE)


tnk.remove <- WhichCells(tnk, idents = 4)
#What % of cells were removed?
100*(length(tnk.remove)/length(Cells(tnk)))

tnk.filter <- subset(tnk, cells = setdiff(Cells(tnk), tnk.remove))
tnk.filter <- SCTransform(tnk.filter, vst.flavor = "v2",
                          vars.to.regress = c("nCount_RNA", "percent_mito"),
                          conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

#B-cell----
b <- subset(obj, idents = "B-cell")
b <- SCTransform(b, vst.flavor = "v2",
                 vars.to.regress = c("nCount_RNA", "percent_mito"),
                 conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
b <- FindClusters(b, resolution = 0.2)

FeaturePlot(b, 
            features = c("Aif1", "Tmem119", "Hexb", "P2ry12"), 
            slot = "counts")

DimPlot_scCustom(b, 
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 7),
                 pt.size = 3)

VlnPlot_scCustom(b, 
                 features = c("Aif1", "Tmem119", "Hexb", "P2ry12", "Cd19", "Cd74"), 
                 slot = "counts",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 7))

FeaturePlot(b, 
            features = "marsh_microglia_id1", 
            cols = c("navy", "gold"))

VlnPlot_scCustom(b, 
                 features = "marsh_microglia_id1",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 7))

b.ave <- AverageExpression(b, assays = "SCT", return.seurat = T)
DoHeatmap(b.ave, 
          features = unlist(TopFeatures(b[["pca"]], 
                                        balanced = TRUE)), 
          size = 3,
          draw.lines = FALSE)

b.remove <- WhichCells(b, idents = c(3, 6))
#What % of cells were removed?
100*(length(b.remove)/length(Cells(b)))

#Neutrophil----
neut <- subset(obj, idents = "Neutrophil")
neut <- SCTransform(neut, vst.flavor = "v2",
                    vars.to.regress = c("nCount_RNA", "percent_mito"),
                    conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
neut <- FindClusters(neut, resolution = 1.5)

FeaturePlot(neut, 
            features = c("Aif1", "Tmem119", "Hexb", "P2ry12"), 
            slot = "counts")

DimPlot_scCustom(neut, 
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 17))

VlnPlot_scCustom(neut, 
                 features = c("Aif1", "Tmem119", "Hexb", "P2ry12", "S100a6", "Camp", "Ngp", "Itgb2"), 
                 slot = "counts",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 17))

FeaturePlot(neut, 
            features = "marsh_microglia_id1", 
            cols = c("navy", "gold"))

VlnPlot_scCustom(neut, 
                 features = "marsh_microglia_id1",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 17))

neut.ave <- AverageExpression(neut, assays = "SCT", return.seurat = T)
DoHeatmap(neut.ave, features = unlist(TopFeatures(neut[["pca"]], balanced = TRUE)), size = 3,
          draw.lines = FALSE)

neut.remove <- WhichCells(neut, idents = c(11, 10))
#What % of cells were removed?
100*(length(neut.remove)/length(Cells(neut)))

#Astrocyte----
astro <- subset(obj, idents = "Astrocyte")
astro <- SCTransform(astro, vst.flavor = "v2",
                     vars.to.regress = c("nCount_RNA", "percent_mito"),
                     conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
astro <- FindClusters(astro, resolution = 0.4)

FeaturePlot(astro, 
            features = c("Aif1", "Tmem119", "Hexb", "P2ry12"), 
            slot = "counts")

DimPlot_scCustom(astro, 
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 5),
                 pt.size = 3)

VlnPlot_scCustom(astro, 
                 features = c("Aif1", "Tmem119", "Hexb", "P2ry12", "S100b"), 
                 slot = "data",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 5))

FeaturePlot(astro, 
            features = "marsh_microglia_id1", 
            cols = c("navy", "gold"))

VlnPlot_scCustom(astro, 
                 features = "marsh_microglia_id1",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 5))

astro.ave <- AverageExpression(astro, assays = "SCT", return.seurat = T)
DoHeatmap(astro.ave, features = unlist(TopFeatures(astro[["pca"]], balanced = TRUE)), size = 3,
          draw.lines = FALSE)

astro.remove <- WhichCells(astro, idents = c(3,4))
#What % of cells were removed?
100*(length(astro.remove)/length(Cells(astro)))

#CNS/Meninges----
cns <- subset(obj, idents = "Meninges/CNS")
cns <- SCTransform(cns, vst.flavor = "v2",
                   vars.to.regress = c("nCount_RNA", "percent_mito"),
                   conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
cns <- FindClusters(cns, resolution = 0.2)

FeaturePlot(cns, 
            features = c("Aif1", "Tmem119", "Hexb", "P2ry12"), 
            slot = "counts")

DimPlot_scCustom(cns, 
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 11))

VlnPlot_scCustom(cns, 
                 features = c("Aif1", "Tmem119", "Hexb", "P2ry12"), 
                 slot = "counts",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 11))

FeaturePlot(cns, 
            features = "marsh_microglia_id1", 
            cols = c("navy", "gold"))

VlnPlot_scCustom(cns, 
                 features = "marsh_microglia_id1",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 11))

cns.ave <- AverageExpression(cns, assays = "SCT", return.seurat = T)
DoHeatmap(cns.ave, features = unlist(TopFeatures(cns[["pca"]], balanced = TRUE)), size = 3,
          draw.lines = FALSE)

cns.remove <- WhichCells(cns, idents = c(0, 4))
#What % of cells were removed?
100*(length(cns.remove)/length(Cells(cns)))

#Endothelial/Vascular1----
endo1 <- subset(obj, idents = "Endothelial/Vascular1")
endo1 <- SCTransform(endo1, vst.flavor = "v2",
                     vars.to.regress = c("nCount_RNA", "percent_mito"),
                     conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
endo1 <- FindClusters(endo1, resolution = 0.4)

FeaturePlot(endo1, 
            features = c("Aif1", "Tmem119", "Hexb", "P2ry12"), 
            slot = "counts")

DimPlot_scCustom(endo1, 
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 6),
                 pt.size = 3)

VlnPlot_scCustom(endo1, 
                 features = c("Aif1", "Tmem119", "Hexb", "P2ry12", "Cldn5", "Acta2"), 
                 slot = "counts",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 6))

FeaturePlot(endo1, 
            features = "marsh_microglia_id1", 
            cols = c("navy", "gold"))

VlnPlot_scCustom(endo1, 
                 features = "marsh_microglia_id1",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 6))

endo1.ave <- AverageExpression(endo1, assays = "SCT", return.seurat = T)
DoHeatmap(endo1.ave, features = unlist(TopFeatures(endo1[["pca"]], balanced = TRUE)), size = 3,
          draw.lines = FALSE)

endo1.remove <- WhichCells(endo1, idents = 5)
#What % of cells were removed?
100*(length(endo1.remove)/length(Cells(endo1)))

#Endothelial/Vascular2----
endo2 <- subset(obj, idents = "Endothelial/Vascular2")
endo2 <- SCTransform(endo2, vst.flavor = "v2",
                     vars.to.regress = c("nCount_RNA", "percent_mito"),
                     conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
endo2 <- FindClusters(endo2, resolution = 0.8)

FeaturePlot(endo2, 
            features = c("Aif1", "Tmem119", "Hexb", "P2ry12"), 
            slot = "counts")

DimPlot_scCustom(endo2, 
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 6),
                 pt.size = 3)

VlnPlot_scCustom(endo2, 
                 features = c("Aif1", "Tmem119", "Hexb", "P2ry12", "Icam1", "Pglyrp1"), 
                 slot = "counts",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 6))

FeaturePlot(endo2, 
            features = "marsh_microglia_id1", 
            cols = c("navy", "gold"))

VlnPlot_scCustom(endo2, 
                 features = "marsh_microglia_id1",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 6))

endo2.ave <- AverageExpression(endo2, assays = "SCT", return.seurat = T)
DoHeatmap(endo2.ave, features = unlist(TopFeatures(endo2[["pca"]], balanced = TRUE)), size = 3,
          draw.lines = FALSE)

endo2.remove <- WhichCells(endo2, idents = 5)
#What % of cells were removed?
100*(length(endo2.remove)/length(Cells(endo2)))

#Choroid plexus----
cp <- subset(obj, idents = "Choroid plexus")
cp <- SCTransform(cp, vst.flavor = "v2",
                  vars.to.regress = c("nCount_RNA", "percent_mito"),
                  conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
cp <- FindClusters(cp, resolution = 0.2)

FeaturePlot(cp, 
            features = c("Aif1", "Tmem119", "Hexb", "P2ry12"), 
            slot = "counts")

DimPlot_scCustom(cp, 
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 3),
                 pt.size = 3)

VlnPlot_scCustom(cp, 
                 features = c("Aif1", "Tmem119", "Hexb", "P2ry12",
                              "Dcn", "Ttr", "Folr1", "Prlr"), 
                 slot = "counts",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 3))

FeaturePlot(cp, 
            features = "marsh_microglia_id1", 
            cols = c("navy", "gold"))

VlnPlot_scCustom(cp, 
                 features = "marsh_microglia_id1",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 3))

cp.ave <- AverageExpression(cp, assays = "SCT", return.seurat = T)
DoHeatmap(cp.ave, features = unlist(TopFeatures(cp[["pca"]], balanced = TRUE)), size = 3,
          draw.lines = FALSE)

cp.remove <- WhichCells(cp, idents = 2)
#What % of cells were removed?
100*(length(cp.remove)/length(Cells(cp)))

#Monocyte/Macrophage1----
mono1 <- subset(obj, idents = "Monocyte/Macrophage1")
mono1 <- SCTransform(mono1, vst.flavor = "v2",
                     vars.to.regress = c("nCount_RNA", "percent_mito"),
                     conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
mono1 <- FindClusters(mono1, resolution = 0.4)

FeaturePlot(mono1, 
            features = c("Hexb", "P2ry12", "Sall1", "Tmem119"), 
            slot = "counts")

DimPlot_scCustom(mono1, 
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 10),
                 pt.size = 2)

VlnPlot_scCustom(mono1, 
                 features = c("Hexb", "P2ry12", "Sall1", "Tmem119"), 
                 slot = "counts",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 10))

FeaturePlot(mono1, 
            features = "marsh_microglia_id1", 
            cols = c("navy", "gold"))

VlnPlot_scCustom(mono1, 
                 features = "marsh_microglia_id1",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 10))

mono1.ave <- AverageExpression(mono1, assays = "SCT", return.seurat = T)
DoHeatmap(mono1.ave, features = unlist(TopFeatures(mono1[["pca"]], balanced = TRUE)), size = 3,
          draw.lines = FALSE)

mono1.remove <- WhichCells(mono1, idents = 8)
#What % of cells were removed?
100*(length(mono1.remove)/length(Cells(mono1)))

#Monocyte/Macrophage2----
mono2 <- subset(obj, idents = "Monocyte/Macrophage2")
mono2 <- SCTransform(mono2, vst.flavor = "v2",
                     vars.to.regress = c("nCount_RNA", "percent_mito"),
                     conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
mono2 <- FindClusters(mono2, resolution = 0.2)

FeaturePlot(mono2, 
            features = c("Hexb", "P2ry12", "Sall1", "Tmem119"), 
            slot = "counts")

DimPlot_scCustom(mono2, 
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 5))

VlnPlot_scCustom(mono2, 
                 features = c("Hexb", "P2ry12", "Sall1", "Tmem119"), 
                 slot = "counts",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 5))

FeaturePlot(mono2, 
            features = "marsh_microglia_id1", 
            cols = c("navy", "gold"))

VlnPlot_scCustom(mono2, 
                 features = "marsh_microglia_id1",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 5))

mono2.ave <- AverageExpression(mono2, assays = "SCT", return.seurat = T)
DoHeatmap(mono2.ave, features = unlist(TopFeatures(mono2[["pca"]], balanced = TRUE)), size = 3,
          draw.lines = FALSE)

mono2.remove <- WhichCells(mono2, idents = 3)
#What % of cells were removed?
100*(length(mono2.remove)/length(Cells(mono2)))

#Monocyte/Macrophage3----
mono3 <- subset(obj, idents = "Monocyte/Macrophage3")
mono3 <- SCTransform(mono3, vst.flavor = "v2",
                     vars.to.regress = c("nCount_RNA", "percent_mito"),
                     conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
mono3 <- FindClusters(mono3, resolution = 0.4)

FeaturePlot(mono3, 
            features = c("Hexb", "P2ry12", "Sall1", "Tmem119"), 
            slot = "counts")

DimPlot_scCustom(mono3, 
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 6))

VlnPlot_scCustom(mono3, 
                 features = c("Hexb", "P2ry12", "Sall1", "Tmem119"), 
                 slot = "counts",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 6))

FeaturePlot(mono3, 
            features = "marsh_microglia_id1", 
            cols = c("navy", "gold"))

VlnPlot_scCustom(mono3, 
                 features = "marsh_microglia_id1",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 6))

mono3.ave <- AverageExpression(mono3, assays = "SCT", return.seurat = T)
DoHeatmap(mono3.ave, features = unlist(TopFeatures(mono3[["pca"]], balanced = TRUE)), size = 3,
          draw.lines = FALSE)

mono3.remove <- WhichCells(mono3, idents = 5)
#What % of cells were removed?
100*(length(mono3.remove)/length(Cells(mono3)))
#Monocyte/Macrophage4----
mono4 <- subset(obj, idents = "Monocyte/Macrophage4")
mono4 <- SCTransform(mono4, vst.flavor = "v2",
                     vars.to.regress = c("nCount_RNA", "percent_mito"),
                     conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
mono4 <- FindClusters(mono4, resolution = 0.2)

FeaturePlot(mono4, 
            features = c("Hexb", "P2ry12", "Sall1", "Tmem119"), 
            slot = "counts")

DimPlot_scCustom(mono4, 
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 6))

VlnPlot_scCustom(mono4, 
                 features = c("Hexb", "P2ry12", "Sall1", "Tmem119"), 
                 slot = "counts",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 6))

FeaturePlot(mono4, 
            features = "marsh_microglia_id1", 
            cols = c("navy", "gold"))

VlnPlot_scCustom(mono4, 
                 features = "marsh_microglia_id1",
                 colors_use = DiscretePalette_scCustomize(palette = "ditto_seq",
                                                          num_colors = 6))

mono4.ave <- AverageExpression(mono4, assays = "SCT", return.seurat = T)
DoHeatmap(mono4.ave, features = unlist(TopFeatures(mono4[["pca"]], balanced = TRUE)), size = 3,
          draw.lines = FALSE)

mono4.remove <- WhichCells(mono4, idents = 4)
#What % of cells were removed?
100*(length(mono4.remove)/length(Cells(mono4)))

#Now the actual subsetting and re-integrating----
remove <- c(tnk.remove, b.remove, neut.remove, astro.remove, 
            endo1.remove, endo2.remove, cns.remove, cp.remove, 
            mono1.remove, mono2.remove, mono3.remove, mono4.remove)
#What % of cells will be removed from the full object?
100*(length(remove)/length(Cells(obj)))

obj2 <- subset(obj, cells = setdiff(Cells(obj), remove))

obj_list <- SplitObject(obj2, split.by = "Batch")
obj_list <- lapply(obj_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                   vars.to.regress = c("percent_mito", "nCount_RNA"))
features <- SelectIntegrationFeatures(obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)
obj_list <- lapply(obj_list, FUN = RunPCA, features = features, npcs = 100)
anchors1 <- FindIntegrationAnchors(obj_list, normalization.method = "SCT",
                                   anchor.features = features, dims = 1:30,
                                   reduction = "rpca", k.anchor = 5)
intg1 <- IntegrateData(anchorset = anchors1, normalization.method = "SCT", 
                       dims = 1:30)
rm(list = c("obj_list", "features", "anchors1"))
#Performing PCA, dimensional reduction on integrated object----
intg1 <- RunPCA(intg1, verbose = F, npcs = 100)

intg1 <- RunUMAP(intg1, reduction = "pca", dims = 1:30, verbose = F) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

saveRDS(intg1, "data_objects/03_cleaned.RDS")
