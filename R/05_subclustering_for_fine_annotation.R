library(glmGamPoi)
library(Seurat)
library(tidyverse)
library(ggplot2)

data <- "data_objects/05_"

obj <- readRDS("data_objects/04_coarse_annotated.RDS")
unique(obj@active.ident)
#Microglia----------------------------------------------------------------------
print("Starting microglia")

microglia <- subset(obj, subset = cell_type == "Microglia")

microglia_list <- SplitObject(microglia, split.by = "Batch")
microglia_list <- lapply(microglia_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                         vars.to.regress = c("percent_mito", "nCount_RNA"))
microglia_features <- SelectIntegrationFeatures(microglia_list, nfeatures = 3000)
microglia_list <- PrepSCTIntegration(microglia_list, anchor.features = microglia_features)
microglia_list <- lapply(microglia_list, FUN = RunPCA, features = microglia_features, npcs = 100)
microglia_anchors <- FindIntegrationAnchors(microglia_list, normalization.method = "SCT",
                                            anchor.features = microglia_features, dims = 1:30,
                                            reduction = "rpca", k.anchor = 5)
microglia <- IntegrateData(anchorset = microglia_anchors, normalization.method = "SCT", 
                           dims = 1:30)

microglia <- RunPCA(microglia, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(microglia, file = paste0(data, "clean_microglia_sub.RDS"))

print("Done with microglia")
print("Starting monocytes & macrophages")
#Monocytes, macrophages---------------------------------------------------------
monomacs <- subset(obj, cell_type == "Monocytes/Macrophages")

monomacs_list <- SplitObject(monomacs, split.by = "Batch")
monomacs_list <- lapply(monomacs_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                        vars.to.regress = c("percent_mito", "nCount_RNA"))
monomacs_features <- SelectIntegrationFeatures(monomacs_list, nfeatures = 3000)
monomacs_list <- PrepSCTIntegration(monomacs_list, anchor.features = monomacs_features)
monomacs_list <- lapply(monomacs_list, FUN = RunPCA, features = monomacs_features, npcs = 100)
monomacs_anchors <- FindIntegrationAnchors(monomacs_list, normalization.method = "SCT",
                                           anchor.features = monomacs_features, dims = 1:30,
                                           reduction = "rpca", k.anchor = 5)
monomacs <- IntegrateData(anchorset = monomacs_anchors, normalization.method = "SCT", 
                          dims = 1:30)

monomacs <- RunPCA(monomacs, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(monomacs, file = paste0(data, "clean_mono-macs_sub.RDS"))

print("Done with monocytes & macrophages")
print("Starting T-cells and NK cells")
#T-cells, NK cells--------------------------------------------------------------
tnk <- subset(obj, cell_type == "T/NK cells")

tnk <- SCTransform(tnk, vst.flavor = "v2",
                   vars.to.regress = c("nCount_RNA", "percent_mito"),
                   conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(tnk, file = paste0(data, "clean_tnk_sub.RDS"))

print("Done with T-cells and NK cells")
print("Starting B-cells")
#B-cells------------------------------------------------------------------------
b <- subset(obj, cell_type == "B-cells")

b <- SCTransform(b, vst.flavor = "v2",
                 vars.to.regress = c("nCount_RNA", "percent_mito"),
                 conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

b <- RunPCA(b, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(b, file = paste0(data, "clean_bcells_sub.RDS"))

print("Done with B-cells")
print("Starting neutrophils")
#Neutrophils--------------------------------------------------------------------
neut <- subset(obj, cell_type == "Neutrophils")

neut_list <- SplitObject(neut, split.by = "Batch")
neut_list <- lapply(neut_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                    vars.to.regress = c("percent_mito", "nCount_RNA"))
neut_features <- SelectIntegrationFeatures(neut_list, nfeatures = 3000)
neut_list <- PrepSCTIntegration(neut_list, anchor.features = neut_features)
neut_list <- lapply(neut_list, FUN = RunPCA, features = neut_features, npcs = 100)
neut_anchors <- FindIntegrationAnchors(neut_list, normalization.method = "SCT",
                                       anchor.features = neut_features, dims = 1:30,
                                       reduction = "rpca", k.anchor = 5)
neut <- IntegrateData(anchorset = neut_anchors, normalization.method = "SCT", 
                      dims = 1:30)

neut <- RunPCA(neut, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(neut, file = paste0(data, "clean_neutrophils_sub.RDS"))

print("Done with neutrophils")
print("Starting non-immune cells")
#Non-immune cells---------------------------------------------------------------
cd45neg <- subset(obj, cell_type != "Microglia" & 
                    cell_type != "Neutrophils" &
                    cell_type != "B-cells" &
                    cell_type != "T/NK cells" & 
                    cell_type != "Monocytes/Macrophages")

cd45neg <- SCTransform(cd45neg, vst.flavor = "v2",
                       vars.to.regress = c("nCount_RNA", "percent_mito"),
                       conserve.memory = T) %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

cd45neg <- RunPCA(cd45neg, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(cd45neg, file = paste0(data, "clean_cd45neg_sub.RDS"))

message("Done!")
