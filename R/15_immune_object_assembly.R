library(Seurat)
library(tidyverse)

microglia <- readRDS("data_objects/08_microglia_annotated.RDS")
macrophages <- readRDS("data_objects/08_macrophages_annotated.RDS")
granulocytes <- readRDS("data_objects/08_granulocytes_annotated.RDS")
mono <- readRDS("data_objects/08_monocytes_annotated.RDS")

microglia <- subset(microglia,
                    subset = artifact_status == "Authentic")
macrophages <- subset(macrophages,
                      subset = artifact_status == "Authentic")
granulocytes <- subset(granulocytes,
                       subset = artifact_status == "Authentic")
mono <- subset(mono,
               subset = artifact_status == "Authentic")

b <- readRDS("data_objects/08_blymph.RDS")
#b <- readRDS("04_data_objects/05_pre-artifact_removal/blymph.RDS")
dc <- readRDS("data_objects/08_dendritic_cells.RDS")
#dc <- readRDS("04_data_objects/05_pre-artifact_removal/dendritic_cells.RDS")
nk <- readRDS("data_objects/08_natural_killer.RDS")
#nk <- readRDS("04_data_objects/05_pre-artifact_removal/natural_killer.RDS")
t <- readRDS("data_objects/08_tcells.RDS")
#t <- readRDS("04_data_objects/05_pre-artifact_removal/tcells.RDS")

obj <- merge(microglia, list(macrophages, granulocytes, 
                             b, dc, mono, nk, t))

message("~~~~~~~~~~~~~~~~~~~~Split object and SCT each~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
obj_list <- SplitObject(obj, split.by = "Batch")
obj_list <- lapply(obj_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                   vars.to.regress = c("percent_mito", "nCount_RNA"),
                   conserve.memory = T)
message("~~~~~~~~~~~~~~~~~~~~Find variable features and prep integration~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
features <- SelectIntegrationFeatures(obj_list, nfeatures = 3000)
obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)
message("~~~~~~~~~~~~~~~~~~~~PCA on each object~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
obj_list <- lapply(obj_list, FUN = RunPCA, features = features, npcs = 100)
message("~~~~~~~~~~~~~~~~~~~~Find integration anchors~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
anchors1 <- FindIntegrationAnchors(obj_list, normalization.method = "SCT",
                                   anchor.features = features, dims = 1:30,
                                   reduction = "rpca", k.anchor = 5)
message("~~~~~~~~~~~~~~~~~~~~Integrate data~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
intg1 <- IntegrateData(anchorset = anchors1, normalization.method = "SCT", 
                       dims = 1:30)
message("~~~~~~~~~~~~~~~~~~~~Perform dimensional reduction and initial clustering~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
intg1 <- RunPCA(intg1, verbose = F, npcs = 100)

intg1 <- RunUMAP(intg1, reduction = "pca", dims = 1:30, verbose = F) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)

Idents(intg1) <- "cell_type"

saveRDS(intg1, file = "data_objects/09_immune_object.RDS")

message("~~~~~~~~~~~~~~~~~~~~Re-integrate microglia since artifacts were removed~~~~~~~~~~~~~~~~~~~~")
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

saveRDS(microglia, file = "data_objects/08_microglia_clean.RDS")

message("~~~~~~~~~~~~~~~~~~~~Re-integrate macrophages since artifacts were removed~~~~~~~~~~~~~~~~~~~~")
mac_list <- SplitObject(macrophages, split.by = "Batch")
mac_list <- lapply(mac_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                   vars.to.regress = c("percent_mito", "nCount_RNA"))
mac_features <- SelectIntegrationFeatures(mac_list, nfeatures = 3000)
mac_list <- PrepSCTIntegration(mac_list, anchor.features = mac_features)
mac_list <- lapply(mac_list, FUN = RunPCA, features = mac_features, npcs = 100)
mac_anchors <- FindIntegrationAnchors(mac_list, normalization.method = "SCT",
                                      anchor.features = mac_features, dims = 1:30,
                                      reduction = "rpca", k.anchor = 5)
mac <- IntegrateData(anchorset = mac_anchors, normalization.method = "SCT", 
                     dims = 1:30)

mac <- RunPCA(mac, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(mac, file = "data_objects/08_macrophages_clean.RDS")

message("~~~~~~~~~~~~~~~~~~~~Re-integrate granulocytes since artifacts were removed~~~~~~~~~~~~~~~~~~~~")
gran_list <- SplitObject(granulocytes, split.by = "Batch")
gran_list <- lapply(gran_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                    vars.to.regress = c("percent_mito", "nCount_RNA"))
gran_features <- SelectIntegrationFeatures(gran_list, nfeatures = 3000)
gran_list <- PrepSCTIntegration(gran_list, anchor.features = gran_features)
gran_list <- lapply(gran_list, FUN = RunPCA, features = gran_features, npcs = 100)
gran_anchors <- FindIntegrationAnchors(gran_list, normalization.method = "SCT",
                                       anchor.features = gran_features, dims = 1:30,
                                       reduction = "rpca", k.anchor = 5)
gran <- IntegrateData(anchorset = gran_anchors, normalization.method = "SCT", 
                      dims = 1:30)

gran <- RunPCA(gran, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(gran, file = "data_objects/08_granulocytes_clean.RDS")

message("~~~~~~~~~~~~~~~~~~~~Re-integrate monocytes since artifacts were removed~~~~~~~~~~~~~~~~~~~~")
mono_list <- SplitObject(mono, split.by = "Batch")
mono_list <- lapply(mono_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                    vars.to.regress = c("percent_mito", "nCount_RNA"))
mono_features <- SelectIntegrationFeatures(mono_list, nfeatures = 3000)
mono_list <- PrepSCTIntegration(mono_list, anchor.features = mono_features)
mono_list <- lapply(mono_list, FUN = RunPCA, features = mono_features, npcs = 100)
mono_anchors <- FindIntegrationAnchors(mono_list, normalization.method = "SCT",
                                       anchor.features = mono_features, dims = 1:30,
                                       reduction = "rpca", k.anchor = 5)
monocytes <- IntegrateData(anchorset = mono_anchors, normalization.method = "SCT", 
                           dims = 1:30)

monocytes <- RunPCA(monocytes, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(monocytes, file = "data_objects/08_monocytes_clean.RDS")
