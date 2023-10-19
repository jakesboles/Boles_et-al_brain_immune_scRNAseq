library(glmGamPoi)
library(Seurat)
library(tidyverse)

obj <- readRDS("data_objects/07_annotated_object.rds")
levels(obj)

message("~~~~~~~~~~~~~~~~~~~~Starting NK cells~~~~~~~~~~~~~~~~~~~~")
nk <- subset(obj, subset = cell_type == "Natural killer cells")

nk_list <- SplitObject(nk, split.by = "Batch")
nk_list <- lapply(nk_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                  vars.to.regress = c("percent_mito", "nCount_RNA"))
nk_features <- SelectIntegrationFeatures(nk_list, nfeatures = 3000)
nk_list <- PrepSCTIntegration(nk_list, anchor.features = nk_features)
nk_list <- lapply(nk_list, FUN = RunPCA, features = nk_features, npcs = 100)
nk_anchors <- FindIntegrationAnchors(nk_list, normalization.method = "SCT",
                                     anchor.features = nk_features, dims = 1:30,
                                     reduction = "rpca", k.anchor = 5)
nk <- IntegrateData(anchorset = nk_anchors, normalization.method = "SCT", 
                    dims = 1:30,
                    k.weight = 50)

nk <- RunPCA(nk, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(nk, file = "data_objects/08_natural_killer.RDS")
message("~~~~~~~~~~~~~~~~~~~~Starting T-cells~~~~~~~~~~~~~~~~~~~~")
t <- subset(obj, subset = cell_type == "αβ T-cells" | 
              cell_type == "γδ T-cells")

t_list <- SplitObject(t, split.by = "Batch")
t_list <- lapply(t_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                 vars.to.regress = c("percent_mito", "nCount_RNA"))
t_features <- SelectIntegrationFeatures(t_list, nfeatures = 3000)
t_list <- PrepSCTIntegration(t_list, anchor.features = t_features)
t_list <- lapply(t_list, FUN = RunPCA, features = t_features, npcs = 100)
t_anchors <- FindIntegrationAnchors(t_list, normalization.method = "SCT",
                                    anchor.features = t_features, dims = 1:30,
                                    reduction = "rpca", k.anchor = 5)
t <- IntegrateData(anchorset = t_anchors, normalization.method = "SCT", 
                   dims = 1:30,
                   k.weight = 50)

t <- RunPCA(t, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(t, file = "data_objects/08_tcells.RDS")
message("~~~~~~~~~~~~~~~~~~~~Starting B lymphocytes~~~~~~~~~~~~~~~~~~~~")
b <- subset(obj, subset = cell_type == "B-cells" | 
              cell_type == "Plasma cells")

b_list <- SplitObject(b, split.by = "Batch")
b_list <- lapply(b_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                 vars.to.regress = c("percent_mito", "nCount_RNA"))
b_features <- SelectIntegrationFeatures(b_list, nfeatures = 3000)
b_list <- PrepSCTIntegration(b_list, anchor.features = b_features)
b_list <- lapply(b_list, FUN = RunPCA, features = b_features, npcs = 100)
b_anchors <- FindIntegrationAnchors(b_list, normalization.method = "SCT",
                                    anchor.features = b_features, dims = 1:30,
                                    reduction = "rpca", k.anchor = 5)
b <- IntegrateData(anchorset = b_anchors, normalization.method = "SCT", 
                   dims = 1:30,
                   k.weight = 50)

b <- RunPCA(b, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(b, file = "data_objects/08_blymph.RDS")

message("~~~~~~~~~~~~~~~~~~~~Starting microglia~~~~~~~~~~~~~~~~~~~~")

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

saveRDS(microglia, file = "data_objects/08_microglia.RDS")

message("~~~~~~~~~~~~~~~~~~~~Starting CD45- cells~~~~~~~~~~~~~~~~~~~~")
cd45neg <- subset(obj, subset = cell_type == "Astrocytes" |
                    cell_type == "Pericytes" | 
                    cell_type == "Endothelial cells" |
                    cell_type == "Oligodendrocytes" | 
                    cell_type == "Choroid plexus cells" | 
                    cell_type == "Leptomeningeal cells" | 
                    cell_type == "Stem cells" | 
                    cell_type == "Erythrocytes")

cd45neg_list <- SplitObject(cd45neg, split.by = "Batch")
cd45neg_list <- lapply(cd45neg_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                       vars.to.regress = c("percent_mito", "nCount_RNA"))
cd45neg_features <- SelectIntegrationFeatures(cd45neg_list, nfeatures = 3000)
cd45neg_list <- PrepSCTIntegration(cd45neg_list, anchor.features = cd45neg_features)
cd45neg_list <- lapply(cd45neg_list, FUN = RunPCA, features = cd45neg_features, npcs = 100)
cd45neg_anchors <- FindIntegrationAnchors(cd45neg_list, normalization.method = "SCT",
                                          anchor.features = cd45neg_features, dims = 1:30,
                                          reduction = "rpca", k.anchor = 5)
cd45neg <- IntegrateData(anchorset = cd45neg_anchors, normalization.method = "SCT", 
                         dims = 1:30)

cd45neg <- RunPCA(cd45neg, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(cd45neg, file = "data_objects/08_cd45neg.RDS")
message("~~~~~~~~~~~~~~~~~~~~Starting granulocytes~~~~~~~~~~~~~~~~~~~~")
gran <- subset(obj, subset = cell_type == "Basophils" | 
                 cell_type == "Neutrophils" | 
                 cell_type == "Mast cells")

gran_list <- SplitObject(gran, split.by = "Batch")
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

saveRDS(gran, file = "data_objects/08_granulocytes.RDS")
message("~~~~~~~~~~~~~~~~~~~~Starting monocytes~~~~~~~~~~~~~~~~~~~~")
mono <- subset(obj, subset = cell_type == "Monocytes")

mono_list <- SplitObject(mono, split.by = "Batch")
mono_list <- lapply(mono_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                    vars.to.regress = c("percent_mito", "nCount_RNA"))
mono_features <- SelectIntegrationFeatures(mono_list, nfeatures = 3000)
mono_list <- PrepSCTIntegration(mono_list, anchor.features = mono_features)
mono_list <- lapply(mono_list, FUN = RunPCA, features = mono_features, npcs = 100)
mono_anchors <- FindIntegrationAnchors(mono_list, normalization.method = "SCT",
                                       anchor.features = mono_features, dims = 1:30,
                                       reduction = "rpca", k.anchor = 5)
mono <- IntegrateData(anchorset = mono_anchors, normalization.method = "SCT", 
                      dims = 1:30)

mono <- RunPCA(mono, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(mono, file = "data_objects/08_monocytes.RDS")
message("~~~~~~~~~~~~~~~~~~~~Starting macrophages~~~~~~~~~~~~~~~~~~~~")
mac <- subset(obj, subset = cell_type == "Macrophages")

mac_list <- SplitObject(mac, split.by = "Batch")
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

saveRDS(mac, file = "data_objects/08_macrophages.RDS")

message("~~~~~~~~~~~~~~~~~~~~Starting dendritic cells~~~~~~~~~~~~~~~~~~~~")
dc <- subset(obj, subset = cell_type == "Dendritic cells")

dc_list <- SplitObject(dc, split.by = "Batch")
dc_list <- lapply(dc_list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2",
                  vars.to.regress = c("percent_mito", "nCount_RNA"))
dc_features <- SelectIntegrationFeatures(dc_list, nfeatures = 3000)
dc_list <- PrepSCTIntegration(dc_list, anchor.features = dc_features)
dc_list <- lapply(dc_list, FUN = RunPCA, features = dc_features, npcs = 100)
dc_anchors <- FindIntegrationAnchors(dc_list, normalization.method = "SCT",
                                     anchor.features = dc_features, dims = 1:30,
                                     reduction = "rpca", k.anchor = 5)
dc <- IntegrateData(anchorset = dc_anchors, normalization.method = "SCT", 
                    dims = 1:30)

dc <- RunPCA(dc, npcs = 100) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(dc, file = "data_objects/08_dendritic_cells.RDS")

message("~~~~~~~~~~~~~~~~~~~~All cells accounted for?~~~~~~~~~~~~~~~~~~~~")
nrow(obj@meta.data) == nrow(microglia@meta.data) + nrow(mono@meta.data) + 
  nrow(mac@meta.data) + nrow(cd45neg@meta.data) + nrow(t@meta.data) + 
  nrow(b@meta.data) + nrow(gran@meta.data) + nrow(nk@meta.data) + nrow(dc@meta.data)