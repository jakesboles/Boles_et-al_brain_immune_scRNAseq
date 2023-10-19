library(Seurat)
library(tidyverse)

#B-cell dataset----
message("~~~~~~~~~~~~~~~~~~~~Preparing B-cells object~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
bcells <- readRDS("data_objects/06_annotated_bcells.RDS")

Idents(bcells) <- bcells@meta.data$pca.SCT_res.0.6 

DefaultAssay(bcells) <- "SCT"

ids <- c("B-cells", "B-cells", "B-cells",
         "B-cells", "B-cells", "Doublets")

levels(bcells)
levels <- as.character(c(0:max(as.numeric(levels(bcells)))))
names(ids) <- levels

bcells <- RenameIdents(bcells, ids)

levels(bcells)

bcells[["cell_type"]] <- Idents(bcells)
#Monocytes/macrophage dataset----
message("~~~~~~~~~~~~~~~~~~~~Preparing monocytes/macrophages object~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
monomacs <- readRDS("data_objects/06_annotated_mono-macs.RDS")

Idents(monomacs) <- monomacs@meta.data$pca.SCT_res.0.8 #greatest median silhouette score

DefaultAssay(monomacs) <- "SCT"

ids <- c("Monocytes", "Macrophages", "Macrophages", "Macrophages", "Monocytes",
         "Monocytes", "Macrophages", "Dendritic cells", "Macrophages", "Macrophages",
         "Monocytes", "Dendritic cells", "Dendritic cells", "Monocytes", "Monocytes",
         "Plasma cells", "Monocytes", "Erythrocytes")

levels(monomacs)
levels <- as.character(c(0:max(as.numeric(levels(monomacs)))))
names(ids) <- levels

monomacs <- RenameIdents(monomacs, ids)

monomacs[["cell_type"]] <- Idents(monomacs)
#Neutrophils----
message("~~~~~~~~~~~~~~~~~~~~Preparing neutrophils object~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
neutrophils <- readRDS("data_objects/06_annotated_neutrophils.RDS")

Idents(neutrophils) <- neutrophils@meta.data$pca.SCT_res.0.6 #greatest median silhouette score

DefaultAssay(neutrophils) <- "SCT"

levels <- as.character(c(0:max(as.numeric(levels(neutrophils)))))
ids <- c(rep("Neutrophils", 8), "Doublets", "Neutrophils")
names(ids) <- levels

neutrophils <- RenameIdents(neutrophils, ids)

neutrophils[["cell_type"]] <- Idents(neutrophils)
#T/NK cells----
message("~~~~~~~~~~~~~~~~~~~~Preparing T-cells/NK cells object~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
tnk <- readRDS("data_objects/06_annotated_tnk.RDS")

Idents(tnk) <- tnk@meta.data$pca.SCT_res.0.6 

DefaultAssay(tnk) <- "SCT"

ids <- c("αβ T-cells", "αβ T-cells", "αβ T-cells", "Natural killer cells",
         "Natural killer cells", "Stem cells", "γδ T-cells", "Basophils", "Macrophages",
         "Stem cells", "Mast cells", "Erythrocytes")


levels <- as.character(c(0:max(as.numeric(levels(tnk)))))
names(ids) <- levels

tnk <- RenameIdents(tnk, ids)
levels(tnk)

tnk[["cell_type"]] <- Idents(tnk)
#CD45neg cells----
message("~~~~~~~~~~~~~~~~~~~~Preparing non-immune object~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
cd45neg <- readRDS("data_objects/06_annotated_cd45neg.RDS")

Idents(cd45neg) <- cd45neg@meta.data$pca.SCT_res.0.6 

DefaultAssay(cd45neg)

ids <- c("Choroid plexus cells", "Endothelial cells", "Astrocytes", "Endothelial cells",
         "Endothelial cells", "Pericytes", "Choroid plexus cells", "Astrocytes",
         "Leptomeningeal cells", "Pericytes", "Leptomeningeal cells", "Astrocytes",
         "Leptomeningeal cells", "Choroid plexus cells", "Oligodendrocytes",
         "Oligodendrocytes", "Astrocytes")

levels <- as.character(c(0:max(as.numeric(levels(cd45neg)))))
names(ids) <- levels

cd45neg <- RenameIdents(cd45neg, ids)

cd45neg[["cell_type"]] <- Idents(cd45neg)
#Load microglia too----
message("~~~~~~~~~~~~~~~~~~~~Preparing microglia object~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
microglia <- readRDS("data_objects/06_annotated_microglia.RDS")

microglia <- FindClusters(microglia,
                          resolution = 0.4) 

Idents(microglia) <- "Microglia"
microglia[["cell_type"]] <- Idents(microglia)

#Create object----
message("~~~~~~~~~~~~~~~~~~~~Loading original full object and extracting metadata~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
obj <- readRDS("data_objects/04_coarse_annotated.RDS")
meta <- obj@meta.data
meta <- meta %>%
  rownames_to_column(var = "cell")
colnames(meta)
meta <- meta[, c(1:14)]
meta_cols <- colnames(meta)

message("~~~~~~~~~~~~~~~~~~~~Creating annotated metadata frame~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
b_meta <- bcells@meta.data %>%
  rownames_to_column(var = "cell") %>%
  dplyr::select(c("cell", "cell_type"))
tnk_meta <- tnk@meta.data %>%
  rownames_to_column(var = "cell") %>%
  dplyr::select(c("cell", "cell_type"))
monomacs_meta <- monomacs@meta.data %>%
  rownames_to_column(var = "cell") %>%
  dplyr::select(c("cell", "cell_type"))
neutrophils_meta <- neutrophils@meta.data %>%
  rownames_to_column(var = "cell") %>%
  dplyr::select(c("cell", "cell_type"))
cd45neg_meta <- cd45neg@meta.data %>%
  rownames_to_column(var = "cell") %>%
  dplyr::select(c("cell", "cell_type"))
microglia_meta <- microglia@meta.data %>%
  rownames_to_column(var = "cell") %>%
  dplyr::select(c("cell", "cell_type"))

meta <- meta %>%
  left_join(b_meta, by = "cell") %>%
  left_join(tnk_meta, by = "cell") %>%
  mutate(cell_type = coalesce(cell_type.x, cell_type.y),
         .keep = "unused") %>%
  left_join(monomacs_meta, by = "cell") %>%
  mutate(cell_type = coalesce(cell_type.x, cell_type.y),
         .keep = "unused") %>%
  left_join(neutrophils_meta, by = "cell") %>%
  mutate(cell_type = coalesce(cell_type.x, cell_type.y),
         .keep = "unused") %>%
  left_join(cd45neg_meta, by = "cell") %>%
  mutate(cell_type = coalesce(cell_type.x, cell_type.y),
         .keep = "unused") %>%
  left_join(microglia_meta, by = "cell") %>%
  mutate(cell_type = coalesce(cell_type.x, cell_type.y),
         .keep = "unused") %>%
  column_to_rownames(var = "cell")

message("~~~~~~~~~~~~~~~~~~~~All cells accounted for?~~~~~~~~~~~~~~~~~~~~")
nrow(meta) == nrow(obj@meta.data)
message("~~~~~~~~~~~~~~~~~~~~Any cells without annotation?~~~~~~~~~~~~~~~~~~~~")
sum(is.na(meta$cell_type))

message("~~~~~~~~~~~~~~~~~~~~Adding metadata to full object~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
obj@meta.data <- meta
message("~~~~~~~~~~~~~~~~~~~~Removing additional doublet cluster~~~~~~~~~~~~~~~~~~~~")
message("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
obj <- subset(obj, subset = cell_type == "Doublets", invert = T)
message("~~~~~~~~~~~~~~~~~~~~TOTAL CELLS:~~~~~~~~~~~~~~~~~~~~")
obj
#Re-integrate, normalize, etc.----
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

saveRDS(intg1,
        file = "data_objects/07_annotated_object.rds")