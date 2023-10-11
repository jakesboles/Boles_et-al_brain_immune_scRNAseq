library(Seurat)
library(scCustomize)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(janitor)
library(kableExtra)

#file path to save plots
dir.create("plots/01_qc_normalization")
plots <- "plots/01_qc_normalization"

#file path to save stats
stats <- "tabular_output/"

#Load raw data into Seurat objects and assemble ---- 
s1.data <- Read10X(data.dir = "raw_data/sample1/")
s1 <- CreateSeuratObject(counts = s1.data, min.cells = 3, min.features = 200,
                         project = "sample1")
s2.data <- Read10X(data.dir = "raw_data/sample2/")
s2 <- CreateSeuratObject(counts = s2.data, min.cells = 3, min.features = 200,
                         project = "sample2")
s3.data <- Read10X(data.dir = "raw_data/sample3/")
s3 <- CreateSeuratObject(counts = s3.data, min.cells = 3, min.features = 200,
                         project = "sample3")
s4.data <- Read10X(data.dir = "raw_data/sample4/")
s4 <- CreateSeuratObject(counts = s4.data, min.cells = 3, min.features = 200,
                         project = "sample4")
s5.data <- Read10X(data.dir = "raw_data/sample5/")
s5 <- CreateSeuratObject(counts = s5.data, min.cells = 3, min.features = 200,
                         project = "sample5")
s6.data <- Read10X(data.dir = "raw_data/sample6/")
s6 <- CreateSeuratObject(counts = s6.data, min.cells = 3, min.features = 200,
                         project = "sample6")
s7.data <- Read10X(data.dir = "raw_data/sample7/")
s7 <- CreateSeuratObject(counts = s7.data, min.cells = 3, min.features = 200,
                         project = "sample7")
s8.data <- Read10X(data.dir = "raw_data/sample8/")
s8 <- CreateSeuratObject(counts = s8.data, min.cells = 3, min.features = 200,
                         project = "sample8")
obj <- merge(s1, y = c(s2, s3, s4, s5, s6, s7, s8),
             add.cell.ids = c("sample1", "sample2", "sample3", "sample4",
                              "sample5", "sample6", "sample7", "sample8"),
             project = "PIPseq1")

#Adding necessary sample metadata, including treatment groups ----
meta <- data.frame(obj@meta.data) %>%
  rownames_to_column(var = "cell_name") %>%
  as_tibble()

BATCH <- c("sample1", "sample2", "sample3", "sample4")
LPS <- c("sample2", "sample4", "sample6", "sample8")
INHIBITORS <- c("sample3", "sample4", "sample7", "sample8")

meta <- meta %>%
  mutate(Dissociation = if_else(orig.ident %in% INHIBITORS, "Inhibitors", "DMSO"),
         Treatment = if_else(orig.ident %in% LPS, "LPS", "Saline"),
         Batch = if_else(orig.ident %in% BATCH, "Group1", "Group2")) %>%
  mutate(Group = case_when(Dissociation == "DMSO" & Treatment == "Saline" ~ "Saline_DMSO",
                           Dissociation == "DMSO" & Treatment == "LPS" ~ "LPS_DMSO",
                           Dissociation == "Inhibitors" & Treatment == "Saline" ~ "Saline_Inhibitors",
                           Dissociation == "Inhibitors"  & Treatment == "LPS" ~ "LPS_Inhibitors")) %>%
  as.data.frame() %>%
  column_to_rownames(var = "cell_name")

obj <- AddMetaData(obj, metadata = meta)

#QC on unfiltered data----
obj <- Add_Mito_Ribo_Seurat(obj, species = "mouse")
obj <- Add_Cell_Complexity_Seurat(obj)

p1 <- QC_Plots_Genes(obj, low_cutoff = 250, high_cutoff = 6500, y_axis_log = T,
                     plot_title = "Genes per cell") #nFeature
p2 <- QC_Plots_UMIs(obj, low_cutoff = 800, high_cutoff = 45000, y_axis_log = T,
                    plot_title = "UMIs per cell") #nCount
p3 <- QC_Plots_Mito(obj, high_cutoff = 5,
                    plot_title = "Mito genes per cell")
p4 <- QC_Plots_Complexity(obj, high_cutoff = 0.8)
(QC1 <- wrap_plots(p1, p2, p3, p4, ncol = 4))


#Filter on percent_mito and recheck QC metrics----
obj_clean <- subset(obj, percent_mito < 5)

#Filter on nCount and recheck QC metrics----
obj_clean <- subset(obj_clean, nCount_RNA > 800 & nCount_RNA < 45000)

#Filter on nFeature and recheck QC metrics----
obj_clean <- subset(obj_clean, nFeature_RNA > 250 & nFeature_RNA < 6500)

#Filter on cell complexity and recheck QC metrics----
obj_clean <- subset(obj_clean, log10GenesPerUMI > 0.8)

#Summary stats on final cleaned object----
med_stats_raw <- Median_Stats(obj, group_by_var = "orig.ident")
med_stats_clean <- Median_Stats(obj_clean, group_by_var = "orig.ident")

stats_raw <- obj@meta.data %>%
  as.data.frame() %>%
  group_by(orig.ident) %>%
  summarise(Cell_Count = n()) %>%
  adorn_totals("row") %>%
  dplyr::select(-1) %>%
  cbind(med_stats_raw)

stats_clean <- obj_clean@meta.data %>%
  as.data.frame() %>%
  group_by(orig.ident) %>%
  summarize(Cell_Count = n()) %>%
  adorn_totals("row") %>%
  dplyr::select(-1) %>%
  cbind(med_stas_clean)

stats_raw[, c(2, 1, 3:ncol(stats_raw))] %>%
  kbl(format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "bordered")) %>%
  save_kable(file = paste0(stats, "qc_raw.html"))

stats_clean[, c(2, 1, 3:ncol(stats_clean))] %>%
  kbl(format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "bordered")) %>%
  save_kable(file = paste0(stats, "qc_filtered.html"))

#Integrate samples on BATCH with RPCA for normalization purposes according to 
##https://satijalab.org/seurat/articles/integration_rpca.html
obj_list <- SplitObject(obj_clean, split.by = "Batch")
rm("obj_clean")
obj_list <- lapply(obj_list, FUN = SCTransform, method = "glmGamPoi")
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

ElbowPlot(intg1, ndims = 100)
ggsave(paste0(plots, "initial_PCA_elbow.png"), units = "in",
       dpi = 600, width = 6, height = 6)

intg1 <- RunUMAP(intg1, reduction = "pca", dims = 1:30, verbose = F) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

saveRDS(intg1, file = "data_objects/01_batch_integrated.RDS")