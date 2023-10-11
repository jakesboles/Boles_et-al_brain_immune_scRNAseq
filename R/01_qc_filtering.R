library(Seurat)
library(scCustomize)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(janitor)
library(kableExtra)

#file path to save plots
plots <- "plots/"

#file path to save stats
stats <- "tabular_output/"

#ADD DATA LOADING OBJECT AFTER 00 IS DONE 

#Adding necessary sample metadata, including treatment groups----
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
  save_kable(file = "stats/raw_qc.html")

stats_clean[, c(2, 1, 3:ncol(stats_clean))] %>%
  kbl(format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "bordered")) %>%
  save_kable(file = "stats/filtered_qc.html")

#Integrate samples on batch with RPCA for normalization purposes