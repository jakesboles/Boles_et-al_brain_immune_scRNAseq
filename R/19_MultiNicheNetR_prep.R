library(Seurat)
library(SingleCellExperiment)
library(multinichenetr)
library(ggplot2)
library(tidyverse)

obj <- readRDS("data_objects/09_immune_object.RDS")

#Load NicheNet database
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
lr_network <- lr_network %>% 
  dplyr::rename(ligand = from, receptor = to) %>% 
  distinct(ligand, receptor) %>% 
  mutate(ligand = make.names(ligand), receptor = make.names(receptor))

ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
colnames(ligand_target_matrix) <- colnames(ligand_target_matrix) %>% 
  make.names()
rownames(ligand_target_matrix) <- rownames(ligand_target_matrix) %>% 
  make.names()

obj$cell_type <- make.names(obj$cell_type)
unique(obj$cell_type)

sce <- Seurat::as.SingleCellExperiment(obj,
                                       assay = "SCT")
sce <- alias_to_symbol_SCE(sce, "mouse") %>%
  makenames_SCE()

sample_id <- "orig.ident"
group_id <- "Treatment"
celltype_id <- "cell_type"
covariates <- NA
batches <- NA

senders_oi <- SummarizedExperiment::colData(sce)[,celltype_id] %>% 
  unique() -> receivers_oi

min_cells <- 5

contrasts_oi <- c("'LPS-Saline','Saline-LPS'")
contrast_tbl <- tibble(
  contrast = c("LPS-Saline","Saline-LPS"),
  group = c("LPS", "Saline"))

logFC_threshold = 0.50
p_val_threshold = 0.05
fraction_cutoff = 0.05
top_n_target = 250
p_val_adj = FALSE 
empirical_pval = FALSE
cores_system = 12
n.cores = min(cores_system, union(senders_oi, receivers_oi) %>% length())

prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)

prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)

multinichenet_output = multi_nichenet_analysis(sce = sce, 
                                               celltype_id = celltype_id, 
                                               sample_id = sample_id, 
                                               group_id = group_id, 
                                               lr_network = lr_network, 
                                               ligand_target_matrix = ligand_target_matrix, 
                                               contrasts_oi = contrasts_oi, 
                                               contrast_tbl = contrast_tbl, 
                                               batches = batches, 
                                               covariates = covariates,
                                               prioritizing_weights = prioritizing_weights, 
                                               min_cells = min_cells, 
                                               logFC_threshold = logFC_threshold, 
                                               p_val_threshold = p_val_threshold,  
                                               fraction_cutoff = fraction_cutoff, 
                                               p_val_adj = p_val_adj, 
                                               empirical_pval = empirical_pval, 
                                               top_n_target = top_n_target, 
                                               n.cores = n.cores, 
                                               sender_receiver_separate = FALSE, 
                                               verbose = TRUE)

saveRDS(multinichenet_output,
        file = "data_objects/11_multinichenetr.RDS")
