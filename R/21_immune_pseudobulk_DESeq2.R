library(tidyverse)
library(cowplot)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(Seurat)
library(scCustomize)
library(patchwork)
library(paletteer)

obj <- readRDS("data_objects/09_immune_object.RDS")

bulk <- AggregateExpression(obj,
                            assays = "RNA",
                            group.by = c("orig.ident", "cell_type"),
                            return.seurat = F,
                            slot = "counts")

exp <- bulk$RNA

cluster_names <- unique(obj$cell_type)
sample_names <- unique(obj$orig.ident)

tstrsplit(colnames(exp), "_") %>%
  str()

counts_list <- list()

for (i in 1:length(cluster_names)){
  column_idx <- which(tstrsplit(colnames(exp), "_")[[2]] == cluster_names[i])
  
  counts_list[[i]] <- exp[, column_idx]
  names(counts_list)[i] <- cluster_names[i]
}

meta <- obj@meta.data %>%
  as.data.frame() %>%
  dplyr::select(orig.ident, Treatment, Dissociation) %>%
  distinct()

rownames(meta) <- meta$orig.ident

t <- table(obj@meta.data$orig.ident,
           obj@meta.data$cell_type)
t

meta_list <- list()

for (i in seq_along(counts_list)){
  df <- data.frame(cluster_sample_id = colnames(counts_list[[i]]))
  
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  df$orig.ident <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  cell_counts <- cell_counts[cell_counts > 0]
  
  sample_order <- match(df$orig.ident, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  df$cell_count <- cell_counts
  
  df <- left_join(df, meta,
                  by = c("orig.ident"))
  
  rownames(df) <- df$cluster_sample_id
  
  meta_list[[i]] <- df
  names(meta_list)[i] <- unique(df$cluster_id)
  
  meta_list[[i]]$Treatment <- factor(meta_list[[i]]$Treatment,
                                     levels = c("Saline", "LPS"))
}

str(meta_list)

names(counts_list) == names(meta_list)

up_degs <- list()
down_degs <- list()
deseq_results <- list()
lfc_shrunk_results <- list()
for (i in seq_along(cluster_names)){
  idx <- which(names(counts_list) == cluster_names[i])
  cluster_counts <- counts_list[[idx]]
  cluster_metadata <- meta_list[[idx]]
  
  dds <- DESeqDataSetFromMatrix(cluster_counts,
                                colData = cluster_metadata,
                                design = ~ Treatment)
  dds <- DESeq(dds)
  
  res <- results(dds)
  
  res %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene") %>%
    filter(log2FoldChange > 1 & padj < 0.05) -> up_degs[[i]]
  
  res %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene") %>%
    filter(log2FoldChange < -1 & padj < 0.05) -> down_degs[[i]]
  
  names(up_degs)[i] <- cluster_names[i] -> names(down_degs)[i]
  
  deseq_results[[i]] <- res
  names(deseq_results)[i] <- cluster_names[i]
  
  lfc_shrunk_results[[i]] <- as.data.frame(lfcShrink(dds, res = res,
                                                     coef = "Treatment_LPS_vs_Saline",
                                                     type = "apeglm"))
  names(lfc_shrunk_results)[i] <- cluster_names[i]
}

save(deseq_results, lfc_shrunk_results, 
     file = "data_objects/10_pseudobulk_deseq_results.RData")