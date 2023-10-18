library(Seurat)
library(tidyverse)

reduction <- "pca"
assay <- "SCT"
choice <- 0.4
results_path <- "chooseR/microglia/"

# Load in the object containing the clustered results
obj <- readRDS(paste0(results_path, "clustered_data.rds"))
Idents(obj) <- obj@meta.data$pca.SCT_res.0.4#UPDATE THIS

# First is a cluster average co-clustering heatmap
# Read the data
grp <- readRDS(paste0(results_path, "frequency_grouped_", choice, ".rds"))

# As the data is symmetrical, we do not need the upper triangle
grp <- grp %>%
  pivot_wider(names_from = "cell_2", values_from = "avg_percent") %>%
  select(str_sort(colnames(.), numeric = T)) %>%
  mutate(cell_1 = as.numeric(as.character(cell_1))) %>%
  arrange(cell_1) %>%
  column_to_rownames("cell_1")
grp[lower.tri(grp)] <- NA
grp <- grp %>%
  as_tibble(rownames = "cell_1") %>%
  pivot_longer(-cell_1, names_to = "cell_2", values_to = "avg_percent") %>%
  mutate_at("cell_2", ordered, levels = unique(.$cell_1)) %>%
  mutate_at("cell_1", ordered, levels = unique(.$cell_1))

# And plot!
plot <- ggplot(grp, aes(cell_1, cell_2, fill = avg_percent)) +
  geom_tile(color = "white") +
  scale_x_discrete("Cluster", expand = c(0, 0)) +
  scale_y_discrete(
    "Cluster",
    limits = rev(levels(grp$cell_2)),
    expand = c(0, 0)
  ) +
  scale_fill_distiller(
    "% co-clustering",
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    palette = "RdYlBu",
    na.value = "white"
  ) +
  ggtitle("Average co-clustering heatmap") +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.position = c(0.9, 0.5),
    legend.title = element_text(size = 8, hjust = 0),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  ) +
  guides(fill = guide_colorbar(barheight = 3, barwidth = 1))
plot

ggsave(
  plot = plot,
  filename = paste0(results_path, "Fig_S6C.png"),
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)

# Let's add the silhouette scores to the Seurat object!
sil_scores <- readRDS(paste0(results_path, "silhouette_", choice, ".rds"))
sil_scores <- as.data.frame(sil_scores[, 3], row.names = Seurat::Cells(obj))
colnames(sil_scores) <- c("sil_score")
obj <- AddMetaData(obj, metadata = sil_scores)

# We also find it useful to visualise the silhouette scores on the UMAP!
plot <- FeaturePlot(
  obj,
  "sil_score",
  reduction = "umap",
  pt.size = 1,
  min.cutoff = -1,
  max.cutoff = 1
) +
  scale_colour_distiller(
    palette = "RdYlBu",
    labels = c(-1, 0, 1),
    breaks = c(-1, 0, 1),
    limits = c(-1, 1)
  ) + 
  ggtitle("Silhouette score") +
  theme(axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 axis.line = element_line(arrow = arrow(angle = 15, 
                                                        length = unit(0.5, "cm"), 
                                                        type = "closed")))
plot

ggsave(
  plot = plot,
  filename = paste0(plots, "Fig_S6B.png"),
  dpi = 300,
  height = 5,
  width = 5.5,
  units = "in"
)

#Below we'll generate a few plots for Supplemental Figure 6
theme <- theme(axis.text = element_blank(),   #thematic elements for plots
               axis.ticks = element_blank(),
               axis.line = element_line(arrow = arrow(angle = 15, 
                                                      length = unit(0.5, "cm"), 
                                                      type = "closed")))

#Labeled DimPlot in downsampled data
downsamp <- obj
Idents(downsamp) <- downsamp@meta.data$pca.SCT_res.0.4

DimPlot_scCustom(downsamp,
                 label = T,
                 colors_use = ColorBlind_Pal(),
                 label.size = 8) +
  theme
ggsave(paste0(results_path, "Fig_S6A.png"),
       units = "in", dpi = 600,
       height = 5, width = 5.5)

#Labeled DimPlot in full data
obj <- readRDS("data_objects/05_clean_microglia_sub.RDS")

DefaultAssay(obj) <- "integrated"

obj <- FindClusters(obj, resolution = 0.4)

DimPlot_scCustom(obj, label = T,
                 pt.size = 1,
                 label.size = 8,
                 repel = T) +
  theme
ggsave(paste0(results_path, "Fig_S6E.png"),
       units = "in", dpi = 600,
       height = 5, width = 5.5)

#Compute dice coefficient between cluster assignment in full and downsampled data
obj@meta.data -> full_meta
downsamp@meta.data -> sampled_meta

full_meta %>%
  dplyr::select(c(integrated_snn_res.0.4)) %>%
  rownames_to_column(var = "cell") -> full_meta
sampled_meta %>%
  dplyr::select(pca.SCT_res.0.4) %>%
  rownames_to_column(var = "cell") -> sampled_meta

sampled_meta %>%
  left_join(full_meta, by = "cell") %>%
  group_by(pca.SCT_res.0.4, integrated_snn_res.0.4) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  complete(pca.SCT_res.0.4, integrated_snn_res.0.4, fill = list("count" = 0)) %>%
  group_by(pca.SCT_res.0.4) %>%
  mutate(n_sampled = sum(count)) %>%
  ungroup() %>%
  group_by(integrated_snn_res.0.4) %>%
  mutate(n_full = sum(count)) %>%
  ungroup() %>%
  mutate(dice = (2 * count) / (n_full + n_sampled)) %>%
  mutate(
    pca.SCT_res.0.4 = factor(pca.SCT_res.0.4, 
                             levels = c(0:6)),
    integrated_snn_res.0.4 = factor(integrated_snn_res.0.4, 
                                    levels = c(0:13))) %>%
  ggplot(aes(y = pca.SCT_res.0.4, x = integrated_snn_res.0.4)) + 
  geom_tile(color = "white",
            aes(fill = dice)) + 
  scale_y_discrete("Down-sampled chooseR clusters", expand = c(0, 0)) +
  scale_x_discrete(
    "Equivalent resolution in full dataset",
    expand = c(0, 0)) +
  scale_fill_paletteer_c("grDevices::Inferno",
                         direction = 1) +
  labs(fill = "Dice\ncoefficient") +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8, hjust = 0),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )
ggsave(filename = paste0(results_path, "Fig_S6D.png"),
       dpi = 300,
       height = 2.5,
       width = 3.5,
       units = "in")
