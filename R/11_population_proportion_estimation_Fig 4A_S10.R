#This script uses the full, cleaned, annotated object and the package Speckle 
#to estimate the differences in cell proportions between fixed variables of interest.
#Plots generated here accompany our pre-print's Figure 4 and Figure S10. 

library(speckle)
library(SingleCellExperiment)
library(CellBench)
library(limma)
library(ggplot2)
library(scater)
library(patchwork)
library(edgeR)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(Seurat)
library(tidyverse)
library(ggbeeswarm)
library(scCustomize)
library(camcorder)

obj <- readRDS("data_objects/07_annotated_object.rds")
Idents(obj)

obj$sample <- obj$orig.ident
obj$group <- obj$Treatment
obj$Group <- NULL

stats <- propeller(obj, 
                   sample = obj$orig.ident,
                   group = obj$Treatment,
                   transform = "asin")
meta <- obj@meta.data

df <- meta %>%
  group_by(orig.ident) %>%
  mutate(total_cells = n()) %>%
  ungroup() %>%
  group_by(orig.ident, cell_type) %>%
  mutate(fraction = (n() / total_cells) * 100) %>%
  dplyr::select(c(orig.ident, Treatment, cell_type, fraction)) %>%
  distinct(orig.ident, cell_type, .keep_all = T) %>%
  pivot_wider(names_from = "cell_type",
              values_from = "fraction")

df[is.na(df)] <- 0

pal <- DiscretePalette_scCustomize(num_colors = 20,
                                   palette = "varibow")

nonimmune_pal <- setdiff(pal, immune_pal)

immune_pal <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF',
                '#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A')

df %>%
  dplyr::select(-matches("Astro|Choroid|Endot|Eryt|Lepto|Oligo|Peri|Stem")) %>%
  pivot_longer(!c("orig.ident", "Treatment"),
               names_to = "BaselineProp.clusters",
               values_to = "fraction") %>%
  left_join(stats, by = "BaselineProp.clusters") %>%
  .[, c(1:4, 11)] %>%
  mutate(title = if_else(FDR > 0.001,
                         paste0(BaselineProp.clusters, "\nFDR = ", round(FDR, digits = 3)),
                         paste0(BaselineProp.clusters, "\nFDR < 0.001"))) %>%
  ggplot(aes(x = factor(Treatment,
                        levels = c("Saline", "LPS")),
             y = fraction)) +
  stat_summary(color = "grey60",
               fun = "mean", geom = "crossbar") +
  stat_summary(color = "grey60",
               fun.data = "mean_se", geom = "errorbar",
               linewidth = 1.2,
               width = 0.5) +
  geom_quasirandom(aes(fill = BaselineProp.clusters),
                   color = "black",
                   shape = 21,
                   size = 5,
                   alpha = 0.8,
                   show.legend = F) + 
  scale_fill_manual(values = immune_pal) + 
  ylim(0, NA) +
  labs(y = "% of sample") +
  facet_wrap(vars(title),
             ncol = 6,
             scales = "free_y") + 
  theme_bw(base_size = 14) + 
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    strip.text = element_text(color = "white", face = "bold"),
    strip.background = element_rect(fill = "black")
  )
ggsave("plots/immune_cell_frequencies_Fig_4A.png",
       units = "in", dpi = 600,
       height = 4, width = 12)

df %>%
  dplyr::select(matches("Astro|Choroid|Endot|Eryt|Lepto|Oligo|Peri|Stem|orig|Treat")) %>%
  pivot_longer(!c("orig.ident", "Treatment"),
               names_to = "BaselineProp.clusters",
               values_to = "fraction") %>%
  left_join(stats, by = "BaselineProp.clusters") %>%
  .[, c(1:4, 11)] %>%
  mutate(title = if_else(FDR > 0.001,
                         paste0(BaselineProp.clusters, "\nFDR = ", round(FDR, digits = 3)),
                         paste0(BaselineProp.clusters, "\nFDR < 0.001"))) %>%
  ggplot(aes(x = factor(Treatment,
                        levels = c("Saline", "LPS")),
             y = fraction)) +
  stat_summary(color = "grey60",
               fun = "mean", geom = "crossbar") +
  stat_summary(color = "grey60",
               fun.data = "mean_se", geom = "errorbar",
               linewidth = 1.2,
               width = 0.5) +
  geom_quasirandom(aes(fill = BaselineProp.clusters),
                   color = "black",
                   shape = 21,
                   size = 5,
                   alpha = 0.6,
                   show.legend = F) + 
  scale_fill_manual(values = nonimmune_pal) + 
  ylim(0, NA) +
  labs(y = "% of sample") +
  facet_wrap(vars(title),
             ncol = 4,
             scales = "free_y") + 
  theme_bw(base_size = 14) + 
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    strip.text = element_text(color = "white", face = "bold"),
    strip.background = element_rect(fill = "black")
  )
ggsave("plots/nonimmune_cell_frequencies_Fig_S10.png",
       units = "in", dpi = 600,
       height = 4, width = 10)