library(Seurat)
library(scCustomize)
library(ggplot2)
library(tidyverse)
library(AnnotationDbi)
library(Orthology.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(nichenetr)
library(rlist)
library(camcorder)
library(readxl)

load("data_objects/pseudobulk_deseq_results.RData")

ad <- c("APOE", "BIN1", "MS4A6A", "PICALM", "CR1", "TREM2", "ABCA7", "PTK2B",
        "PLCG2", "SPI1", "SORL1", "CD2AP", "SLC24A4", "RIN3", "CASS4", "ADAM10",
        "HAVCR2", "SCIMP", "CLNK", "TNIP1", "ABCA1", "USP6NL", "INPP5D",
        "CD33", "WWOX", "ABI3", "APH1B", "HS3ST1", "CHRNE", "CCDC6", "KAT8")

gwas <- read_xlsx("Nalls_et_al_2019_snps.xlsx",
                  range = "A1:AE108")

pd_ms <- gwas %>%
  arrange(`P, random effects, all studies`) %>%
  pull(5) %>%
  convert_human_to_mouse_symbols() %>%
  na.omit() %>%
  unique()

ad_ms <- convert_human_to_mouse_symbols(ad) %>%
  na.omit() 

for (i in seq_along(deseq_results)){
  deseq_results[[i]] <- deseq_results[[i]] %>%
    as.data.frame() %>%
    mutate(cell = names(deseq_results)[i]) %>%
    rownames_to_column(var = "gene")
}

df <- list.rbind(deseq_results)

gg_record(device = "png",
          units = "in",
          dpi = 600,
          height = 4, width = 9)

df %>%
  filter(gene %in% ad_ms) %>%
  mutate(sig = case_when(padj >= 0.05 ~ "",
                         padj < 0.05 & padj >= 0.01 ~ "*",
                         padj < 0.01 & padj >= 0.001 ~ "**",
                         padj < 0.001 ~ "***")) %>%
  ggplot(aes(x = factor(gene, levels = ad_ms), y = cell)) +
  geom_tile(aes(fill = log2FoldChange),
            color = "grey60") + 
  geom_text(aes(label = sig),
            color = "black") +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(high = "orange3", low = "midnightblue", mid = "white") + 
  guides(fill = guide_colorbar(title.position = "left")) +
  theme_void() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "italic"),
        axis.title = element_blank(),
        axis.text.y = element_text(color = "black", hjust = 1),
        legend.position= "right",
        legend.title = element_text(angle = 90,
                                    hjust = 0.5))
ggsave("05_plots/Fig_5D_gwas_hits_ad.png",
       units = "in", dpi = 600,
       height = 4, width = 9) 

gg_record(device = "png",
          units = "in",
          dpi = 600,
          height = 4, width = 16)

df %>%
  filter(gene %in% pd_ms) %>%
  mutate(sig = case_when(padj >= 0.05 ~ "",
                         padj < 0.05 & padj >= 0.01 ~ "*",
                         padj < 0.01 & padj >= 0.001 ~ "**",
                         padj < 0.001 ~ "***")) %>%
  ggplot(aes(x = factor(gene, levels = pd_ms), y = cell)) +
  geom_tile(aes(fill = log2FoldChange),
            color = "grey60") + 
  geom_text(aes(label = sig),
            color = "black") +
  scale_y_discrete(limits = rev) +
  guides(fill = guide_colorbar(title.position = "left")) +
  scale_fill_gradient2(high = "orange3", low = "midnightblue", mid = "white") + 
  theme_void() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "italic"),
        axis.title = element_blank(),
        axis.text.y = element_text(color = "black", hjust = 1),
        legend.position= "right",
        legend.title = element_text(angle = 90,
                                    hjust = 0.5))
ggsave("plots/Fig_5E_gwas_hits_pd.png",
       units = "in", dpi = 600,
       height = 4, width = 16)        

