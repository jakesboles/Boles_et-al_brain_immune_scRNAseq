library(Seurat)
library(clusterProfiler)
library(tidyverse)
library(ggplot2)
library(scCustomize)
library(org.Mm.eg.db)
library(enrichplot)
library(ggupset)
library(ggtext)
library(camcorder)
library(stringr)
library(ggpubr)
library(ggrepel)
library(MAST)
library(ggtext)
library(scales)

load("data_objects/10_pseudobulk_deseq_results.RData")

list <- list()

for (i in seq_along(lfc_shrunk_results)){
  lfc_shrunk_results[[i]] %>%
    arrange(desc(log2FoldChange)) %>%
    rownames_to_column(var = "X") %>%
    mutate(Entrez = mapIds(org.Mm.eg.db, X, 'ENTREZID', 'SYMBOL')) %>%
    na.omit() %>%
    pull(Entrez) -> gene_ids
  
  lfc_shrunk_results[[i]] %>%
    arrange(desc(log2FoldChange)) %>%
    rownames_to_column(var = "X") %>%
    mutate(Entrez = mapIds(org.Mm.eg.db, X, 'ENTREZID', 'SYMBOL')) %>%
    na.omit() %>%
    pull(log2FoldChange) -> log2fc
  
  names(log2fc) <- gene_ids
  
  list[[i]] <- log2fc
  
  names(list)[i] <- names(lfc_shrunk_results)[i]
  
}

res <- compareCluster(list,
                      fun = "gseGO",
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      minGSSize = 30,
                      maxGSSize = 300,
                      pvalueCutoff = 1)
names(list)

res <- arrange(res, p.adjust)
res[grep("Microglia", res@compareClusterResult$Cluster), 1:8] %>% head(30L)
pathways1 <- c("GO:0002181", "GO:0015986", "GO:0000398", "GO:0030593", "GO:0140053")

res[grep("Macrophages", res@compareClusterResult$Cluster), 1:8] %>% head(30L)
pathways2 <- c("GO:0001909", "GO:0070555", "GO:0032496", "GO:0001906")

res[grep("Neutrophils", res@compareClusterResult$Cluster), 1:8] %>% head(30L)
pathways3 <- c("GO:0000819", "GO:0071706", "GO:0002221")

res[grep("Mast cells", res@compareClusterResult$Cluster), 1:8] %>% head(30L)
#no significant pathways

res[grep("Basophils", res@compareClusterResult$Cluster), 1:8] %>% head(30L)
#no significant pathways

res[grep("B-cells", res@compareClusterResult$Cluster), 1:8] %>% head(30L)
pathways6 <- c("GO:0042113", "GO:0007015", "GO:0050863")

res[grep("Plasma cells", res@compareClusterResult$Cluster), 1:8] %>% head(30L)
pathways7 <- c("GO:0016064", "GO:0002250")

res[grep("Dendritic cells", res@compareClusterResult$Cluster), 1:8] %>% head(30L)
pathways8 <- c("GO:0006909", "GO:0042060")

res[grep("Monocytes", res@compareClusterResult$Cluster), 1:8] %>% head(30L)
pathways9 <- c("GO:0034341", "GO:0072593", "GO:0034612")

res[grep("Natural killer cells", res@compareClusterResult$Cluster), 1:8] %>% head(30L)
#no new pathways from top 30

res[grep("αβ T-cells", res@compareClusterResult$Cluster), 1:8] %>% head(30L)
pathways11 <- c("GO:0035456", "GO:0007159", "GO:0032606")

res[grep("γδ T-cells", res@compareClusterResult$Cluster), 1:8] %>% head(30L)
pathways12 <- c("GO:0098609")

pathways <- c(pathways1, pathways2, pathways3, pathways6, pathways7, 
              pathways8, pathways9, pathways11, pathways12)

firstup <- function(x){
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

res@compareClusterResult %>%
  as.data.frame() %>%
  filter(ID %in% pathways) %>%
  pivot_wider(names_from = "Cluster",
              values_from = "p.adjust",
              id_cols = c("Description")) %>%
  column_to_rownames(var = "Description") %>%
  mutate_at(vars(Microglia:`Mast cells`),
            function(x) if_else(x < 0.05, 1, 0)) %>%
  mutate(rank = rowSums(., na.rm = T)) %>%
  arrange(desc(rank)) %>%
  rownames() %>%
  firstup() -> order

res@compareClusterResult %>%
  as.data.frame() %>%
  filter(ID %in% pathways) %>%
  mutate(size = if_else(p.adjust < 0.05, 1, 2) %>%
           factor(levels = c(2, 1),
                  labels = c("BH-adjusted\np > 0.05", "BH-adjusted\np < 0.05"))) %>%
  ggplot(aes(x = Cluster, y = factor(firstup(Description),
                                     levels = order), color = NES, size = size)) +
  geom_point() + 
  scale_color_gradient2(low = "midnightblue",
                        high = "orange3") + 
  scale_y_discrete(limits = rev) +
  scale_x_discrete(limits = sort(names(list))) +
  labs(size = "",
       color = "Normalized\nenrichment\nscore") +
  theme_bw(base_size = 10) + 
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9),
    axis.title = element_blank()
  ) 
ggsave("plots/Fig_5C_multi_gsea.png",
       units = "in", dpi = 600,
       height = 6, width = 7)
