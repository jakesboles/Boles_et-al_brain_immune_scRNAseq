library(UpSetR)
library(tidyverse)
library(scCustomize)
library(DESeq2)
library(camcorder)
library(rstudioapi)
library(ComplexUpset)
library(ggplot2)
library(rlist)

load("data_objects/10_pseudobulk_deseq_results.RData")

up <- list()
down <- list()

padj_threshold <- 0.01

for (i in seq_along(deseq_results)){
  deseq_results[[i]] %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene") %>%
    dplyr::filter(log2FoldChange > 2 & padj < padj_threshold) %>%
    dplyr::pull(gene) -> up[[i]]
  
  deseq_results[[i]] %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene") %>%
    dplyr::filter(log2FoldChange < -2 & padj < padj_threshold) %>%
    dplyr::pull(gene) -> down[[i]]
  
  names(up)[i] <- names(deseq_results)[i] -> names(down)[i]
}

unique_degs_plot <- function(list){
  unique <- list()
  shared <- list()
  
  for (i in seq_along(list)){
    oi <- list[[i]]
    bg <- list.rbind(list[-i])
    
    unique[[i]] <- setdiff(oi, bg)
    shared[[i]] <- intersect(oi, bg)
  }
  
  length.total <- list()
  length.unique <- list()
  length.shared <- list()
  
  fraction.unique <- list()
  fraction.shared <- list()
  
  for (i in seq_along(list)){
    length.total[[i]] <- length(list[[i]])
    length.unique[[i]] <- length(unique[[i]])
    length.shared[[i]] <- length(shared[[i]])
    
    fraction.unique[[i]] <- length.unique[[i]] / length.total[[i]]
    fraction.shared[[i]] <- length.shared[[i]] / length.total[[i]]
  }
  
  names(length.total) <- names(list)
  names(fraction.unique) <- names(list)
  names(fraction.shared) <- names(list)
  
  d <- data.frame(num_degs = unlist(length.total),
                  frac_unique = unlist(fraction.unique),
                  frac_shared = unlist(fraction.shared))
  
  p <- d %>%
    arrange(num_degs) %>%
    rownames_to_column(var = "cell") %>%
    mutate(cell = fct_inorder(cell)) %>%
    pivot_longer(c(frac_unique, frac_shared),
                 names_to = "uniqueness",
                 values_to = "fraction") %>%
    mutate(uniqueness = factor(uniqueness,
                               level = c("frac_shared", "frac_unique"),
                               labels = c("Shared", "Unique"))) %>%
    mutate(fraction = fraction * num_degs) %>%
    ggplot(aes(y = cell, x = fraction)) +
    geom_bar(aes(fill = uniqueness),
             stat = "identity") + 
    labs(x = paste0("Total ", deparse(substitute(list)), "-regulated genes")) +
    scale_x_continuous(expand = c(0, 0),
                       position = "top",
                       limits = c(0, ymax),
                       breaks = seq(0, ymax, by = 200)) +
    theme_bw(base_size = 16) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.8, 0.25),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.margin = margin(0, 1, 0, 0, "lines"))
  
  return(p)
  
  
}

unique_degs_plot(up)
ggsave("plots/Fig_5A_upreg_uniqueness.png",
       units = "in", dpi = 600,
       height = 3, width = 5.5)

unique_degs_plot(down)
ggsave("plots/Fig_5B_downreg_uniqueness.png",
       units = "in", dpi = 600,
       height = 3, width = 5.5)