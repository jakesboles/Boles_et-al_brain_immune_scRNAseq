library(Seurat)
library(SingleCellExperiment)
library(multinichenetr)
library(ggplot2)
library(tidyverse)
library(scCustomize)
library(ggpubr)
library(cowplot)
library(grid)
library(rstudioapi)
library(camcorder)
library(paletteer)

obj <- readRDS("data_objects/11_multinichenetr.RDS")

contrast_tbl <- tibble(
  contrast = c("LPS-Saline","Saline-LPS"),
  group = c("LPS", "Saline"))

ligand_target_matrix <- readRDS("04_data_objects/07_multinichenetr/ligand_target_matrix_nsga2r_final_mouse.rds")
colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
  make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
  make.names()

#Microglia to cells that are changing with LPS
prioritized_tbl_oi_all <- get_top_n_lr_pairs(obj$prioritization_tables, 75, 
                                             rank_per_group = F,
                                             senders_oi = "Microglia",
                                             receivers_oi = c("B.cells", "Dendritic.cells",
                                                              "Monocytes", "Neutrophils"))

prioritized_tbl_oi <- obj$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>%
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% 
                            unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

plot_oi <- make_sample_lr_prod_activity_plots(obj$prioritization_tables, prioritized_tbl_oi_all)

gg_record(device = "png",
          units = "in", dpi = 600,
          height = 10, width = 5)
p <- plot_oi[[1]] + 
  scale_y_discrete(position = "right",
                   limits = rev) + 
  guides(size = F) +
  labs(color = "Scaled ligand-receptor\npseudobulked expression") + 
  theme(axis.text.x = element_blank(),
        legend.position = "bottom")
p
p$data %>%
  ggplot(aes(x = sample, y = lr_interaction)) + 
  geom_tile(aes(fill = scaled_LR_pb_prod),
            color = "white",
            linejoin = "bevel") +
  scale_fill_gradient2(low = "midnightblue",
                       high = "orange3") +
  scale_y_discrete(limits = rev,
                   position = "right",
                   expand = c(0,0)) +
  facet_grid(rows = vars(sender_receiver),
             cols = vars(factor(group, levels = c("Saline", "LPS"))),
             switch = "y",
             scales = "free",
             space = "free") + 
  scale_x_discrete(expand = c(0, 0)) +
  labs(fill = "Scaled ligand-receptor\npseudobulked expression") + 
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold.italic", color = "black", size = 9, hjust = 0),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(color = "black", face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "grey15"),
        axis.title = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        strip.text.x.top = element_text(size = 11),
        legend.position = "bottom",
        panel.spacing.x = unit(1, "pt"),
        legend.text = element_text(size = 9),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, 0, -10))
ggsave("plots/Fig_4D.png",
       units = "in", dpi = 600,
       height = 10, width = 5)

#Microglia-microglia communication----
prioritized_tbl_oi_all <- get_top_n_lr_pairs(obj$prioritization_tables, 30, 
                                             rank_per_group = F,
                                             senders_oi = "Microglia",
                                             receivers_oi = "Microglia")

prioritized_tbl_oi <- obj$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>%
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% 
                            unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

plot_oi <- make_sample_lr_prod_activity_plots(obj$prioritization_tables, prioritized_tbl_oi_all)
plot_oi

gg_record(device = "png",
          units = "in", dpi = 600,
          height = 5, width = 5)
p <- plot_oi[[1]] + 
  scale_y_discrete(position = "right",
                   limits = rev) + 
  guides(size = F) +
  labs(color = "Scaled\nligand-receptor\npseudobulked\nexpression") + 
  theme(axis.text.x = element_blank(),
        legend.position = "right")
p$data %>%
  ggplot(aes(x = sample, y = lr_interaction)) + 
  geom_tile(aes(fill = scaled_LR_pb_prod),
            color = "white",
            # shape = 21,
            # width = 2,
            linejoin = "bevel") +
  scale_fill_gradient2(low = "midnightblue",
                       high = "orange3",
                       breaks = c(-1, 0, 1)) +
  scale_y_discrete(limits = rev,
                   position = "right",
                   expand = c(0,0)) +
  facet_grid(rows = vars(sender_receiver),
             cols = vars(factor(group, levels = c("Saline", "LPS"))),
             switch = "y",
             scales = "free",
             space = "free") + 
  scale_x_discrete(expand = c(0, 0)) +
  labs(fill = "Scaled ligand-receptor\npseudobulked expression") + 
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold.italic", color = "black", size = 9, hjust = 0),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(color = "black", face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "grey15"),
        axis.title = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        strip.text.x.top = element_text(size = 11),
        legend.position = "bottom",
        panel.spacing.x = unit(1, "pt"),
        legend.text = element_text(size = 9),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, 0, -10))
ggsave("plots/Fig_S12.png",
       units = "in", dpi = 600,
       height = 5, width = 5)
