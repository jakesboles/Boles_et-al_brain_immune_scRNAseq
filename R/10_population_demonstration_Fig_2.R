library(Seurat)
library(tidyverse)
library(scCustomize)
library(janitor)
library(ggbeeswarm)
library(dittoSeq)
library(rstatix)
library(ggpubr)
library(camcorder)
library(ggh4x)

obj <- readRDS("data_objects/07_annotated_object.RDS")

plots <- "plots/Fig_2"

obj@active.ident <- factor(obj@active.ident,
                           levels = sort(levels(obj))) 

theme <- theme(axis.text = element_blank(),
               axis.ticks = element_blank(),
               axis.line = element_line(arrow = arrow(angle = 15, 
                                                      length = unit(0.5, "cm"), type = "closed")))

DimPlot_scCustom(obj, 
                 label = F,
                 colors_use = DiscretePalette_scCustomize(num_colors = 20,
                                                          palette = "varibow"),
                 figure_plot = F) +
  theme

ggsave(paste0(plots, "B.png"),
       units = "in", dpi = 600,
       height = 7, width = 8)

#Highlight rare cell types that we've captured ----
gd <- WhichCells(obj, idents = "γδ T-cells")
mast <- WhichCells(obj, idents = "Mast cells")
baso <- WhichCells(obj, idents = "Basophils")
plasma <- WhichCells(obj, idents = "Plasma cells")

pal <- DiscretePalette_scCustomize(num_colors = 20, palette = "varibow")

highlight <- function(cells, color, title, panel){
  l1 <- length(cells)
  l2 <- nrow(obj@meta.data)
  
  frac <- (l1 / l2) * 100
  frac <- round(frac, digits = 3)
  
  
  p <- Cell_Highlight_Plot(obj,
                           cells_highlight = list(title = cells),
                           highlight_color = color) +
    ggtitle(title,
            subtitle = paste0(frac, "% of total cells")) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    NoLegend() + 
    theme
  
  file <- deparse(substitute(cells))
  
  ggsave(p,
         filename = paste0(plots, panel, ".png"),
         units = "in", dpi = 600,
         height = 4, width = 4)
}

highlight(gd, pal[20], "γδ T-cells", "F")
highlight(mast, pal[10], "Mast cells", "D")
highlight(baso, pal[3], "Basophils", "C")
highlight(plasma, pal[17], "Plasma cells", "E")

#Donut plot ----
df <- obj@meta.data

donut <- df %>%
  mutate(total = n()) %>%
  group_by(cell_type) %>%
  mutate(fraction = (n() / total) * 100,
         count = n()) %>%
  distinct(cell_type, .keep_all = T) %>%
  dplyr::select(cell_type, total, fraction, count) %>%
  mutate(cell_type = factor(cell_type,
                            levels = sort(cell_type))) %>%
  arrange(cell_type)

donut$ymax <- cumsum(donut$fraction)
donut$ymin <- c(0, head(donut$ymax, n=-1))

pal2 <- DiscretePalette_scCustomize(num_colors = 20,
                                    palette = "polychrome")

donut %>%
  mutate(cell_type_2 = paste0(count, cell_type, sep = " ")) %>%
  ggplot(aes(ymax = ymax, ymin = ymin,
             xmax = 21, xmin = 20,
             fill = cell_type)) +
  geom_rect() + 
  coord_polar(theta = "y") + 
  scale_fill_manual(values = pal2, labels = paste(donut$count, donut$cell_type)) +
  theme_void(base_size = 16) + 
  ggtitle("Cell frequency") + 
  labs(fill = paste0(donut$total, " total cells")) +
  xlim(c(19.5, NA)) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, vjust = -5)
  )
ggsave(paste0(plots, "G.png"),
       units = "in", dpi = 600,
       height = 8, width = 8)

obj@meta.data <- obj@meta.data %>%
  mutate(immune = if_else(cell_type %in% c("Astrocytes", "Choroid plexus cells",
                                           "Endothelial cells", "Erythrocytes",
                                           "Leptomeningeal cells",
                                           "Oligodendrocytes", "Pericytes", "Stem cells"),
                          "Non-immune cells", "Immune cells"))
dittoBarPlot(obj,
             var = "immune", 
             group.by = "orig.ident",
             color.panel = c("dodgerblue", "grey30"),
             main = "Immune cell purity according to RNA",
             xlab = "Sample",
             ylab = "Fraction of total sample",
             scale = "percent") + 
  scale_x_discrete(limits = rev) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  coord_flip() +
  theme_classic(base_size = 16) +
  theme(
    axis.text = element_text(color = "black"),
    #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
ggsave(paste0(plots, "H.png"),
       units = "in", dpi = 600,
       height = 4, width = 8)
