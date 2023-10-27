library(Seurat)
library(CellChat)
library(ComplexHeatmap)
library(rstudioapi)
library(patchwork)
library(RColorBrewer)
library(scales)

obj_list <- readRDS("data_objects/11_cellchat_list.RDS")
obj_list
names(obj_list) <- c("Saline", "LPS")

cellchat <- mergeCellChat(obj_list, add.names = names(obj_list))
cellchat

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),
                           color.use = c("navy", "orange"), width = 1) + 
  scale_y_continuous(expand = c(0, 0))  + 
  theme(
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), 
                           measure = "weight",
                           color.use = c("navy", "orange"), width = 1) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )
gg1$layers[[2]] <- NULL
gg2$layers[[2]] <- NULL
gg1 + gg2
ggsave("plots/Fig_S11A-B_total_interactions_comparison.png",
       units = "in", dpi = 600,
       height = 4, width = 4.5)

#Show overall change in flow of information
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T, do.flip = F,
        color.use = c("navy", "orange"),
        measure = "count",
        font.size = 16,
        show.raw = F) + 
  scale_y_continuous(expand = c(0, 0)) + 
  ggtitle("Relative number of interactions within CellChat signaling classes") +
  labs(y = "Relative number\nof interactions") +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.text = element_text(size = 16),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 18)
  )
ggsave("plots/Fig_S11C_info_flow.png",
       units = "in", dpi = 600,
       width = 15, height = 4)

#Get the maximum edge weight in the list of objects
weight.max <- getMaxWeight(obj_list, attribute = c("idents", "count"))
#Show number of interactions, colored by sender cell, grouped by treatment group
#This makes the classic net bubble diagram with colored lines between sender-receiver pairs
par(mai = c(0.1, 1, 0.1, 1))
netVisual_circle(obj_list[[1]]@net$count, weight.scale = T, label.edge= F, 
                 edge.weight.max = weight.max[2], edge.width.max = 12, 
                 vertex.label.cex = 1.5)
savePlotAsImage(file = "plots/Fig_4B_chords_saline.png",
                format = "png",
                width = 800,
                height = 400)

par(mai = c(0.1, 1, 0.1, 1))
netVisual_circle(obj_list[[2]]@net$count, weight.scale = T, label.edge= F, 
                 edge.weight.max = weight.max[2], edge.width.max = 12, 
                 vertex.label.cex = 1.5)
savePlotAsImage(file = "plots/Fig_4B_chords_lps.png",
                format = "png",
                width = 800,
                height = 400)

#Heatmap of differential number of interactions for all cell pairs
#This provides more detail for the scatter plot above
netVisual_heatmap(cellchat,
                  font.size = 12,
                  font.size.title = 16,
                  color.heatmap = c("midnightblue", "orange3"),
                  width = 3)
savePlotAsImage(file = "plots/Fig_4C_change_in_signaling_heatmap.png",
                format = "png",
                width = 500, height = 500)

#This makes heatmaps of the relative strength of INCOMING signals from all cell types,
#classified into signaling systems
i <- 1
ht1 = netAnalysis_signalingRole_heatmap(obj_list[[i]], pattern = "incoming", 
                                        signaling = pathway.union,  
                                        width = 10, height = 25, 
                                        color.heatmap = "YlOrRd",
                                        font.size = 16)
ht2 = netAnalysis_signalingRole_heatmap(obj_list[[i+1]], pattern = "incoming", 
                                        signaling = pathway.union, 
                                        width = 10, height = 25, 
                                        color.heatmap = "YlOrRd",
                                        font.size = 16)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
savePlotAsImage(file = "plots/Fig_S11D_incoming_signals_by_treatment.png",
                format = "png",
                width = 1100,
                height = 1200)

#This makes heatmaps of the relative strength of OUTGOING signals from all cell types,
#classified into signaling systems
i <- 1
pathway.union <- union(obj_list[[i]]@netP$pathways, obj_list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(obj_list[[i]], pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        width = 10, 
                                        height = 25,
                                        font.size = 16,
                                        color.heatmap = "PuBuGn")
ht2 = netAnalysis_signalingRole_heatmap(obj_list[[i+1]], pattern = "outgoing", 
                                        signaling = pathway.union, 
                                        width = 10, 
                                        height = 25,
                                        font.size = 16,
                                        color.heatmap = "PuBuGn")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
savePlotAsImage(file = "plots/Fig_S11E_outgoing_signals_by_treatment.png",
                format = "png",
                width = 1100,
                height = 1200)

#LEGENDS:----
#For differential interactions heatmap Fig. 4C
ggplot(mtcars, aes(cyl, mpg, color = NA_real_)) + 
  geom_point() + 
  scale_color_gradientn(colours = c("midnightblue", "white", "orange3"),
                        limits = c(-10, 15),
                        values = rescale(c(-10, 0, 15)),
                        breaks = c(-10, -5, 0, 5, 10, 15),
                        name = "Relative values") + 
  guides(colour = guide_colorbar(title.position = "left")) +
  theme(legend.title = element_text(angle = 90,
                                    hjust = 0.5),
        legend.text = element_text(size = 12))
ggsave(file = "plots/Fig_4C_legend.png",
       units = "in", dpi = 600,
       width = 2, height = 2)

#For incoming signals heatmap SFig 11D
ggplot(mtcars, aes(cyl, mpg, color = NA_real_)) + 
  geom_point() + 
  scale_color_gradientn(colours = c("white", brewer.pal(name = "YlOrRd", n = 9)),
                        limits = c(0, 1),
                        breaks = c(0, 1),
                        name = "Relative strength") + 
  guides(colour = guide_colorbar(title.position = "left")) +
  theme(legend.title = element_text(angle = 90,
                                    hjust = 0.5),
        legend.text = element_text(size = 12))
ggsave(file = "plots/Fig_S11D_legend.png",
       units = "in", dpi = 600,
       width = 2, height = 2)

#For outgoing signals heatmap SFig 11E
ggplot(mtcars, aes(cyl, mpg, color = NA_real_)) + 
  geom_point() + 
  scale_color_gradientn(colours = c("white", brewer.pal(name = "PuBuGn", n = 9)),
                        limits = c(0, 1),
                        breaks = c(0, 1),
                        name = "Relative strength") + 
  guides(colour = guide_colorbar(title.position = "left")) +
  theme(legend.title = element_text(angle = 90,
                                    hjust = 0.5),
        legend.text = element_text(size = 12))
ggsave(file = "plots/Fig_S11E_legend.png",
       units = "in", dpi = 600,
       width = 2, height = 2)
