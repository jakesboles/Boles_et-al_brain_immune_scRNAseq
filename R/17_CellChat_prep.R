library(CellChat)
library(Seurat)

obj <- readRDS("data_objects/09_immune_object.RDS")

sal <- subset(obj, subset = Treatment == "Saline")
lps <- subset(obj, subset = Treatment == "LPS")

objs <- list(Saline = sal, 
             LPS = lps)
cellchat_list <- list()

for (i in seq_along(objs)){
  labels <- Idents(objs[[i]])
  
  input <- GetAssayData(objs[[i]], assay = "SCT", slot = "data")
  
  meta <- objs[[i]]@meta.data
  meta <- meta[, c(1:14, 24, 25)]
  
  message(paste0("---------------------------------------Creating ", names(objs)[i], " CellChat object---------------------------------------"))
  
  cellchat_list[[i]] <- createCellChat(input, meta, group.by = "cell_type")
  
  message(paste0("---------------------------------------Loading DB and adding to ", names(objs)[i], " CellChat object---------------------------------------"))
  
  cellchatdb <- CellChatDB.mouse
  cellchat_list[[i]]@DB <- cellchatdb
  
  cellchat_list[[i]] <- subsetData(cellchat_list[[i]])
  
  message(paste0("---------------------------------------Identifying over-expressed genes in ", names(objs)[i], " CellChat object---------------------------------------"))
  
  cellchat_list[[i]] <- identifyOverExpressedGenes(cellchat_list[[i]])
  
  message(paste0("---------------------------------------Identifying over-expressed interactions in ", names(objs)[i], " CellChat object---------------------------------------"))
  
  cellchat_list[[i]] <- identifyOverExpressedInteractions(cellchat_list[[i]])
  cellchat_list[[i]]
  
  message(paste0("---------------------------------------Calculating community probability in ", names(objs)[i], " CellChat object---------------------------------------"))
  
  cellchat_list[[i]] <- computeCommunProb(cellchat_list[[i]],
                                          population.size = T)
  message(paste0("---------------------------------------Calculating pathway probability in ", names(objs)[i], " CellChat object---------------------------------------"))
  
  cellchat_list[[i]] <- computeCommunProbPathway(cellchat_list[[i]])
  
  message(paste0("---------------------------------------Aggregating network in ", names(objs)[i], " CellChat object---------------------------------------"))
  
  cellchat_list[[i]] <- aggregateNet(cellchat_list[[i]])
  
  message(paste0("---------------------------------------Computing network centrality in ", names(objs)[i], " CellChat object---------------------------------------"))
  
  cellchat_list[[i]] <- netAnalysis_computeCentrality(cellchat_list[[i]],
                                                      slot.name = "netP")
  
  message(paste0("---------------------------------------Done with ", names(objs)[i], "!---------------------------------------"))
}

saveRDS(cellchat_list,
        file = "data_objects/11_cellchat_list.RDS")