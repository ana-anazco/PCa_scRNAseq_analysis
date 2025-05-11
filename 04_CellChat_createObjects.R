# CellChat analysis per sample ----
# Author: Ana Añazco 

library(CellChat)
library(Seurat)
library(patchwork)
library(dplyr)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE)

# Function to run CellChat on multiple samples within a Seurat object ----
run_cellchat_by_sample <- function(seurat_obj, sample_column = "sample", cluster_column = "clusters", assay = "SCT", output_prefix = "CellChat") {
  samples <- unique(seurat_obj[[sample_column]][, 1])
  
  for (s in samples) {
    obj_sample <- subset(seurat_obj, subset = get(sample_column) == s)
    
    # Create CellChat object
    cellchat <- createCellChat(object = obj_sample, group.by = cluster_column, assay = assay)
    cellchat@meta$samples <- cellchat@meta[[sample_column]]
    cellchat@DB <- subsetDB(CellChatDB.mouse)
    
    # Run CellChat pipeline
    cellchat <- subsetData(cellchat)
    future::plan("multisession", workers = 4)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    # Save CellChat object
    saveRDS(cellchat, file = paste0(output_prefix, "_", s, ".rds"))
  }
  
  message("CellChat analysis completed for all samples.")
}

# Optional: set custom order for cell types in plots (if needed) ----
# new_levels <- c("Luminal", "Basal", "PrU", "Stromal", "Endothelial", "DC", "Myeloid", "T_cells", "B cells")
# cellchat@meta$ident <- factor(cellchat@meta$ident, levels = new_levels)

# Example usage ----
# run_cellchat_by_sample(seurat_obj = Harmony, sample_column = "sample", cluster_column = "clusters")
