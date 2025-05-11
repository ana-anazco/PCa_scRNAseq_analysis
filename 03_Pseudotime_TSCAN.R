# Pseudotime analysis with TSCAN ----
# Author: Ana Añazco | Balanced sampling and trajectory inference

library(Seurat)
library(TSCAN)
library(scuttle)
library(scater)
library(scran)
library(ggplot2)
library(ggraph)

# Function to balance cell numbers per cluster by sample ----
subset_balanced_cells <- function(seurat_obj, ident_column = "seurat_clusters", max_cells = 250) {
  samples <- unique(seurat_obj$sample)
  balanced_list <- list()
  
  for (s in samples) {
    obj_sample <- subset(seurat_obj, subset = sample == s)
    Idents(obj_sample) <- ident_column  # <- ESTA LÍNEA ES CLAVE
    
    list_subset <- list()
    for (cl in unique(Idents(obj_sample))) {
      idx <- WhichCells(obj_sample, idents = cl)
      if (length(idx) > max_cells) {
        idx <- sample(idx, max_cells)
      }
      list_subset[[as.character(cl)]] <- idx
    }
    
    sample_name <- paste0("balanced_", s)
    balanced_obj <- obj_sample[, unlist(list_subset)]
    balanced_obj <- SetIdent(balanced_obj, value = ident_column)
    saveRDS(balanced_obj, paste0(sample_name, ".rds"))
  }
  
  message("Balanced Seurat objects saved per sample.")
}


# Function to perform TSCAN pseudotime and plot trajectory ----
run_tscan_pseudotime <- function(rds_path, cluster_column = "seurat_clusters", start_cluster, output_prefix = "TSCAN") {
  obj <- readRDS(rds_path)
  sce <- as.SingleCellExperiment(obj, assay = "SCT")
  
  clusters <- as.factor(colData(sce)[[cluster_column]])
  aggregated <- aggregateAcrossCells(sce, ids = clusters, use.assay.type = "logcounts")
  centroids <- reducedDim(aggregated, "PCA")
  
  mst <- TSCAN::createClusterMST(centroids, clusters = NULL)
  edge_data <- reportEdges(aggregated, mst = mst, clusters = NULL, use.dimred = "UMAP")
  
  colLabels(sce) <- clusters
  mapped <- mapCellsToEdges(sce, mst = mst, use.dimred = "PCA")
  ordered <- orderCells(mapped, mst, start = as.character(start_cluster))
  pseudo_avg <- averagePseudotime(pathStat(ordered))
  
  plot <- plotUMAP(sce, colour_by = I(pseudo_avg), text_by = cluster_column, text_colour = "black") +
    geom_line(data = edge_data, mapping = aes(x = UMAP_1, y = UMAP_2, group = edge))
  
  pdf(paste0(output_prefix, "_", start_cluster, ".pdf"))
  print(plot)
  dev.off()
  
  message("Pseudotime trajectory saved: ", output_prefix, "_", start_cluster, ".pdf")
  return(invisible(list(pseudotime = pseudo_avg, mst = mst)))
}

# Function to plot gene expression along MST trajectory ----
Gene_trajectory <- function(SeuratObject, gene_list, cluster_column = "sub_", output_prefix = "Trajectory") {
  set.seed(5054)
  
  sce <- as.SingleCellExperiment(SeuratObject)
  clusters_factor <- as.factor(sce@colData[[cluster_column]])
  by.cluster <- scuttle::aggregateAcrossCells(sce, ids = clusters_factor, use.assay.type = "logcounts")
  centroids <- reducedDim(by.cluster, "PCA")
  mst <- TSCAN::createClusterMST(centroids, clusters = NULL)
  line.data <- TSCAN::reportEdges(by.cluster, mst = mst, clusters = NULL, use.dimred = "UMAP")
  
  colLabels(sce) <- clusters_factor
  pseudo <- TSCAN::quickPseudotime(sce, use.dimred = "PCA", outgroup = TRUE)
  
  plot1 <- DotPlot(SeuratObject, features = gene_list, cols = "RdBu", group.by = cluster_column, assay = "SCT")
  
  for (g in unique(plot1$data$features.plot)) {
    gene_data <- subset(plot1$data, features.plot == g)
    gene_data <- gene_data[order(gene_data$id), ]
    pct_exp <- gene_data$pct.exp
    avg_exp <- gene_data$avg.exp.scaled
    
    temp_plot <- ggraph(pseudo$mst, layout = "igraph", algorithm = 'kk') +
      geom_edge_fan(show.legend = FALSE) +
      geom_node_point(aes(fill = avg_exp, size = pct_exp), alpha = 1, shape = 21, color = "black") +
      geom_node_text(aes(label = name), nudge_y = -0.1, nudge_x = -0.1) +
      scale_fill_gradient2(name = "Avg. Expression", limits = c(-2, 2), low = "#0066ff", high = "red") +
      scale_size_continuous(name = "% Expressed", limits = c(0, 100), range = c(2, 8), breaks = c(0, 25, 50, 75, 100)) +
      labs(title = g) +
      theme_void()
    
    ggsave(filename = paste0(output_prefix, "_", g, ".pdf"), plot = temp_plot, width = 7, height = 5, units = 'in')
  }
}

# Example usage ----

# Step 1: Create balanced objects per sample
# subset_balanced_cells(seurat_obj = LumBasal, ident_column = "sub_", max_cells = 250)

# Step 2: Run TSCAN pseudotime from one of the balanced samples
# run_tscan_pseudotime(rds_path = "balanced_Healthy_Cast.rds", cluster_column = "sub_", start_cluster = 8, output_prefix = "TSCAN_HC")

# Step 3: Plot gene expression on MST trajectory
# genes <- c("Krt8", "Krt14", "Krt15")
# Gene_trajectory(SeuratObject = balanced_Healthy, gene_list = genes, cluster_column = "sub_", output_prefix = "Trajectory_H")

