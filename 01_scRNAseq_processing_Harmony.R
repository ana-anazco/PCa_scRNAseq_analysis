# scRNA-seq analysis pipeline using Seurat and Harmony
# Author: Ana Añazco

# Load libraries ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(harmony)
library(writexl)
library(Nebulosa)
library(readxl)

# Define input/output directories ----
base_dir <- "./Matrix/" # Change to your data path
output_dir <- "./results/"
dir.create(output_dir, showWarnings = FALSE)

# Load 10X data automatically from folders ----
dirs <- list.dirs(path = base_dir, recursive = FALSE, full.names = TRUE)
sample_names <- basename(dirs)

# Function to create Seurat object from 10X folder ----
create_seurat_object <- function(dir, sample_name) {
  counts <- Read10X(data.dir = dir)
  CreateSeuratObject(counts = counts, project = "PCproject", 
                     min.cells = 3, min.features = 200) %>% 
    RenameCells(add.cell.id = sample_name)
}

# Create named list and merge Seurat objects ----
seurat_list <- mapply(create_seurat_object, dir = dirs, sample_name = sample_names, SIMPLIFY = FALSE)
names(seurat_list) <- sample_names
merged <- Reduce(function(x, y) merge(x, y), seurat_list)
merged <- Reduce(function(x, y) merge(x, y), seurat_list)

# QC & filtering ----
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^mt-")
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Adaptive thresholds ----
minCov <- 1000
countLOW <- if (min(merged$nCount_RNA) >= minCov) min(merged$nCount_RNA) else quantile(merged$nCount_RNA, 0.01)
countHIGH <- quantile(merged$nCount_RNA, 0.99)
featureLOW <- quantile(merged$nFeature_RNA, 0.01)

# Subset based on thresholds ----
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & 
                   percent.mt < 5 & nCount_RNA > countLOW & nCount_RNA < countHIGH)

# Sample annotations ----
# NOTE: the following mapping is specific to this dataset. Modify if sample folder names differ.
merged$sample <- sapply(str_split_fixed(Cells(merged), "_", 2)[,1], function(x) {
  recode(x,
         "S01" = "Cast_KO_1", "S02" = "Cast_WT_1",
         "S03" = "NC_KO_1", "S04" = "Ctrl_1",
         "S05" = "Ctrl_2", "S06" = "NC_KO_2",
         "S07" = "Cast_KO_2", "S08" = "Cast_WT_2")
})

# Collapse biological replicates ----
# This step groups technical replicates under a common biological condition. Adapt as needed.
merged$sample_group <- recode(merged$sample,
                              Ctrl_1 = "Healthy", Ctrl_2 = "Healthy",
                              NC_KO_1 = "Tumour", NC_KO_2 = "Tumour",
                              Cast_WT_1 = "Healthy_Cast", Cast_WT_2 = "Healthy_Cast",
                              Cast_KO_1 = "Tumour_Cast", Cast_KO_2 = "Tumour_Cast")

# SCTransform & Harmony integration ----
merged <- SCTransform(merged, verbose = FALSE) %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunHarmony(group.by.vars = "sample", reduction.save = "harmony") %>%
  RunUMAP(reduction = "harmony", dims = 1:40) %>%
  FindNeighbors(reduction = "harmony", dims = 1:40) %>%
  FindClusters(resolution = 1.2)

# Save intermediate object ----
saveRDS(merged, file = file.path(output_dir, "harmony_integrated.rds"))

# Top variable features & PCA plot ----
VariableFeaturePlot(merged) + LabelPoints(points = head(VariableFeatures(merged), 10))

# Marker genes ----
markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_xlsx(markers, file.path(output_dir, "FindAllMarkers.xlsx"))

# Heatmap of top markers ----
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(merged, features = top10$gene) + NoLegend()
pdf(file.path(output_dir, "heatmap_top10.pdf"))
DoHeatmap(merged, features = top10$gene) + NoLegend()
dev.off()

# Cell type annotation (optional) ----
db_file <- "ScTypeDB_short.xlsx"  # update if needed
tissue <- "Prostate"
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

gs_list <- gene_sets_prepare(db_file, tissue)
es.max <- sctype_score(scRNAseqData = merged[["SCT"]]@scale.data, scaled = TRUE,
                       gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# Assign cluster labels ----
sctype_results <- do.call("rbind", lapply(unique(merged$seurat_clusters), function(cl){
  es.cl <- sort(rowSums(es.max[, WhichCells(merged, idents = cl)]), decreasing = TRUE)
  head(data.frame(cluster = cl, type = names(es.cl), score = es.cl, 
                  ncells = length(WhichCells(merged, idents = cl))), 1)
}))
sctype_results$type[sctype_results$score < sctype_results$ncells / 4] <- "Unknown"

# Save annotation and update metadata ----
write_xlsx(sctype_results, file.path(output_dir, "annotated_clusters.xlsx"))
merged$celltype <- sctype_results$type[match(merged$seurat_clusters, sctype_results$cluster)]

# Remove Seminal Vesicle cluster (SV) ----
merged <- subset(merged, idents = "SV", invert = TRUE)
saveRDS(merged, file = file.path(output_dir, "harmony_no_SV.rds"))

# Cell counts per cluster/sample ----
count_df <- table(merged$celltype, merged$sample_group) %>% as.data.frame()
colnames(count_df) <- c("Cluster", "Sample", "Count")
count_df <- count_df %>%
  group_by(Sample) %>%
  mutate(Percentage = Count / sum(Count) * 100)
write_xlsx(count_df, file.path(output_dir, "cluster_frequencies.xlsx"))

message("scRNA-seq analysis complete. Objects saved in results/")

# Optional: re-subclustering specific populations ----
## Example: recluster 'Basal' cells (uncomment to use)

## Check graph names to identify available graph
# names(merged@graphs)

## If no graph exists, build neighbor graph first
# merged <- FindNeighbors(merged, dims = 1:10, graph.name = "SCT_snn")

## Perform subclustering on the 'Basal' cluster
# merged <- FindSubCluster(
#   object = merged,
#   cluster = "Basal",
#   graph.name = "SCT_snn",
#   subcluster.name = "Basal_",
#   resolution = 0.05,
#   algorithm = 1
# )

