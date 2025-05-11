# Visualization of scRNA-seq results using Seurat ----
# Author: Ana Añazco | Companion script to reproduce thesis figures

# Load libraries ----
library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)
library(Nebulosa)
library(gridExtra)
library(readxl)
library(openxlsx)
library(ggpubr)
library(writexl)

# Load Seurat object ----
# Replace with your actual object path
seurat_obj <- readRDS("seurat_object.rds")

# Set identity column for plotting
Idents(seurat_obj) <- "clusters"

# Feature plots ----
plot_density(seurat_obj, "Pten")
plot_density(seurat_obj, "Pten") + facet_grid(. ~ seurat_obj$sample)

# Violin plots ----
VlnPlot(seurat_obj, features = c("Pcif1"), assay = "SCT", log = TRUE)

# Custom ggplot violin
data <- FetchData(seurat_obj, vars = c("Pcif1", "ident"))
ggplot(data, aes(x = ident, y = Pcif1, fill = ident)) + 
  geom_violin() + 
  geom_jitter(size = 0.3, alpha = 1) + 
  theme_minimal() +
  labs(y = "Expression Level", x = "Identity")

# UMAPs by cluster and sample ----
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'clusters')
DimPlot(seurat_obj, reduction = "umap", repel = TRUE, group.by = 'sample')

# RMP DotPlot (RMP_annotated.xlsx is provided) ----
markers <- read_xlsx("RMP_annotated.xlsx", sheet = 1)
genes <- markers$Genes
DotPlot(seurat_obj, features = genes, cols = "RdBu", group.by = "clusters",
        assay = "SCT", dot.scale = 6) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# RMP DotPlot, filter by pct ----
dotdata <- DotPlot(seurat_obj, features = genes)$data
filtered_genes <- dotdata %>% filter(pct.exp > 35) %>% distinct(features.plot)
genes_ordered <- sort(as.character(filtered_genes$features.plot))
DotPlot(seurat_obj, features = genes_ordered, cols = "RdBu", group.by = "clusters",
        dot.scale = 4.5, scale = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# Export DotPlot data
write.xlsx(dotdata, "DotPlot_data.xlsx")

# Stackbar plot ----
# This data comes from the output of script 01 (cluster_frequencies.xlsx)
frequencies <- read_xlsx("cluster_frequencies.xlsx")
frequencies$Cell_type <- factor(frequencies$Cluster, levels = c("Myeloid", "DC", 
                                                                "B cells", "T_cells", 
                                                                "Stromal", "Endothelial", 
                                                                "Luminal", "Basal"))

Set3 <- c("#5f91a8", "#B2C69AFF", "#E6BBC7FF", "#ffd182", "#ee8bd5", "#69e287", 
          "#bb9de7", "#83afa6" ,"#f3aaa0", "#76e5f8", "#c79666", "#e75555")
ggplot(frequencies, aes(fill = Cell_type, y = Percentage, x = Sample)) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  scale_fill_brewer(palette = "Set3")

# Marker analysis in one cluster----
Idents(seurat_obj) <- "clusters"
luminal_markers <- FindMarkers(seurat_obj, ident.1 = "Luminal")
write.xlsx(luminal_markers, file = "FindMarkers_Luminal.xlsx")

# FindAllMarkers and export by cluster ----
all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
wb <- createWorkbook()
for (cl in unique(all_markers$cluster)) {
  addWorksheet(wb, sheetName = as.character(cl))
  writeData(wb, sheet = as.character(cl), x = all_markers[all_markers$cluster == cl, ])
}
saveWorkbook(wb, "FindAllMarkers_by_cluster.xlsx", overwrite = TRUE)

# Top markers heatmap ----
top10 <- all_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()
write_xlsx(top10, "Top10_harmony_annot.xlsx")
