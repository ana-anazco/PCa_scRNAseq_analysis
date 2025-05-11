# GO term analysis functions ----
# Author: Ana Añazco | Modular functions to perform GO enrichment from Seurat or Excel

library(clusterProfiler)
library(org.Mm.eg.db) # Use org.Hs.eg.db if human
library(readxl)
library(writexl)
library(ggplot2)
library(enrichplot)

# Function 1: From Excel with FindMarkers output format ----
run_GO_from_excel <- function(excel_path, sheet, direction = c("up", "down"), plot_type = c("dotplot", "barplot", "goplot", "emapplot"), output_dir = "./GO_results/") {
  dir.create(output_dir, showWarnings = FALSE)
  
  df <- read_xlsx(excel_path, sheet = sheet)
  if (!"avg_log2FC" %in% colnames(df)) {
    stop("Column 'avg_log2FC' not found. Make sure your Excel is in FindMarkers format.")
  }
  
  if ("gene" %in% colnames(df)) {
    df$genes <- df$gene
  } else {
    df$genes <- rownames(df)
  }
  FC <- data.frame(log2FC = df$avg_log2FC, genes = df$genes)
  
  if (direction == "up") {
    selected_genes <- FC$genes[FC$log2FC > 1]
  } else {
    selected_genes <- FC$genes[FC$log2FC < -1]
  }
  
  genes_mapped <- bitr(selected_genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db)
  
  enrich_result <- enrichGO(gene = genes_mapped$ENTREZID,
                            OrgDb = org.Mm.eg.db,
                            ont = "BP",
                            pvalueCutoff = 0.05)
  
  plot_type <- match.arg(plot_type)
  sheet_name <- gsub(" ", "_", sheet)
  plot_file <- paste0(output_dir, plot_type, "_", sheet_name, ".pdf")
  
  pdf(plot_file)
  if (plot_type == "dotplot") print(dotplot(enrich_result, showCategory = 20, font.size = 7))
  if (plot_type == "barplot") print(barplot(enrich_result, showCategory = 20))
  if (plot_type == "goplot") print(goplot(enrich_result))
  if (plot_type == "emapplot") {
    enrich_result2 <- pairwise_termsim(enrich_result)
    print(emapplot(enrich_result2, showCategory = 25))
  }
  dev.off()
  
  return(enrich_result)
}

# Function 2: From SeuratObject using FindMarkers directly ----
run_GO_from_seurat <- function(seurat_obj, idents_column, ident1, direction = c("up", "down"), plot_type = c("dotplot", "barplot", "goplot", "emapplot"), output_name = "GO_plot", output_dir = "./GO_results/") {
  dir.create(output_dir, showWarnings = FALSE)
  
  Idents(seurat_obj) <- idents_column
  markers <- FindMarkers(seurat_obj, ident.1 = ident1, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  
  markers$genes <- rownames(markers)
  FC <- data.frame(log2FC = markers$avg_log2FC, genes = markers$genes)
  
  if (direction == "up") {
    selected_genes <- FC$genes[FC$log2FC > 1]
  } else {
    selected_genes <- FC$genes[FC$log2FC < -1]
  }
  
  genes_mapped <- bitr(selected_genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db)
  
  enrich_result <- enrichGO(gene = genes_mapped$ENTREZID,
                            OrgDb = org.Mm.eg.db,
                            ont = "BP",
                            pvalueCutoff = 0.05)
  
  plot_type <- match.arg(plot_type)
  plot_file <- paste0(output_dir, plot_type, "_", output_name, ".pdf")
  
  pdf(plot_file)
  if (plot_type == "dotplot") print(dotplot(enrich_result, showCategory = 20, font.size = 7))
  if (plot_type == "barplot") print(barplot(enrich_result, showCategory = 20))
  if (plot_type == "goplot") print(goplot(enrich_result))
  if (plot_type == "emapplot") {
    enrich_result2 <- pairwise_termsim(enrich_result)
    print(emapplot(enrich_result2, showCategory = 25))
  }
  dev.off()
  
  return(enrich_result)
}

# Examples of usage ----

# Example 1: Using an Excel file with a sheet named 'Basal'
# result_excel <- run_GO_from_excel(
#   excel_path = "FindMarkers_SeuratObject.xlsx",
#   sheet = "Basal",
#   direction = "up",
#   plot_type = "dotplot",
#   output_dir = "./GO_results/"
# )

# Example 2: Using a Seurat object directly
# result_seurat <- run_GO_from_seurat(
#   seurat_obj = harmony_no_SV,
#   idents_column = "clusters",
#   ident1 = "Tumour_Cast",
#   direction = "down",
#   plot_type = "emapplot",
#   output_name = "Tumour_Cast_down",
#   output_dir = "./GO_results/"
# )
