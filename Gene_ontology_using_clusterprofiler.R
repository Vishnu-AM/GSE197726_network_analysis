# Load necessary libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(tidyverse)

# Define input and output paths
top_25_file <- "/home/vishnu/transcriptomics_for_thesis_with_shriyansh/network_analysis_for_v2_files/turquoise/degree/top_25_nodes_by_degree.csv"
output_folder <- "/home/vishnu/transcriptomics_for_thesis_with_shriyansh/GO_gse197726/final_GO/turquoise"

# Create output folder if it doesn't exist
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# Read gene symbols
top_genes <- read.csv(top_25_file, stringsAsFactors = FALSE,sep = "\t")
gene_symbols <- na.omit(top_genes[[2]])  # adjust index if column name is different

# Run enrichment analysis
ego <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff = 0.1,
  readable      = TRUE
)

  
  
# Skip if no enrichment
if (!is.null(ego) && nrow(ego) > 0) {
  
  # Save as TSV
  write.table(as.data.frame(ego),
              file = file.path(output_folder, "enrichGO_top25_turquoise.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  
  # Save dotplot - TIFF
  tiff(file.path(output_folder, "dotplot_top25.tiff"),
       width = 12, height = 15, units = "in", res = 300, bg = "transparent")
  print(dotplot(ego, showCategory =30))
  dev.off()
  
  # Save dotplot - SVG
  svg(file.path(output_folder, "dotplot_top25.svg"),
      width = 12, height = 16, bg = "transparent")
  print(dotplot(ego, showCategory = 30))
  dev.off()
  
  # Save cnetplot - TIFF with legend at top right
  tiff(file.path(output_folder, "cnetplot_top25.tiff"),
       width = 10, height = 10, units = "in", res = 300, bg = "transparent")
  
  print(
    cnetplot(ego, showCategory = 30, node_label = "all",
             cex_label_category = 1.5, cex_label_gene = 1, layout = "kk") +
      theme(
        legend.position = c(0.95, 0.5),
        legend.justification = c(1, 1),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      )
  )
  
  dev.off()
  
  # Save cnetplot - SVG with legend at top right
  svg(file.path(output_folder, "cnetplot_top25.svg"),
      width = 10, height = 10, bg = "transparent")
  
  print(
    cnetplot(ego, showCategory = 30, node_label = "all",
             cex_label_category = 1.5, cex_label_gene = 1, layout = "kk") +
     theme(
        legend.position = c(0.95, 0.5),
        legend.justification = c(1, 1),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      )
  )
  
  dev.off()
  
  
} else {
  cat("No enrichment results found.\n")
}


##########################COMMUNITY BASED GENE ONTOLOGY#########################################################

community_data <- read.csv("/home/vishnu/transcriptomics_for_thesis_with_shriyansh/network_analysis_for_v2_files/turquoise/community/louvain_communities.csv", stringsAsFactors = FALSE,sep = "\t")


# Pivot longer to get gene + community pairs
long_df <- pivot_longer(community_data, 
                        cols = everything(), 
                        names_to = "community", 
                        values_to = "gene_symbol") %>%
  filter(!is.na(gene_symbol) & gene_symbol != "")


gene_communities <- long_df %>%
  group_by(community) %>%
  summarise(genes = list(unique(gene_symbol))) %>%
  deframe()
str(gene_communities)

ck <- compareCluster(
  geneCluster = gene_communities,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  
  readable = TRUE
)

head(ck)

dotplot(ck, showCategory = 3) +
  ggtitle("GO Turquoise") 
 



write.table(ck, "/home/vishnu/transcriptomics_for_thesis_with_shriyansh/GO_gse197726/final_GO/turquoise/community_GO/GO_turquoise_Community_Enrichment_Results.csv", row.names = FALSE, sep="\t")


