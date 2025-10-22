library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(dplyr)
library(stringr)
library(biomaRt)
library(dryR)
library(ggplot2)  # Ensure ggplot2 is loaded for ggsave
library(biomaRt)
library(reshape2)
library(ComplexHeatmap)
library(circlize)


# Load expression data
countData <- read.delim("/home/vishnu/transcriptomics_for_thesis_with_shriyansh/geo_datasets/GSE197726/subsetted_expression_module_wise/turquoise_exp/turquoise_expression.tsv")

# Set 'ensembl_id' as row names and remove the column
rownames(countData) <- countData$ensembl_id
countData$ensembl_id <- NULL

# Format the column names
colnames(countData) <- gsub("\\.", "_", colnames(countData))

# Define time points (0, 4, 8, 12, 16, 20 repeated for 3 replicates across 3 conditions)
time_points <- c(0, 4, 8, 12, 16, 20)
time_vector <- rep(time_points, each = 3)
time <- rep(time_vector, times = 3)
print(time)

# Define group factor from column names
columns <- colnames(countData)
group <- sub("_ZT.*", "", columns)
group <- factor(group)
print(group)

# Connect to Ensembl for mouse data
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Map Ensembl IDs to gene symbols and Entrez IDs
ensembl_ids <- rownames(countData)
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# Add ensembl_gene_id back to countData and merge with gene mapping
countData$ensembl_gene_id <- rownames(countData)
countData_merged <- merge(gene_map, countData, by = "ensembl_gene_id")

# Replace empty gene symbols with "unknown"
countData_merged$mgi_symbol[countData_merged$mgi_symbol == ""] <- "unknown"

# Make gene symbols unique (handle duplicates)
countData_merged$mgi_symbol <- make.unique(countData_merged$mgi_symbol, sep = "_")

# Set gene symbol as row names
rownames(countData_merged) <- countData_merged$mgi_symbol

# Round expression counts (excluding mapping columns)
expr_cols <- setdiff(names(countData_merged), c("ensembl_gene_id", "entrezgene_id", "mgi_symbol"))
countData_merged[, expr_cols] <- round(countData_merged[, expr_cols])

# Prefilter: Keep rows where row sum (only expression) > 0
keep <- rowSums(countData_merged[, expr_cols]) > 0
countData_final <- countData_merged[keep, ]

# Done!
head(countData_final)

# Remove annotation columns and assign to countData_dryr
countData_dryr <- countData_final[, !(colnames(countData_final) %in% c("ensembl_gene_id", "entrezgene_id", "mgi_symbol"))]

# run the analysis for count data (e.g. RNA-Seq data)
dryList   = dryseq(countData_dryr,group,time)

# explore the results
dryList[["results"]]     # data frame summarizing results
dryList[["parameters"]]  # coefficients: phase, amplitude and mean for each group
dryList[["ncounts"]]     # normalized counts
dryList[["counts"]]      # raw counts
dryList[["cook"]]        # cook's distance for outlier detection
dryList[["BICW_rhythm"]] # BICW for each rhythmic model
dryList[["BICW_mean"]]   # BICW for each mean model

# generate a pdf with a global summary of all models
plot_models_rhythm(dryList, "/home/vishnu/transcriptomics_for_thesis_with_shriyansh/dryR_output/gse_197726_turquoise_dry_r")

# plot a feature of interest
dry_plot(dryList, "ENSMUSG00000002489")
# Extract the results dataframe
results_df <- dryList$results

models<- dryList$results
# Convert row names to an explicit column
results_df <- data.frame(ensembl_id = rownames(results_df), results_df, check.names = FALSE)

# Save as a TSV file
write.table(results_df, "/home/vishnu/transcriptomics_for_thesis_with_shriyansh/dryR_output/turquoise/results_file/dryR_turquoise_module_results.tsv", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#plot the expression for top 25 genes by degree

# Define the input CSV file
csv_file <- "/home/vishnu/transcriptomics_for_thesis_with_shriyansh/network_analysis_for_v2_files/turquoise/degree/top_25_nodes_by_degree.csv"  # Change this to your actual file path

# Read the CSV file (assuming it has a 'Node' column with ENSEMBL IDs)
gene_data <- read.delim(csv_file, stringsAsFactors = FALSE)

# Extract the list of genes from the 'Node' column
genes <- gene_data$gene_id

gene_symbol<- gene_data$gene_id

# Define the base output directory
output_dir <- "/home/vishnu/transcriptomics_for_thesis_with_shriyansh/dryR_output/turquoise/expression_plot/top_25_by_degree/"

# Loop through each gene
for (gene in genes) {
  
  # Create a subfolder for the gene
  gene_dir <- file.path(output_dir, gene)
  if (!dir.exists(gene_dir)) {
    dir.create(gene_dir, recursive = TRUE)
  }
  
  # Generate the plot
  plot_obj <- dry_plot(dryList, gene)  # Assumes dry_plot returns a ggplot object
  
  # Save the plot in different formats
  ggsave(filename = file.path(gene_dir, paste0(gene, ".svg")), plot = plot_obj, device = "svg")
  ggsave(filename = file.path(gene_dir, paste0(gene, ".png")), plot = plot_obj, device = "png")
  ggsave(filename = file.path(gene_dir, paste0(gene, ".tiff")), plot = plot_obj, device = "tiff", dpi = 300)
  ggsave(filename = file.path(gene_dir, paste0(gene, ".pdf")), plot = plot_obj, device = "pdf")
  
  message(paste("Saved plots for", gene, "in", gene_dir))
}

# Filter dryR results to keep only the genes in the top 25 list
top_25_dryr_results <- results_df[results_df$ensembl_id %in% genes, ]

# Print a preview of the filtered data
head(top_25_dryr_results)

# Save the subset as a new TSV file
write.table(top_25_dryr_results, 
            "/home/vishnu/transcriptomics_for_thesis_with_shriyansh/dryR_output/turquoise/results_file/dryR_turquoise_top_25.tsv", 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Function to determine colors based on amplitude and phase
get_color <- function(ref_amp, ref_phase, test_amp, test_phase) {
  if (is.na(ref_amp) | is.na(ref_phase)) {
    return("white")  # WT is NA, everything should be white
  }
  if (is.na(test_amp) | is.na(test_phase)) {
    return("white")  # If test (KO/mRE) has NA, keep it white
  }
  if (test_amp == ref_amp & test_phase == ref_phase) {
    return("orange")  # Matches WT completely
  } else if (test_amp != ref_amp & test_phase == ref_phase) {
    return("darkgreen")  # Different amp, same phase
  } else {
    return("blue")  # Both amp and phase different
  }
}

# Create a heatmap matrix
heatmap_data <- data.frame(
  ensembl_id = top_25_dryr_results$ensembl_id,
  WT = "orange",
  KO = mapply(get_color, top_25_dryr_results$amp_WT, top_25_dryr_results$phase_WT, 
              top_25_dryr_results$amp_KO, top_25_dryr_results$phase_KO),
  mRE = mapply(get_color, top_25_dryr_results$amp_WT, top_25_dryr_results$phase_WT, 
               top_25_dryr_results$amp_mRE, top_25_dryr_results$phase_mRE)
)

# If WT has NA, set all three columns (WT, KO, mRE) to "white"
heatmap_data[is.na(top_25_dryr_results$amp_WT) | is.na(top_25_dryr_results$phase_WT), 2:4] <- "white"

# Melt the data for ggplot2
heatmap_melt <- melt(heatmap_data, id.vars = "ensembl_id")
# Correct the factor order for the 'value' column
heatmap_melt$value <- factor(heatmap_melt$value, levels = c("orange", "blue", "white"))

plot <- ggplot(heatmap_melt, aes(x = variable, y = ensembl_id, fill = value)) +
  geom_tile(color = "black", width = 0.95, height = 0.98, linewidth = 0.2) +  # Nearly full cell coverage
  scale_fill_manual(
    values = c(
      "orange" = "orange",        # Normal
      "blue" = "blue",            # Alt. Amp & Phase
      "white" = "white"           # Non-rhythmic
    ),
    name = "dryR",  # Legend title
    labels = c(
      "Normal", 
      "Alt. Amp & Phase", 
      "Non-rhythmic"
    )  # Legend labels in the desired order
  ) +
  theme_minimal() +
  labs(title = "Gene Comparison Across WT, KO, mRE", x = "", y = "Gene ID") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",  # Place the legend on the right
    legend.title = element_text(size = 12, color = "black"),  # Increase legend title size and bold
    legend.text = element_text(size = 12, color = "black"),   # Increase legend text size and bold
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_fixed(ratio = 1)  # ratio 1 gives square tiles
# Check if the plot object is valid by printing it
print(plot)

# Save as SVG with transparent background
ggsave("/home/vishnu/transcriptomics_for_thesis_with_shriyansh/dryR_output/turquoise/heatmap_plot.svg", plot = plot, device = "svg", bg = "transparent", width = 8, height = 6)

# Save as TIFF with transparent background and 300 dpi
ggsave("/home/vishnu/transcriptomics_for_thesis_with_shriyansh/dryR_output/turquoise/heatmap_plot.tiff", plot = plot, device = "tiff", bg = "transparent", dpi = 300, width = 8, height = 6)


#SAME PLOT WITHOUT LEGEND
plot_1 <- ggplot(heatmap_melt, aes(x = variable, y = ensembl_id, fill = value)) +
  geom_tile(color = "black", width = 0.95, height = 0.98, linewidth = 0.2) +  # Nearly full cell coverage
  scale_fill_manual(
    values = c(
      "orange" = "orange",        # Normal
      "blue" = "blue",            # Alt. Amp & Phase
      "white" = "white"           # Non-rhythmic
    )
  ) +
  theme_minimal() +
  labs(title = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none",  # <---- Remove legend
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_fixed(ratio = 1)

# Check the plot
print(plot_1)

# Save without legend
ggsave("/home/vishnu/transcriptomics_for_thesis_with_shriyansh/dryR_output/turquoise/heatmap_plot_no_legend.svg", plot = plot_1, device = "svg", bg = "transparent", width = 8, height = 6)
ggsave("/home/vishnu/transcriptomics_for_thesis_with_shriyansh/dryR_output/turquoise/heatmap_plot_no_legend.tiff", plot = plot_1, device = "tiff", bg = "transparent", dpi = 300, width = 8, height = 6)

