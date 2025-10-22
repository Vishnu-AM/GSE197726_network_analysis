library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(dplyr)
library(stringr)
library(biomaRt)
library(patchwork)
library(gplots)
library(dendextend)
allowWGCNAThreads()          # allow multi-threading (optional)

# Set a common theme with a specific font family
theme_set(
  theme_minimal(base_size = 14, base_family = "Times New Roman")
)

#set the working directory
setwd("/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper")

#save the workspace image
save.image("/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/R_studio_codes/proper_codes/wgcna_salmon_counts_v2.RData")

#load the workspace image
load("/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/R_studio_codes/proper_codes/wgcna_salmon_counts_v2.RData")


# 1. Fetch Data ------------------------------------------------

data <- read.delim("/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/geo_datasets/GSE197726/raw_counts_from_salmon/salmon_gse197726_raw_counts.tsv", header = TRUE)

print(colnames(data))




#get metadata
geo_id <- "GSE197726"
gse_metadata <- getGEO (geo_id, GSEMatrix = TRUE)
phenoData<- pData (phenoData(gse_metadata[[1]]))
head (phenoData)
phenoData<- phenoData[c(1,2,11,12)]

# change the inconsistencies in PhenoData title column
phenoData$title <- gsub("_", "-", phenoData$title)
phenoData$title <- gsub(" ", "-", phenoData$title)


# Replace column names of data (excluding first column)
colnames(data)[-1] <- phenoData$title
colnames(data)

# prepare data
data[1:10,1:10]
print(colnames(data))

data <- data %>%
  pivot_longer(cols = -"ensembl_id", names_to = "samples", values_to = "counts") %>%
  inner_join(phenoData, by = c("samples" = "title")) %>%
  dplyr::select(1,2,3) %>%
  pivot_wider(names_from = samples, values_from = counts, names_sort = FALSE) %>%
  column_to_rownames(var = "ensembl_id")

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]


### NOTE: If there are batch effects observed, correct for them before moving ahead


# exclude outlier samples (NO SAMPLES WERE EXCLUDED)
#samples.to.be.excluded <- c()
#data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]


# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

# exclude outlier samples
colData <- phenoData


# fixing column names in colData
colnames(colData) <- c("title", "geo_accession", "genotype", "timepoint")

#cleaning up colData
colData <- colData %>%
  mutate(
    genotype = gsub("genotype: ", "", genotype),
    timepoint = gsub("time point: ", "", timepoint)
  )


#change the rownames
rownames(colData) <- colData$title

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data))
all(rownames(colData) == colnames(data))

#round the values of the data variable as deseq expects integers
countData<-round(data)

# create dds
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ 1) # not specifying model


## remove all genes with counts < 15 in more than 75% of samples (54*0.75=40.5 )
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 41,]
nrow(dds75) # 13918 genes


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

###### check for outliers
###### detect outlier samples - hierarchical clustering - method 1

outpath_hclust<- "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/outputs/01_htree"

###### Perform hierarchical clustering
htree <- hclust(dist(norm.counts), method = "average")

dend <- as.dendrogram(htree)

labels <- labels(dend)

# Extract group
group <- gsub("-.*", "", labels)

# Define colors (as you want)
group_colors <- c("WT" = "orange", "KO" = "darkblue", "mRE" = "brown")

# Map colors to each label
labels_colors(dend) <- group_colors[group]

# Plot as multiple file formats

# --- PNG ---
png(file.path(outpath_hclust, "htree_plot.png"), width = 6, height = 5, units = "in", res = 300, bg = "transparent")
par(mar = c(10,4,4,2), cex = 0.7)
plot(hang.dendrogram(dend), main = "Hierarchical Clustering", ylab = "Height")
legend("topright", legend = names(group_colors), fill = group_colors, border = "black", cex = 0.8)
par(cex = 1)
dev.off()

# --- PDF ---
pdf(file.path(outpath_hclust, "htree_plot.pdf"), width = 6, height = 5)
par(mar = c(10,4,4,2), cex = 0.7)
plot(hang.dendrogram(dend), main = "Hierarchical Clustering", ylab = "Height")
legend("topright", legend = names(group_colors), fill = group_colors, border = "black", cex = 0.8)
par(cex = 1)
dev.off()

# --- SVG ---
svg(file.path(outpath_hclust, "htree_plot.svg"), width = 6, height = 5)
par(mar = c(10,4,4,2), cex = 0.7)
plot(hang.dendrogram(dend), main = "Hierarchical Clustering", ylab = "Height")
legend("topright", legend = names(group_colors), fill = group_colors, border = "black", cex = 0.8)
par(cex = 1)
dev.off()

# --- TIFF ---
tiff(file.path(outpath_hclust, "htree_plot.tiff"), width = 6, height = 5, units = "in", res = 300)
par(mar = c(10,4,4,2), cex = 0.7)
plot(hang.dendrogram(dend), main = "Hierarchical Clustering", ylab = "Height")
legend("topright", legend = names(group_colors), fill = group_colors, border = "black", cex = 0.8)
par(cex = 1)
dev.off()


###### detect outlier samples - Principal Component Analysis (PCA) - method 2
output_dir_pca <- "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/outputs/01.1_pca"

pca <- prcomp(norm.counts)
pca.dat <- as.data.frame(pca$x)
pca.dat$samples <- rownames(pca.dat)


# Replace this with your real metadata
metadata_df <- data.frame(samples = rownames(norm.counts), genotype = c("WT", "KO", "mRE"))

# Merge PCA data with metadata
pca.dat <- pca.dat %>%
  left_join(metadata_df, by = "samples")

# Define publication-friendly colors
genotype_colors <- c("WT" = "orange", "KO" = "darkblue", "mRE" = "brown")


# Create the PCA plot and assign to a variable
pca_plot <- ggplot(pca.dat, aes(x = PC1, y = PC2, color = genotype)) +
  geom_point(size = 3) +
  geom_text(aes(label = samples), size = 4, vjust = -1, show.legend = FALSE) +
  scale_color_manual(values = genotype_colors) +
  labs(
    x = paste0("PC1: ", pca.var.percent[1], "%"),
    y = paste0("PC2: ", pca.var.percent[2], "%"),
    color = "Genotype"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )

# Save in multiple formats
ggsave(file.path(output_dir_pca, "pca_plot.png"), plot = pca_plot, width = 10, height = 10, dpi = 300, bg = "transparent")
ggsave(file.path(output_dir_pca, "pca_plot.pdf"), plot = pca_plot, width = 10, height = 10, bg = "transparent")
ggsave(file.path(output_dir_pca, "pca_plot.svg"), plot = pca_plot, width = 10, height = 10, bg = "transparent")
ggsave(file.path(output_dir_pca, "pca_plot.tiff"), plot = pca_plot, width = 10, height = 10, dpi = 300, bg = "transparent")


# --- Scree Plot Data ---
cumulative_var <- round(cumsum(pca.var.percent), 2)
scree_data <- data.frame(
  PC = factor(paste0("PC", 1:length(pca.var.percent)), levels = paste0("PC", 1:length(pca.var.percent))),
  Variance = pca.var.percent,
  Cumulative = cumulative_var
)

# Scree bar plot
scree_bar_plot <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  labs(x = NULL, y = "Variance Explained (%)") +  # remove x-axis title
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 12),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

# Cumulative variance line plot
scree_line_plot <- ggplot(scree_data, aes(x = PC, y = Cumulative, group = 1)) +
  geom_line(color = "darkblue", size = 1) +
  geom_point(color = "darkblue", size = 2) +
  labs(x = NULL, y = "Cumulative Variance (%)") +  # keep only bottom plot x-axis
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.title = element_text(size = 12),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

# Stack vertically
stacked_plot <- scree_bar_plot / scree_line_plot
print(stacked_plot)

# --- Save stacked plot in multiple formats ---
ggsave(file.path(output_dir_pca, "scree_stacked.png"), plot = stacked_plot, width = 10, height = 10, dpi = 300, bg = "transparent")
ggsave(file.path(output_dir_pca, "scree_stacked.pdf"), plot = stacked_plot, width = 10, height = 10, bg = "transparent")
ggsave(file.path(output_dir_pca, "scree_stacked.svg"), plot = stacked_plot, width = 10, height = 10, bg = "transparent")
ggsave(file.path(output_dir_pca, "scree_stacked.tiff"), plot = stacked_plot, width = 10, height = 10, dpi = 300, bg = "transparent")



###### Plot the VST normalized counts

output_dir_vst_exp<- "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/outputs/02_vst_expression_plot"
# Your norm.counts: rows = samples, columns = genes
samples <- rownames(norm.counts)  # 54 samples

# Define genotype for each sample
genotype <- rep(c("WT", "KO", "mRE"), each = 18)  # 18 samples per genotype

# Create metadata dataframe
metadata_df <- data.frame(
  Sample = samples,
  Genotype = genotype
)

# Convert norm.counts to long format for ggplot
expr_normalized_df <- norm.counts %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(
    -Sample, 
    names_to = "Gene", 
    values_to = "Expression"
  ) %>%
  left_join(metadata_df, by = "Sample") %>%  # Add genotype info
  mutate(Sample = factor(Sample, levels = samples))  # Keep original order

# Violin plot with legend
norm_exp_violin_plot <- ggplot(expr_normalized_df, aes(x = Sample, y = Expression, fill = Genotype)) +
  geom_violin() +
  geom_point(alpha = 0.5, shape = 22, size = 2, color = "black") +
  scale_fill_manual(values = c("WT" = "orange", "KO" = "darkblue", "mRE" = "brown")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(
    title = "Normalized expression",
    x = "Samples",
    y = "Expression"
  )

# Create output directory if it doesn't exist
dir.create(output_dir_vst_exp, showWarnings = FALSE, recursive = TRUE)

# Save in PNG
ggsave(
  filename = file.path(output_dir_vst_exp, "vst_expression_violin.png"),
  plot = norm_exp_violin_plot,
  width = 12, height = 10, dpi = 300
)

# Save in PDF
ggsave(
  filename = file.path(output_dir_vst_exp, "vst_expression_violin.pdf"),
  plot = norm_exp_violin_plot,
  width = 12, height = 10
)

# Save in SVG
ggsave(
  filename = file.path(output_dir_vst_exp, "vst_expression_violin.svg"),
  plot = norm_exp_violin_plot,
  width = 12, height = 10
)

# Save in TIFF
ggsave(
  filename = file.path(output_dir_vst_exp, "vst_expression_violin.tiff"),
  plot = norm_exp_violin_plot,
  width = 12, height = 10, dpi = 300
)

# 4. Network Construction by WGCNA ---------------------------------------------------

# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# Assuming sft.data and log_mean_k are already defined
log_mean_k <- log10(sft.data$mean.k.)

# Plot 1: Scale-free topology fit
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.035) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()+
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

print(a1)

# Plot 2: Mean connectivity
a2 <- ggplot(sft.data, aes(Power, log_mean_k, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.11) +
  labs(x = 'Power', y = 'Mean Connectivity (log10)') +
  theme_classic()+
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

print(a2)

#plot the sft
combined_sft_plot<- a1/a2

output_sft_plot<- "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/outputs/03_sft_plot"

# Create output directory if it doesn't exist
dir.create(output_sft_plot, showWarnings = FALSE, recursive = TRUE)

# PNG
ggsave(
  filename = file.path(output_sft_plot, "sft_plot.png"),
  plot = combined_sft_plot,
  width = 10, height = 12, dpi = 300
)

# PDF
ggsave(
  filename = file.path(output_sft_plot, "sft_plot.pdf"),
  plot = combined_sft_plot,
  width = 10, height = 12
)

# SVG
ggsave(
  filename = file.path(output_sft_plot, "sft_plot.svg"),
  plot = combined_sft_plot,
  width = 10, height = 12
)

# TIFF
ggsave(
  filename = file.path(output_sft_plot, "sft_plot.tiff"),
  plot = combined_sft_plot,
  width = 10, height = 12, dpi = 300
)


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          deepSplit = 2,
                          minModuleSize = 50,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor

## Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# PLOT THE NUMBER OF GENES PER MODULE
# Convert table to data frame and sort
genes_per_module <- as.data.frame(table(bwnet$colors)) %>%
  rename(Module = Var1, Count = Freq) %>%
  arrange(desc(Count))  
print(genes_per_module)

# Define output folder
out_folder_module_plot <- "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/outputs/04_module_size"  # Replace with your actual path

# Filter out grey module
genes_per_module_filtered <- genes_per_module %>%
  filter(Module != "grey")

# Update color mapping (remove grey)
module_colors <- setNames(unique(bwnet$colors[bwnet$colors != "grey"]), 
                          unique(bwnet$colors[bwnet$colors != "grey"]))

# Plot
module_plot <- ggplot(genes_per_module_filtered, aes(x = reorder(Module, -Count), y = Count, fill = Module)) +
  geom_bar(stat = "identity", color = "black") +  # black outline for bars
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)  # adds box
  ) +
  scale_fill_manual(values = module_colors) +
  scale_y_continuous(breaks = seq(0, max(genes_per_module_filtered$Count), by = 1000)) +
  labs(
    x = "Modules",
    y = "No. of Genes"
  ) +
  guides(fill = "none") +
  geom_text(aes(label = Count), vjust = -0.5, size = 3)+theme(
    axis.title.x = element_text(face = "bold", size = 12),  # bold x label
    axis.title.y = element_text(face = "bold", size = 12)   # bold y label
  )

print(module_plot)



# Save the plot in multiple formats
ggsave(file.path(out_folder_module_plot, "module_sizes_plot.png"), plot = module_plot, bg = "transparent", width = 12, height = 10, dpi = 300)
ggsave(file.path(out_folder_module_plot, "module_sizes_plot.svg"), plot = module_plot, bg = "transparent", width = 12, height = 10)
ggsave(file.path(out_folder_module_plot, "module_sizes_plot.pdf"), plot = module_plot, bg = "transparent", width = 12, height = 10)
ggsave(file.path(out_folder_module_plot, "module_sizes_plot.tiff"), plot = module_plot, bg = "transparent", width = 12, height = 10, dpi = 300, device = "tiff")





# Define output folder
out_folder_dendro_plot <- "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/outputs/05_cluster_dendro"

# Define the formats and corresponding functions
formats <- list(
  png = function(file) png(file, width = 1800, height = 1500, res = 300, bg = "transparent"),
  tiff = function(file) tiff(file, width = 1800, height = 1500, res = 300, bg = "transparent"),
  svg = function(file) svg(file, width = 6, height = 5, bg = "transparent"),
  pdf = function(file) pdf(file, width = 6, height = 5, bg = "transparent")
)
# Loop over different file formats and save the dendrogram plot
for (ext in names(formats)) {
  
  # Create the file path for the output
  file_path <- file.path(out_folder_dendro_plot, paste0("module_dendrogram_plot.", ext))
  
  # Open the correct graphics device for the format
  formats[[ext]](file_path)
  
  # Plot the dendrogram with module colors
  plotDendroAndColors(
    bwnet$dendrograms[[1]],
    bwnet$colors,
    "Module colors",
    dendroLabels = FALSE,
    addGuide = TRUE,
    hang = 0.03,
    guideHang = 0.05,
    main = ""
  )
  
  # Close the device to save the file
  dev.off()
}

# Also print the dendrogram directly in RStudio
plotDendroAndColors(
  bwnet$dendrograms[[1]],
  bwnet$colors,
  "Module colors",
  dendroLabels = FALSE,
  addGuide = TRUE,
  hang = 0.03,
  guideHang = 0.05,
  main = ""
)

print(bwnet$dendrograms)
# grey module = all genes that doesn't fall into other modules were assigned to the grey module


# Create a single data frame
module_df <- data.frame(
  ensembl_id = names(bwnet$colors),
  Module = bwnet$colors
)

# Define output path
output_file <- "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/output_text_files/03_module_info/all_modules.tsv"

# Save as TSV
write_tsv(module_df, file = output_file)


###### Visualize the adjacency matrix

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)
moduleColors = bwnet$colors
geneTree = bwnet$dendrograms[[1]]
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(norm.counts, power = 18);

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^18;

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;

#set color
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')



tomplot<- TOMplot(plotTOM, geneTree, moduleColors, main = "", col = myheatcol)



###### Plot the module exoression profile
# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(norm.counts, bwnet$colors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$genotype = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-genotype) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )


# Convert genotype into a factor to preserve its order
mME$genotype <- factor(mME$genotype, levels = unique(mME$genotype))

# Define output folder
out_mME <- "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/outputs/06_module_trait_association"  # Replace with your actual path

# Create the plot
mME_plot <- mME %>%
  ggplot(aes(x = genotype, y = name, fill = value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1)
  ) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(
    y = "Modules",
    fill = NULL
  )

print(mME_plot)
# Save in four formats with transparent background
ggsave(file.path(out_mME, "module_trait_heatmap.png"), plot = mME_plot, bg = "transparent", width = 9, height = 5, dpi = 300)
ggsave(file.path(out_mME, "module_trait_heatmap.tiff"), plot = mME_plot, bg = "transparent", width = 9, height = 5, dpi = 300)
ggsave(file.path(out_mME, "module_trait_heatmap.svg"), plot = mME_plot, bg = "transparent", width = 9, height = 5)
ggsave(file.path(out_mME, "module_trait_heatmap.pdf"), plot = mME_plot, bg = "transparent", width =9, height = 5)


#PLOT THE EXPRESSION PROFILE OF THE SELECTED MODULES##############################################33

# Define your output base directory
output_dir <- "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/outputs/07_exp_of_select_modules"
dir.create(output_dir, showWarnings = FALSE)

# Get all unique modules
modules <- unique(module_df$colors)

# Set row names and transpose expression matrix
row.names(module_df) <- module_df$ensembl_id
norm.counts_r_t <- t(norm.counts)

# Loop through each module and create plots
for (mod in modules) {
  
  # Filter genes of the current module
  submod <- module_df %>% filter(colors == mod)
  subexpr <- norm.counts_r_t[submod$ensembl_id, ]
  
  # Tidy data
  submod_df <- data.frame(subexpr) %>%
    mutate(ensembl_id = row.names(.)) %>%
    pivot_longer(-ensembl_id) %>%
    mutate(module = module_df[ensembl_id, ]$colors)
  
  submod_df$module <- factor(submod_df$module, levels = mod)
  submod_df$name <- factor(submod_df$name, levels = unique(submod_df$name))
  
  # Color mapping for this module only
  color_mapping <- setNames(mod, mod)
  
  # Create plot
  p <- ggplot(submod_df, aes(x = name, y = value, group = ensembl_id, color = module)) +
    geom_line(alpha = 0.2) +
    theme_minimal(base_family = "sans") +
    theme(
      axis.text.x = element_text(angle = 90),
      legend.position = "none",
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA)
    ) +
    scale_color_manual(values = color_mapping) +
    labs(
      x = "Genotype",
      y = "Normalized Expression"
    )
  
  # Make output subdir for this module
  module_path <- file.path(output_dir, mod)
  dir.create(module_path, showWarnings = FALSE)
  
  # Define file paths
  file_base <- file.path(module_path, paste0("module_", mod))
  ggsave(paste0(file_base, ".png"), plot = p, width = 8, height = 6, dpi = 300, bg = "transparent")
  ggsave(paste0(file_base, ".tiff"), plot = p, width = 8, height = 6, dpi = 300, bg = "transparent", device = "tiff")
  ggsave(paste0(file_base, ".pdf"), plot = p, width = 8, height = 6, bg = "transparent")
}
# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

head(colData)
# binarize categorical variables

# Convert genotype to numeric indicators
colData$genotype_WT <- ifelse(colData$genotype == "WT", 1, 0)
colData$genotype_KO <- ifelse(colData$genotype == "KO", 1, 0)
colData$genotype_mRE <- ifelse(colData$genotype == "mRE", 1, 0)

colnames(colData) <- gsub("genotype_", "Genotype-", colnames(colData))


traits <- colData[, c("Genotype-WT", "Genotype-KO", "Genotype-mRE")]

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)



# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
print (heatmap.data)
head(heatmap.data)
colnames(heatmap.data)
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

# Remove "ME" and "Genotype-" prefixes from column names
colnames(heatmap.data) <- gsub("^ME", "", colnames(heatmap.data))
colnames(heatmap.data) <- gsub("^Genotype-", "", colnames(heatmap.data))

# Check the result
print(colnames(heatmap.data))
print(MEs0$MEblack)


# Define output directory
out_dir_mod_trait_bi <- "/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/outputs/08_module_trait_relationship_binarized"
dir.create(out_dir_mod_trait_bi, showWarnings = FALSE)

# Define file formats
formats <- c("svg", "png", "tiff", "pdf")

# Create the plot once
p <- CorLevelPlot(
  heatmap.data,
  x = names(heatmap.data)[9:11],
  y = names(heatmap.data)[1:8],
  col = c("#2C7BB6","#ABD9E9","#FFFFBF","#FDAE61","#D7191C")
)

print(p)
# Loop through formats and save
for (fmt in formats) {
  file_path <- file.path(out_dir_mod_trait_bi, paste0("correlation_plot.", fmt))
  
  if (fmt == "svg") svg(file_path, width = 8, height = 6, bg = "transparent")
  if (fmt == "png") png(file_path, width = 800, height = 600, res = 150, bg = "transparent")
  if (fmt == "tiff") tiff(file_path, width = 1000, height = 800, res = 150, bg = "transparent")
  if (fmt == "pdf") pdf(file_path, width = 8, height = 6, bg = "transparent")
  
  print(p)  
  
  dev.off()
}
######################TESTING THE BINARIZATION#################################
MEturquoise_test <- c(
  -0.06352805, -0.10264834, -0.06417784, -0.08946219, -0.08775759, -0.08419830,
  -0.10248717, -0.10754553, -0.11406159, -0.13785260, -0.12766519, -0.13276890,
  -0.10798781, -0.12005800, -0.12101435, -0.07962311, -0.09353559, -0.07927985,
  0.15246734,  0.19248927,  0.20918906,  0.16245766,  0.19360947,  0.18729832,
  0.16931182,  0.19595992,  0.20439261,  0.17874269,  0.17562046,  0.21102987,
  0.19411664,  0.20039505,  0.18896512,  0.21105819,  0.18198761,  0.21562181,
  -0.05363609, -0.05655736, -0.10023877, -0.08664692, -0.07587580, -0.07880171,
  -0.13755735, -0.08505800, -0.11902263, -0.08342245, -0.10083649, -0.08063795,
  -0.12503861, -0.08642848, -0.09109815, -0.08373223, -0.08618085, -0.07829106
)

genotype_test <- c(rep(1, 18), rep(0, 18), rep(0, 18))

correlation_test <- cor(MEturquoise_test, genotype_test)
print(correlation_test)
###############################################################################

#extract the genes in selected modules
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()



# 6B. Intramodular analysis: Identifying driver genes ---------------

traits <- colData[, c("Genotype-WT", "Genotype-KO", "Genotype-mRE")]

# Get module names
modNames <- substring(names(module_eigengenes), 3)

for (trait in names(traits)) {
  
  # Extract the trait data
  trait_data <- as.data.frame(colData[[trait]])
  names(trait_data) <- trait
  
  # Compute module membership and significance
  geneModuleMembership <- as.data.frame(cor(norm.counts, module_eigengenes, use = 'p'))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste("MM", modNames, sep = "")
  names(MMPvalue) <- paste("p.MM", modNames, sep = "")
  
  geneTraitSignificance <- as.data.frame(cor(norm.counts, trait_data, use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) <- paste("GS.", trait, sep = "")
  names(GSPvalue) <- paste("p.GS.", trait, sep = "")
  
  # Combine all p-values to find the smallest non-zero value
  all_pvals <- c(as.matrix(GSPvalue), as.matrix(MMPvalue))
  min_nonzero <- min(all_pvals[all_pvals > 0])
  
  # Replace 0s with the smallest non-zero value
  GSPvalue[GSPvalue == 0] <- min_nonzero
  MMPvalue[MMPvalue == 0] <- min_nonzero
  

  
  # Create a dataframe with annotation and modules
  probes <- colnames(norm.counts)
  #probes2annot <- match(probes, annot$ILMN_ID)
  
  geneInfo0 <- data.frame( probes,
                          moduleColor = moduleColors,
                          geneTraitSignificance,
                          GSPvalue)
  
  # Order modules by significance for the trait
  modOrder <- order(-abs(cor(module_eigengenes, trait_data, use = 'p')))
  
  # Add module membership information
  for (mod in 1:ncol(geneModuleMembership)) {
    oldnames <- names(geneInfo0)
    geneInfo0 <- data.frame(geneInfo0,
                            geneModuleMembership[, modOrder[mod]],
                            MMPvalue[, modOrder[mod]])
    names(geneInfo0) <- c(oldnames,
                          paste("MM.", modNames[modOrder[mod]], sep = ""),
                          paste("p.MM.", modNames[modOrder[mod]], sep = ""))
  }
  
  # Rearrange data and save to CSV
  geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0[[paste0("GS.", gsub(" ", ".", trait))]]))
  geneInfo <- geneInfo0[geneOrder,]
  
  write.csv(geneInfo, file = paste0("/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/output_text_files/04_GS_vs_MM_gene_info", gsub(" ", "_", tolower(trait)), ".csv"), row.names = FALSE)
  write.table(geneInfo, 
              file = paste0("/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/output_text_files/04_GS_vs_MM_gene_info", 
                            gsub(" ", "_", tolower(trait)), ".tsv"), 
              row.names = FALSE, 
              sep = "\t",   # Use tab as separator
              quote = FALSE)  # Remove quotes from strings
}



#################################################################################

selectedModules <- c("turquoise", "green", "blue", "yellow")
selectedTraits <- c("Genotype-WT", "Genotype-KO", "Genotype-mRE")

xValsList <- list()
yValsList <- list()
geneTraitSignificanceList <- list()
GSPvalueList <- list()

for (trait in selectedTraits) {
  if (!trait %in% names(colData)) {
    stop(paste("Trait", trait, "not found in datTraits"))
  }
  
  trait_data <- as.data.frame(colData[[trait]])
  names(trait_data) <- trait
  
  geneTraitSignificance <- as.data.frame(cor(norm.counts, trait_data, use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  
  names(geneTraitSignificance) <- paste("GS.", trait, sep = "")
  names(GSPvalue) <- paste("p.GS.", trait, sep = "")
  
  geneTraitSignificanceList[[trait]] <- geneTraitSignificance
  GSPvalueList[[trait]] <- GSPvalue
  
  xValsList[[trait]] <- c()
  yValsList[[trait]] <- c()
  
  for (module in selectedModules) {
    column <- match(module, modNames)
    if (is.na(column)) next
    
    moduleGenes <- moduleColors == module
    xValsList[[trait]] <- c(xValsList[[trait]], abs(geneModuleMembership[moduleGenes, column]))
    yValsList[[trait]] <- c(yValsList[[trait]], abs(geneTraitSignificance[moduleGenes, 1]))
  }
}

# Set plot limits
xMin <- min(unlist(xValsList), na.rm = TRUE)
xMax <- max(unlist(xValsList), na.rm = TRUE)
yMin <- min(unlist(yValsList), na.rm = TRUE)
yMax <- max(unlist(yValsList), na.rm = TRUE)

xTicks <- pretty(c(xMin, xMax))
yTicks <- pretty(c(yMin, yMax))

# Plotting for selected traits and modules
sizeGrWindow(14, 3 * length(selectedTraits))
par(mfrow = c(length(selectedTraits), length(selectedModules)), mar = c(4, 4, 3, 1), oma = c(1, 1, 1, 1))

for (trait in selectedTraits) {
  for (module in selectedModules) {
    column <- match(module, modNames)
    if (is.na(column)) next
    
    moduleGenes <- moduleColors == module
    
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificanceList[[trait]][moduleGenes, 1]),
                       xlab = "Module membership",
                       ylab = "Gene significance",
                       main = paste("MM vs. GS\n(", trait,")"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,
                       xlim = c(0, 1), ylim = c(0, 1))
  }
}



########################################################################33333
#identify the genes that significantly correlating with WT
gene.signf.corr.mRE <- cor(norm.counts, traits$genotype_mRE, use = 'p')
gene.signf.corr.pvals.mRE <- corPvalueStudent(gene.signf.corr.mRE, nSamples)


top_25_genes_associated_with_mRE<-gene.signf.corr.pvals.mRE %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)

#save the file
write.table(top_25_genes_associated_with_mRE, "/home/vishnu/transcriptomics_for_thesis_with_shriyansh/wgcna_output_salmon/wgcna_binarized_implementation_v2/top_25_genes_v2/mRE_v2/top_25_genes_mRE_v2.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

#EXTRACT THE EDGELIST FOR SELECTED MODULES-----------------------------------------------------------------------------------------------------------------------------------
# pick out a few modules of interest here
modules_of_interest = c("turquoise")

#extracting networks
genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = t(norm.counts)[genes_of_interest$ensembl_id,]
expr_of_interest[1:5,1:5]

# Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = 18)

# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

#get the edgelist
edge_list_turqoise = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

edge_list_unique_turqoise <- edge_list_turqoise %>%
  mutate(
    gene1_fixed = pmin(gene1, gene2),
    gene2_fixed = pmax(gene1, gene2)
  ) %>%
  dplyr::select(gene1_fixed, gene2_fixed, correlation, module1, module2) %>%
  distinct()

# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list_turqoise,
            file = "/home/vishnu/transcriptomics_for_thesis_with_shriyansh/wgcna_output_salmon/wgcna_binarized_implementation_v2/edge_list_v2/turqoise/turqoise_edgelist_v2.tsv",
            delim = "\t")

write_delim(edge_list_unique_turqoise,
            file = "/home/vishnu/transcriptomics_for_thesis_with_shriyansh/wgcna_output_salmon/wgcna_binarized_implementation_v2/edge_list_v2/turqoise/unique_turqoise_edgelist_v2.tsv",
            delim = "\t")

