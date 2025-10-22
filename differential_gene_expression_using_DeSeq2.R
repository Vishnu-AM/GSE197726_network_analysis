library(DESeq2)
library(ggplot2)
library(tidyverse)
library(airway)
library(biomaRt)
library(tibble)
library(ggrepel)
library(EnhancedVolcano)

# Load counts
counts <- read.csv("/home/shriyansh/Documents/rna_seq_analysis/GSE197726/gse197726_gene_symbol_counts.csv", 
                   header = TRUE, row.names = 1)


data <- read.delim("/home/shriyansh/Documents/IGF2BP2_PROJECT/clean_files_for_paper/geo_datasets/GSE197726/raw_counts_from_salmon/salmon_gse197726_raw_counts.tsv", header = TRUE)

print(colnames(data))


# Filter for WT and KO samples
counts_wt_ko <- counts[, 1:36]
counts_wt_ko <- round(counts_wt_ko)


# Extract sample names
sample_names <- colnames(counts_wt_ko)

# Assign condition
condition <- ifelse(grepl("^WT_", sample_names), "WT", "KO")
condition <- factor(condition, levels = c("WT", "KO"))

# Extract ZT timepoints
timepoint <- sub(".*ZT\\.([0-9]+).*", "ZT_\\1", sample_names)
timepoint <- factor(timepoint, levels = unique(timepoint))

# Metadata
coldata <- data.frame(
  row.names = sample_names,
  condition = condition,
  timepoint = timepoint
)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_wt_ko, colData = coldata, design = ~ condition)

# Filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


# Run DESeq2
dds$condition <- relevel(dds$condition, ref = "WT")
dds <- DESeq(dds)
res <- results(dds)
summary(res)

# Create output folder if it doesn't exist
output_dir_deseq <- "gse197726_salmon_deseq2_result"
if (!dir.exists(output_dir_deseq)) {
  dir.create(output_dir_deseq)
}

# Save DESeq2 results as CSV
output_file <- file.path(output_dir_deseq, "gse197726_salmon_deseq2_results.csv")
write.csv(as.data.frame(res), file = output_file)

# Convert results to dataframe
res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)
rownames(res_df) <- toupper(rownames(res_df))
res_df <- rownames_to_column(res_df, var = "gene")
head(res_df)

all_clock_related_genes <- c(
  "NR1D1", "NR1D2","IGF2BP2"
)

intersect(all_clock_related_genes, res_df$gene)
plot <- EnhancedVolcano(res_df,
                        x = "log2FoldChange",
                        y = "padj",
                        lab = res_df$gene,
                        pCutoff = 0.01,
                        FCcutoff = 1,
                        selectLab = all_clock_related_genes,
                        drawConnectors = TRUE,
                        boxedLabels = TRUE,
                        max.overlaps = Inf,
                        title = "Wildtype vs Knockout",
                        subtitle = NULL,
                        xlim = c(-5, 5)
)
plot + 
  theme(plot.title = element_text(hjust = 0.5))

clock_controlled_genes <- c("MYDO1"
)

intersect(all_clock_related_genes, res_df$gene)

plot <- EnhancedVolcano(res_df,
                        x = "log2FoldChange",
                        y = "padj",
                        lab = res_df$gene,
                        pCutoff = 0.01,
                        FCcutoff = 1,
                        selectLab = clock_controlled_genes,
                        drawConnectors = TRUE,
                        boxedLabels = TRUE,
                        max.overlaps = Inf,
                        title = "Wildtype vs Knockout",
                        subtitle = NULL,
                        xlim = c(-5, 5)
)
plot + 
  theme(plot.title = element_text(hjust = 0.5))

