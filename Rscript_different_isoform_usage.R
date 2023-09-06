
# LRonly

setwd("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/diffIsoformUsage/data_LRonly/")

# Prepare data

diff_isos.cytosolVSnucleus.filtered <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/diffIsoformUsage/data_LRonly/known_diffIsoformUsage/diff_isos.cytosolVSnucleus.filtered.txt", header=FALSE)
colnames(diff_isos.cytosolVSnucleus.filtered) <- c("gene name", "isoform name", "pvalue", "sample1 isoform count", 
                    "sample2 isoform count", "sample1 alternative isoforms for gene count",
                    "sample2 alternative isoforms for gene count")
# Filter data for p-value < 0.05
filtered_data <- diff_isos.cytosolVSnucleus.filtered[diff_isos.cytosolVSnucleus.filtered$pvalue < 0.05, ]
# Filter out rows with NA values
filtered_data <- filtered_data[complete.cases(filtered_data), ]
colnames(filtered_data) <- c("gene name", "isoform name", "pvalue", "sample1 isoform count", 
                                                   "sample2 isoform count", "sample1 alternative isoforms for gene count",
                                                   "sample2 alternative isoforms for gene count")

# Compute adjusted p-values using the column 'pvalue'
filtered_data$adjusted_pvalue <- p.adjust(filtered_data$pvalue, method = "bonferroni")
filtered_data_adjustedPval <- filtered_data[filtered_data$adjusted_pvalue < 0.05, ]

# Specify the file path to save the TSV file in current directory
file_path <- "filtered_data_adjustedPval.tsv"

# Save the data with headers to a TSV file
write.table(filtered_data_adjustedPval, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# GO enrichment analysis:

# Load necessary packages
library(clusterProfiler)
library(org.Hs.eg.db)

# Create a data frame with gene symbols, p-values, and adjusted p-values

selected_columns <- c("gene name", "pvalue", "adjusted_pvalue")
gene_data <- filtered_data_adjustedPval[selected_columns]
colnames(gene_data)[1] <- "gene"

ensembl_genes <- gene_data$gene
# Remove the ".XX" part from Ensembl gene IDs
corrected_ensembl_genes <- sub("\\.\\d+", "", ensembl_genes)

# Perform GO enrichment analysis using corrected gene IDs
enrich_result <- enrichGO(gene = corrected_ensembl_genes, 
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENSEMBL",
                          ont = "BP")
# Open a PNG graphics device
png("enrichment_plot.png", width = 800, height = 600)

# Create the barplot
barplot(enrich_result, showCategory = 10)

# Close the graphics device
dev.off()

# Extract the list of GO terms from the enrich_result data frame
go_terms <- enrich_result$ID

# Display the extracted GO terms
print(go_terms)

file_path <- "go_ids.txt"

# Write the GO IDs to the text file
write.table(go_terms, file = file_path, row.names = FALSE, col.names = FALSE, quote = FALSE)

# Print a message to confirm the file has been saved
cat("GO IDs saved to:", file_path, "\n")


# SRsupport

setwd("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/diffIsoformUsage/data_SRsupport/")

# Prepare data

diff_isos.cytosolVSnucleus.filtered <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/diffIsoformUsage/data_SRsupport/known_diffIsoformUsage/diff_isos.cytosolVSnucleus.filtered.txt", header=FALSE)
colnames(diff_isos.cytosolVSnucleus.filtered) <- c("gene name", "isoform name", "pvalue", "sample1 isoform count", 
                                                   "sample2 isoform count", "sample1 alternative isoforms for gene count",
                                                   "sample2 alternative isoforms for gene count")
# Filter data for p-value < 0.05
filtered_data <- diff_isos.cytosolVSnucleus.filtered[diff_isos.cytosolVSnucleus.filtered$pvalue < 0.05, ]
# Filter out rows with NA values
filtered_data <- filtered_data[complete.cases(filtered_data), ]
colnames(filtered_data) <- c("gene name", "isoform name", "pvalue", "sample1 isoform count", 
                             "sample2 isoform count", "sample1 alternative isoforms for gene count",
                             "sample2 alternative isoforms for gene count")

# Compute adjusted p-values using the column 'pvalue'
filtered_data$adjusted_pvalue <- p.adjust(filtered_data$pvalue, method = "bonferroni")
filtered_data_adjustedPval <- filtered_data[filtered_data$adjusted_pvalue < 0.05, ]

# Specify the file path to save the TSV file in current directory
file_path <- "filtered_data_adjustedPval.tsv"

# Save the data with headers to a TSV file
write.table(filtered_data_adjustedPval, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# GO enrichment analysis:

# Load necessary packages
library(clusterProfiler)
library(org.Hs.eg.db)

# Create a data frame with gene symbols, p-values, and adjusted p-values

selected_columns <- c("gene name", "pvalue", "adjusted_pvalue")
gene_data <- filtered_data_adjustedPval[selected_columns]
colnames(gene_data)[1] <- "gene"

ensembl_genes <- gene_data$gene
# Remove the ".XX" part from Ensembl gene IDs
corrected_ensembl_genes <- sub("\\.\\d+", "", ensembl_genes)

# Perform GO enrichment analysis using corrected gene IDs
enrich_result <- enrichGO(gene = corrected_ensembl_genes, 
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENSEMBL",
                          ont = "BP")
# Open a PNG graphics device
png("enrichment_plot.png", width = 800, height = 600)

# Create the barplot
barplot(enrich_result, showCategory = 10)

# Close the graphics device
dev.off()

# Extract the list of GO terms from the enrich_result data frame
go_terms <- enrich_result$ID

go_terms <- enrich_result$Description

# Display the extracted GO terms
print(go_terms)

file_path <- "go_ids.txt"

# Write the GO IDs to the text file
write.table(go_terms, file = file_path, row.names = FALSE, col.names = FALSE, quote = FALSE)

# Print a message to confirm the file has been saved
cat("GO IDs saved to:", file_path, "\n")
