#LOAD libraries
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(UpSetR)

rm(list = ls())
# gencode_annotation <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/gencode.v43.primary_assembly.attrs.tsv")
# # Create the gencode list
# 
# #gencode_filtered <- gencode_annotation[gencode_annotation$transcriptType %in% c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "lncRNA", "protein_coding", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene"), ]
# gencode_filtered <- gencode_annotation[gencode_annotation$geneType %in% c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "lncRNA", "protein_coding", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene"), ]
# id_list_gencode <- paste(gencode_filtered$transcriptId, gencode_filtered$geneId, sep = "_")

### LRonly
# Where to save the ouput 
setwd("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_LRonly/UPSET_PLOTS")

# Loading data
counts_matrix <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_LRonly/counts_matrix.tsv", header = TRUE)
isoform_sizes <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_LRonly/flair.collapse.isoforms.sizes", header=TRUE, col.names=c("ID", "Length"))
gencode_annotation <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/gencode.v43.primary_assembly.attrs.tsv")


# Create the gencode list
#gencode_filtered <- gencode_annotation[gencode_annotation$transcriptType %in% c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "lncRNA", "protein_coding", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene"), ]
gencode_filtered <- gencode_annotation[gencode_annotation$geneType %in% c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "lncRNA", "protein_coding", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene"), ]
id_list_gencode <- paste(gencode_filtered$transcriptId, gencode_filtered$geneId, sep = "_")

counts_repl <- reshape2::melt(counts_matrix, 
                              measure.vars = c("H0CRep1", "H0TRep1", "H0NRep1", "H0CRep2", "H0TRep2", "H0NRep2"), variable.name = "Sample")

# Log2-transform the values
counts_repl$log2_value <- log2(counts_repl$value + 1)  # Adding 1 to avoid log(0)

# Create a new column to indicate the group
counts_repl$Group <- ifelse(counts_repl$Sample %in% c("H0CRep1", "H0CRep2"), "H0C",
                            ifelse(counts_repl$Sample %in% c("H0TRep1", "H0TRep2"), "H0T",
                                   ifelse(counts_repl$Sample %in% c("H0NRep1", "H0NRep2"), "H0N", "Other")))

# Create the box plot to explore the data
# Define the correct order of samples and associated colors
sample_order <- c("H0CRep1", "H0NRep1", "H0TRep1", "H0CRep2", "H0NRep2", "H0TRep2")
colors <- c("#d62728", "#1f77b4", "#2ca02c")

# Create the collapsed labels based on the Group column
collapsed_labels <- unique(counts_repl$Group)

# Create the box plot
p1 <- ggplot(counts_repl, aes(x = Sample, y = log2_value, fill = Group)) +
  geom_boxplot() +
  labs(x = "Samples", y = "Log2 Transformed Values", title = "Read Counts by Replicate") +
  scale_fill_manual(values = colors) +  # Assign custom colors
  theme_minimal()

# Adjust x-axis labels without collapsing
p1 + scale_x_discrete(breaks = sample_order, labels = sample_order)

# Save as PNG
png_filename <- "Read_Counts_by_Replicate_LRonly.png"
ggsave(png_filename, plot = p1, width = 8, height = 6, units = "in", dpi = 300)

# Save as PDF
pdf_filename <- "Read_Counts_by_Replicate_LRonly.pdf"
ggsave(pdf_filename, plot = p1, width = 8, height = 6)

#merge replicates from counts matrix

# Convert counts_matrix to a data frame
counts_matrix <- as.data.frame(counts_matrix)

merge_replicates <- counts_matrix <- counts_matrix %>%
  mutate(
    H0C = H0CRep1 + H0CRep2,
    H0T = H0TRep1 + H0TRep2,
    H0N = H0NRep1 + H0NRep2
  )

#add sizes column to df merge_replicates
# Inner join based on the 'ID' column

size_merged <- inner_join(merge_replicates, isoform_sizes, by = "ID")

# Assuming your data frame is named counts_matrix_Rep_Merged
# Select the relevant columns (IDs and the sample columns)
sample_columns <- c("ID", "H0C", "H0T", "H0N")
sample_data <- size_merged[sample_columns]

# Create a list of ID codes for each sample filter by 1 because we need at least to be present in the 2 replicates
id_list_h0c <- sample_data$ID[sample_data$H0C > 1]
id_list_h0t <- sample_data$ID[sample_data$H0T > 1]
id_list_h0n <- sample_data$ID[sample_data$H0N > 1]


# Combine the three lists into a list of sets
set_list <- list(
  H0C = id_list_h0c,
  H0T = id_list_h0t,
  H0N = id_list_h0n,
  GENCODE = id_list_gencode
)

upset(fromList(set_list), sets = c("H0C", "H0N", "H0T", "GENCODE"), order.by = "freq", 
      text.scale = c(1.2, 1.2, 1, 1, 1.5, 1.3), keep.order = TRUE)

# Create the UpSet plot
upset_plot <- upset(fromList(set_list), sets = c("H0C", "H0N", "H0T", "GENCODE"), order.by = "freq", 
                    text.scale = c(1.2, 1.2, 1, 1, 1.5, 1.3), keep.order = TRUE)

# Save as PNG
png("upset_plot_cellcompatments_LRonly_GENCODE.png", width = 800, height = 600)
print(upset_plot)
dev.off()

# Save as PDF
pdf("upset_plot_cellcompatments_LRonly_GENCODE.pdf", width = 8, height = 6)
print(upset_plot)
dev.off()

list.files()


##############################
##### CLEAN R space before to RUN next steps
#############################
rm(list = ls())


### SRsupport
# Where to save the ouput 
setwd("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_SRsupport/UPSET_PLOTS")

# Loading data
counts_matrix <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_SRsupport/flair.quantify.counts.tsv", header = TRUE)
isoform_sizes <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_SRsupport/flair.collapse.isoforms.sizes", header=TRUE, col.names=c("ID", "Length"))
gencode_annotation <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/gencode.v43.primary_assembly.attrs.tsv")


# Create the gencode list
#gencode_filtered <- gencode_annotation[gencode_annotation$transcriptType %in% c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "lncRNA", "protein_coding", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene"), ]
gencode_filtered <- gencode_annotation[gencode_annotation$geneType %in% c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "lncRNA", "protein_coding", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene"), ]
id_list_gencode <- paste(gencode_filtered$transcriptId, gencode_filtered$geneId, sep = "_")

counts_repl <- reshape2::melt(counts_matrix, 
                              measure.vars = c("H0CRep1", "H0TRep1", "H0NRep1", "H0CRep2", "H0TRep2", "H0NRep2"), variable.name = "Sample")

# Log2-transform the values
counts_repl$log2_value <- log2(counts_repl$value + 1)  # Adding 1 to avoid log(0)

# Create a new column to indicate the group
counts_repl$Group <- ifelse(counts_repl$Sample %in% c("H0CRep1", "H0CRep2"), "H0C",
                            ifelse(counts_repl$Sample %in% c("H0TRep1", "H0TRep2"), "H0T",
                                   ifelse(counts_repl$Sample %in% c("H0NRep1", "H0NRep2"), "H0N", "Other")))

# Create the box plot to explore the data
# Define the correct order of samples and associated colors
sample_order <- c("H0CRep1", "H0NRep1", "H0TRep1", "H0CRep2", "H0NRep2", "H0TRep2")
colors <- c("#d62728", "#1f77b4", "#2ca02c")

# Create the collapsed labels based on the Group column
collapsed_labels <- unique(counts_repl$Group)

# Create the box plot
p1 <- ggplot(counts_repl, aes(x = Sample, y = log2_value, fill = Group)) +
  geom_boxplot() +
  labs(x = "Samples", y = "Log2 Transformed Values", title = "Read Counts by Replicate") +
  scale_fill_manual(values = colors) +  # Assign custom colors
  theme_minimal()

# Adjust x-axis labels without collapsing
p1 + scale_x_discrete(breaks = sample_order, labels = sample_order)

# Save as PNG
png_filename <- "Read_Counts_by_Replicate_SRsupport.png"
ggsave(png_filename, plot = p1, width = 8, height = 6, units = "in", dpi = 300)

# Save as PDF
pdf_filename <- "Read_Counts_by_Replicate_SRsupport.pdf"
ggsave(pdf_filename, plot = p1, width = 8, height = 6)

#merge replicates from counts matrix

# Convert counts_matrix to a data frame
counts_matrix <- as.data.frame(counts_matrix)

merge_replicates <- counts_matrix <- counts_matrix %>%
  mutate(
    H0C = H0CRep1 + H0CRep2,
    H0T = H0TRep1 + H0TRep2,
    H0N = H0NRep1 + H0NRep2
  )

#add sizes column to df merge_replicates
# Inner join based on the 'ID' column

size_merged <- inner_join(merge_replicates, isoform_sizes, by = "ID")

# Assuming your data frame is named counts_matrix_Rep_Merged
# Select the relevant columns (IDs and the sample columns)
sample_columns <- c("ID", "H0C", "H0T", "H0N")
sample_data <- size_merged[sample_columns]

# Create a list of ID codes for each sample filter by 1 because we need at least to be present in the 2 replicates
id_list_h0c <- sample_data$ID[sample_data$H0C > 1]
id_list_h0t <- sample_data$ID[sample_data$H0T > 1]
id_list_h0n <- sample_data$ID[sample_data$H0N > 1]


# Combine the three lists into a list of sets
set_list <- list(
  H0C = id_list_h0c,
  H0T = id_list_h0t,
  H0N = id_list_h0n,
  GENCODE = id_list_gencode
)

upset(fromList(set_list), sets = c("H0C", "H0N", "H0T", "GENCODE"), order.by = "freq", 
      text.scale = c(1.2, 1.2, 1, 1, 1.5, 1.3), keep.order = TRUE)

# Create the UpSet plot
upset_plot <- upset(fromList(set_list), sets = c("H0C", "H0N", "H0T", "GENCODE"), order.by = "freq", 
                    text.scale = c(1.2, 1.2, 1, 1, 1.5, 1.3), keep.order = TRUE)

# Save as PNG
png("upset_plot_cellcompatments_SRsupport_GENCODE.png", width = 800, height = 600)
print(upset_plot)
dev.off()

# Save as PDF
pdf("upset_plot_cellcompatments_SRsupport_GENCODE.pdf", width = 8, height = 6)
print(upset_plot)
dev.off()

list.files()

##############################
##### CLEAN R space before to RUN next steps
#############################
rm(list = ls())

############################### FILETERED DATA #####################

### LRonly_filtered
# Where to save the ouput 
setwd("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_LRonly/UPSET_PLOTS")

# Loading data
counts_matrix <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_LRonly/known_sample_counts.tsv", header = TRUE)
isoform_sizes <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_LRonly/flair.collapse.isoforms.sizes", header=TRUE, col.names=c("ID", "Length"))
gencode_annotation <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/gencode.v43.primary_assembly.attrs.tsv")


# Create the gencode list
#gencode_filtered <- gencode_annotation[gencode_annotation$transcriptType %in% c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "lncRNA", "protein_coding", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene"), ]
gencode_filtered <- gencode_annotation[gencode_annotation$geneType %in% c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "lncRNA", "protein_coding", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene"), ]
id_list_gencode <- paste(gencode_filtered$transcriptId, gencode_filtered$geneId, sep = "_")

counts_repl <- reshape2::melt(counts_matrix, 
                              measure.vars = c("H0C", "H0T", "H0N"), variable.name = "Sample")

# Log2-transform the values
counts_repl$log2_value <- log2(counts_repl$value + 1)  # Adding 1 to avoid log(0)

# Create the box plot to explore the data
# Define the correct order of samples and associated colors
sample_order <- c("H0C", "H0N", "H0T")
colors <- c("#d62728", "#1f77b4", "#2ca02c")

# Create the collapsed labels based on the Group column
#collapsed_labels <- unique(counts_repl$Group)

# Create the box plot
p1 <- ggplot(counts_repl, aes(x = Sample, y = log2_value, fill = Sample)) +
  geom_boxplot() +
  labs(x = "Samples", y = "Log2 Transformed Values", title = "Read Counts by Replicate") +
  scale_fill_manual(values = colors, breaks = sample_order, labels = sample_order) +
  theme_minimal()

print(p1)

# Save as PNG
png_filename <- "Read_Counts_by_Replicate_LRonly_filtered_known.png"
ggsave(png_filename, plot = p1, width = 8, height = 6, units = "in", dpi = 300)

# Save as PDF
pdf_filename <- "Read_Counts_by_Replicate_LRonly_filtered_known.pdf"
ggsave(pdf_filename, plot = p1, width = 8, height = 6)


#add sizes column to df merge_replicates
# Inner join based on the 'ID' column

size_merged <- inner_join(counts_matrix, isoform_sizes, by = "ID")

# Assuming your data frame is named counts_matrix_Rep_Merged
# Select the relevant columns (IDs and the sample columns)
sample_columns <- c("ID", "H0C", "H0T", "H0N")
sample_data <- size_merged[sample_columns]

# Create a list of ID codes for each sample filter by 1 because we need at least to be present in the 2 replicates
id_list_h0c <- sample_data$ID[sample_data$H0C > 1]
id_list_h0t <- sample_data$ID[sample_data$H0T > 1]
id_list_h0n <- sample_data$ID[sample_data$H0N > 1]


# Combine the three lists into a list of sets
set_list <- list(
  H0C = id_list_h0c,
  H0T = id_list_h0t,
  H0N = id_list_h0n,
  GENCODE = id_list_gencode
)

upset(fromList(set_list), sets = c("H0C", "H0N", "H0T", "GENCODE"), order.by = "freq", 
      text.scale = c(1.2, 1.2, 1, 1, 1.5, 1.3), keep.order = TRUE)

# Create the UpSet plot
upset_plot <- upset(fromList(set_list), sets = c("H0C", "H0N", "H0T", "GENCODE"), order.by = "freq", 
                    text.scale = c(1.2, 1.2, 1, 1, 1.5, 1.3), keep.order = TRUE)

# Save as PNG
png("upset_plot_cellcompatments_LRonly_filter_known_GENCODE.png", width = 800, height = 600)
print(upset_plot)
dev.off()

# Save as PDF
pdf("upset_plot_cellcompatments_LRonly_filter_known_GENCODE.pdf", width = 8, height = 6)
print(upset_plot)
dev.off()

list.files()

##############################
##### CLEAN R space before to RUN next steps
#############################
rm(list = ls())

### SRsupport filtered
# Where to save the ouput 
setwd("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_SRsupport/UPSET_PLOTS")

# Loading data
counts_matrix <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_SRsupport/known_sample_counts.tsv", header = TRUE)
isoform_sizes <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_SRsupport/flair.collapse.isoforms.sizes", header=TRUE, col.names=c("ID", "Length"))
gencode_annotation <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/gencode.v43.primary_assembly.attrs.tsv")


# Create the gencode list
#gencode_filtered <- gencode_annotation[gencode_annotation$transcriptType %in% c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "lncRNA", "protein_coding", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene"), ]
gencode_filtered <- gencode_annotation[gencode_annotation$geneType %in% c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "lncRNA", "protein_coding", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene", "translated_processed_pseudogene", "translated_unprocessed_pseudogene"), ]
id_list_gencode <- paste(gencode_filtered$transcriptId, gencode_filtered$geneId, sep = "_")

counts_repl <- reshape2::melt(counts_matrix, 
                              measure.vars = c("H0C", "H0T", "H0N"), variable.name = "Sample")

# Log2-transform the values
counts_repl$log2_value <- log2(counts_repl$value + 1)  # Adding 1 to avoid log(0)

# Create the box plot to explore the data
# Define the correct order of samples and associated colors
sample_order <- c("H0C", "H0N", "H0T")
colors <- c("#d62728", "#1f77b4", "#2ca02c")

# Create the collapsed labels based on the Group column
#collapsed_labels <- unique(counts_repl$Group)

# Create the box plot
p1 <- ggplot(counts_repl, aes(x = Sample, y = log2_value, fill = Sample)) +
  geom_boxplot() +
  labs(x = "Samples", y = "Log2 Transformed Values", title = "Read Counts by Replicate") +
  scale_fill_manual(values = colors, breaks = sample_order, labels = sample_order) +
  theme_minimal()

print(p1)

# Save as PNG
png_filename <- "Read_Counts_by_Replicate_SRsupport_filtered_known.png"
ggsave(png_filename, plot = p1, width = 8, height = 6, units = "in", dpi = 300)

# Save as PDF
pdf_filename <- "Read_Counts_by_Replicate_SRsupport_filtered_known.pdf"
ggsave(pdf_filename, plot = p1, width = 8, height = 6)


#add sizes column to df merge_replicates
# Inner join based on the 'ID' column

size_merged <- inner_join(counts_matrix, isoform_sizes, by = "ID")

# Assuming your data frame is named counts_matrix_Rep_Merged
# Select the relevant columns (IDs and the sample columns)
sample_columns <- c("ID", "H0C", "H0T", "H0N")
sample_data <- size_merged[sample_columns]

# Create a list of ID codes for each sample filter by 1 because we need at least to be present in the 2 replicates
id_list_h0c <- sample_data$ID[sample_data$H0C > 1]
id_list_h0t <- sample_data$ID[sample_data$H0T > 1]
id_list_h0n <- sample_data$ID[sample_data$H0N > 1]


# Combine the three lists into a list of sets
set_list <- list(
  H0C = id_list_h0c,
  H0T = id_list_h0t,
  H0N = id_list_h0n,
  GENCODE = id_list_gencode
)

upset(fromList(set_list), sets = c("GENCODE", "H0T", "H0N", "H0C"), order.by = "freq", 
      text.scale = c(1.2, 1.2, 1, 1, 1.5, 1.3), keep.order = TRUE)

# Create the UpSet plot
upset_plot <- upset(fromList(set_list), sets = c("H0C", "H0N", "H0T", "GENCODE"), order.by = "freq", 
                    text.scale = c(1.2, 1.2, 1, 1, 1.5, 1.3), keep.order = TRUE)

# Save as PNG
png("upset_plot_cellcompatments_SRsupport_filter_known_GENCODE.png", width = 800, height = 600)
print(upset_plot)
dev.off()

# Save as PDF
pdf("upset_plot_cellcompatments_SRsupport_filter_known_GENCODE.pdf", width = 8, height = 6)
print(upset_plot)
dev.off()

list.files()


