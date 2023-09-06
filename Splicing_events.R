library(tidyr)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(UpSetR)


### BASH files=$(find /users/project/gencode_006070_no_backup/scarbonell/TFM/long_reads/FLAIR/all/data_LRonly/ -type f -name 'out*' -printf '"%f",') && echo "${files%,}"

file_names <- c("out.fishers.H0C.H0N.diffsplice.ir.events.quant.tsv",
                "out.fishers.H0C.H0N.diffsplice.es.events.quant.tsv",
                "out.fishers.H0C.H0T.diffsplice.ir.events.quant.tsv",
                "out.fishers.H0N.H0T.diffsplice.es.events.quant.tsv",
                "out.fishers.H0C.H0N.diffsplice.alt5.events.quant.tsv",
                "out.fishers.H0C.H0T.diffsplice.alt5.events.quant.tsv",
                "out.fishers.H0N.H0T.diffsplice.ir.events.quant.tsv",
                "out.fishers.H0C.H0T.diffsplice.es.events.quant.tsv",
                "out.fishers.H0C.H0N.diffsplice.alt3.events.quant.tsv",
                "out.fishers.H0N.H0T.diffsplice.alt5.events.quant.tsv",
                "out.fishers.H0N.H0T.diffsplice.alt3.events.quant.tsv",
                "out.fishers.H0C.H0T.diffsplice.alt3.events.quant.tsv")

## BASH: files=$(find /users/project/gencode_006070_no_backup/scarbonell/TFM/long_reads/FLAIR/all/data_LRonly/ -type f -name 'out.fishers*' -printf '"%f",')
# Loop through each filename and extract the label
# IFS=',' read -ra file_array <<< "$files"
# for file in "${file_array[@]}"; do
# label=$(echo "$file" | sed -E 's/out\.fishers\.([^.]+)\.([^.]+)\.flair\.diffsplice\.([^.]+)\..*$/\1.\2.\3/')
# echo -n "\"$label\","
# done


labels <- c("H0C.H0N.ir","H0C.H0N.es","H0C.H0T.ir","H0N.H0T.es",
            "H0C.H0N.alt5","H0C.H0T.alt5","H0N.H0T.ir","H0C.H0T.es",
            "H0C.H0N.alt3","H0N.H0T.alt5","H0N.H0T.alt3","H0C.H0T.alt3")

##########
### KNOW filter
#########

## LRonly

diffsplice_list <- lapply(file_names, function(file) {
  file_path <- file.path("/Users/scarbonell/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/DIffSplice/Known/known_FLAIR_input_results/data_LRonly", file)
  diff_isos <- read.delim(file_path, header = TRUE)
  diff_isos$compartment <- labels[file_names == file]
  diff_isos
})

setwd("/Users/scarbonell/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/DiffSplice/known/known_plots_files/data_LRonly/UPSETPLOTS")
list.files()

# Check column names in each data frame
for (i in seq_along(diffsplice_list)) {
  print(colnames(diffsplice_list[[i]]))
}

# Loop through the list and rename specific headers to "pvalue"
for (i in seq_along(diffsplice_list)) {
  colnames(diffsplice_list[[i]]) <- gsub(".*_pval$", "pvalue", colnames(diffsplice_list[[i]]))
}

# Check the modified column names in each data frame
for (i in seq_along(diffsplice_list)) {
  print(colnames(diffsplice_list[[i]]))
}


merged_diffsplice <- do.call(rbind, diffsplice_list)


merged_diffsplice <- merged_diffsplice %>%
  mutate(label = gsub("^(H0\\w+)\\.(H0\\w+)\\..*$", "\\1vs\\2", compartment))

rownames(merged_diffsplice) <- sub("^.*\\.", "", rownames(merged_diffsplice))

# Calculate the corrected p-values using FDR method
merged_diffsplice$corrected_pvalue <- p.adjust(merged_diffsplice$pvalue, method = "fdr")
# threshold for adjusted p-values will be also 0.05, which corresponds to a 5% FDR

# Print the updated dataset
head(merged_diffsplice)

# Split the compartment column by dots and select the last part
merged_diffsplice$event <- sapply(strsplit(merged_diffsplice$compartment, "\\."), function(x) tail(x, 1))

# Print the modified dataframe
head(merged_diffsplice)


########## pvalue ####

# Filter data with p-value < 0.05
significant_data <- merged_diffsplice[merged_diffsplice$pvalue < 0.05, ]

# Calculate the count of events for each label
label_counts <- table(significant_data$label)

# Create a plot
p1 <- ggplot(significant_data, aes(x = label, fill = compartment, weight = label_counts[label])) +
  geom_bar(position = "dodge") +
  labs(x = "Labels", y = "Count of Events", title = "p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 4000000000))

# Print the plot
print(p1)

dev.off()


# Filter data with corrected p-value < 0.05 (5% FDR)
significant_data2 <- merged_diffsplice[merged_diffsplice$corrected_pvalue < 0.05, ]

# Calculate the count of events for each label
label_counts <- table(significant_data2$label)

# Create a plot
p2 <- ggplot(significant_data2, aes(x = label, fill = compartment, weight = label_counts[label])) +
  geom_bar(position = "dodge") +
  labs(x = "Labels", y = "Count of Events", title = "Ajusted p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 4000000000))

# Print the plot
print(p2)

# Create the grid of plots
grid_plots <- grid.arrange(p1, p2, ncol = 2)

#### with ajusted p-values 

# Split the data into a list by label
data_by_label <- split(significant_data2, significant_data2$label)

# Create a list to store isoform_ids for each label
isoform_ids_list <- lapply(data_by_label, function(df) df$isoform_ids)


# Combine isoform_ids into a single character vector
isoform_ids_combined <- unlist(isoform_ids_list)

# Prepare the data for the UpSet plot
upset_data <- data.frame(
  label = rep(names(isoform_ids_list), lengths(isoform_ids_list)),
  isoform_ids = isoform_ids_combined,  # Use the combined character vector
  stringsAsFactors = FALSE
)

# Create a named list with isoform_ids for each label
isoform_ids_named_list <- list()
for (label in unique(upset_data$label)) {
  isoform_ids_named_list[[label]] <- upset_data$isoform_ids[upset_data$label == label]
}

# Create the UpSet plot
upset_plot <- upset(fromList(isoform_ids_named_list), order.by = "freq")
print(upset_plot)


# Save as PNG
png("upset_plot_cellcompatments_LRonly_SPLICING.png", width = 800, height = 600)
print(upset_plot)
dev.off()

# Save as PDF
pdf("upset_plot_cellcompatments_LRonly_SPLICING", width = 8, height = 6)
print(upset_plot)
dev.off()
############

########################################
rm(list = ls())

## Reaload file_manes and labels

## SRsupport

diffsplice_list <- lapply(file_names, function(file) {
  file_path <- file.path("/Users/scarbonell/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/DIffSplice/Known/known_FLAIR_input_results/data_SRsupport", file)
  diff_isos <- read.delim(file_path, header = TRUE)
  diff_isos$compartment <- labels[file_names == file]
  diff_isos
})

setwd("/Users/scarbonell/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/DiffSplice/known/known_plots_files/data_SRsupport/UPSETPLOTS")
list.files()

# Check column names in each data frame
for (i in seq_along(diffsplice_list)) {
  print(colnames(diffsplice_list[[i]]))
}

# Loop through the list and rename specific headers to "pvalue"
for (i in seq_along(diffsplice_list)) {
  colnames(diffsplice_list[[i]]) <- gsub(".*_pval$", "pvalue", colnames(diffsplice_list[[i]]))
}

# Check the modified column names in each data frame
for (i in seq_along(diffsplice_list)) {
  print(colnames(diffsplice_list[[i]]))
}


merged_diffsplice <- do.call(rbind, diffsplice_list)


merged_diffsplice <- merged_diffsplice %>%
  mutate(label = gsub("^(H0\\w+)\\.(H0\\w+)\\..*$", "\\1vs\\2", compartment))

rownames(merged_diffsplice) <- sub("^.*\\.", "", rownames(merged_diffsplice))

# Calculate the corrected p-values using FDR method
merged_diffsplice$corrected_pvalue <- p.adjust(merged_diffsplice$pvalue, method = "fdr")
# threshold for adjusted p-values will be also 0.05, which corresponds to a 5% FDR

# Split the compartment column by dots and select the last part
merged_diffsplice$event <- sapply(strsplit(merged_diffsplice$compartment, "\\."), function(x) tail(x, 1))


# Print the updated dataset
head(merged_diffsplice)


########## pvalue ####

# Filter data with p-value < 0.05
significant_data <- merged_diffsplice[merged_diffsplice$pvalue < 0.05, ]

# Calculate the count of events for each label
label_counts <- table(significant_data$label)

# Define a custom color palette with 4 distinct colors
custom_palette <- rep(c("#588c7e", "#f2e394", "#f2ae72", "#d96459"), 3)  # Repeat 4 colors 3 times

p1 <- ggplot(significant_data, aes(x = label, fill = event, weight = label_counts[label])) +
  geom_bar(position = "dodge") +
  labs(x = NULL, y = "Count of Events", title = "p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 6000000000)) +
  scale_fill_manual(values = custom_palette)  # Set custom color palette

print(p1)


# Filter data with corrected p-value < 0.05 (5% FDR)
significant_data2 <- merged_diffsplice[merged_diffsplice$corrected_pvalue < 0.05, ]

# Calculate the count of events for each label
label_counts <- table(significant_data2$label)


# Define a custom color palette with 4 distinct colors
custom_palette <- rep(c("#588c7e", "#f2e394", "#f2ae72", "#d96459"), 3)  # Repeat 4 colors 3 times

p2 <- ggplot(significant_data2, aes(x = label, fill = event, weight = label_counts[label])) +
  geom_bar(position = "dodge") +
  labs(x = NULL, y = "Count of Events", title = "Ajusted p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 6000000000)) +
  scale_fill_manual(values = custom_palette)  # Set custom color palette

print(p2)


# Create the grid of plots
grid_plots <- grid.arrange(p1, p2, ncol = 2)

# # Save as png
# png("Splicing_pvalues_SRsupport_filter_known_GENCODE.png", width = 8, height = 6)
# print(grid_plots)
# dev.off()

# Save as PDF
pdf("Splicing_pvalues_SRsupport_filter_known_GENCODE.pdf", width = 8, height = 6)
grid.arrange(p1, p2, ncol = 2)
dev.off()

#### with ajusted p-values 


# Split the data into a list by label
data_by_label <- split(significant_data2, significant_data2$label)

# Create a list to store isoform_ids for each label
isoform_ids_list <- lapply(data_by_label, function(df) df$isoform_ids)


# Combine isoform_ids into a single character vector
isoform_ids_combined <- unlist(isoform_ids_list)

# Prepare the data for the UpSet plot
upset_data <- data.frame(
  label = rep(names(isoform_ids_list), lengths(isoform_ids_list)),
  isoform_ids = isoform_ids_combined,  # Use the combined character vector
  stringsAsFactors = FALSE
)

# Create a named list with isoform_ids for each label
isoform_ids_named_list <- list()
for (label in unique(upset_data$label)) {
  isoform_ids_named_list[[label]] <- upset_data$isoform_ids[upset_data$label == label]
}

# Create the UpSet plot
upset_plot <- upset(fromList(isoform_ids_named_list), order.by = "freq")
print(upset_plot)

# Save as PNG
png("upset_plot_cellcompatments_SRsupport_SPLICING.png", width = 800, height = 600)
print(upset_plot)
dev.off()

# Save as PDF
pdf("upset_plot_cellcompatments_SRsupport_SPLICING.pdf", width = 8, height = 6)
print(upset_plot)
dev.off()
############

