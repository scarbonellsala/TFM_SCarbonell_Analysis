# Load libraries
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)

setwd("/Users/scarbonell/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/R_scripts/TFM_SCarbonell/figures")

# Center and scale the data
# Centering data means that the average of a variable is subtracted from the data. 
# Scaling data means that the standard deviation of a variable is divided out of the data.

##########
# LR only
##########

# Load data
flair.quantify.counts_LRonly <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_LRonly/counts_matrix.tsv", header=TRUE)
isoform_sizes <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_LRonly/flair.collapse.isoforms.sizes", header=TRUE, col.names=c("ID", "Length"))
# Convert flair.quantify.counts_LRonly to a data frame if it's not already
flair.quantify.counts_LRonly <- as.data.frame(flair.quantify.counts_LRonly)

# Perform the inner join to add a column containing isoform sizes
data_LRonly_with_length <- inner_join(flair.quantify.counts_LRonly, isoform_sizes, by = "ID")

head(data_LRonly_with_length)

# Select the numeric columns from column 2 to 7
numeric_data <- data_LRonly_with_length[, 2:7]

# Scale the numeric data
scaled_data <- scale(numeric_data, center = TRUE, scale = TRUE)

# Check for values equal to 0 in scaled data
zero_values <- apply(scaled_data == 0, 2, all)
# zero_values will be a logical vector indicating whether each column contains only zeros
# Print columns with all values equal to 0
colnames(scaled_data)[zero_values]

# Check for NA values in scaled data
na_values <- apply(is.na(scaled_data), 2, any)
# na_values will be a logical vector indicating whether each column contains any NA values
# Print columns with NA values
colnames(scaled_data)[na_values]

# Rename the scaled columns and make absolute value
colnames(scaled_data) <- paste0(colnames(numeric_data), "_scaled")
abs_scaled_data <- abs(scaled_data)


head(abs_scaled_data)

counts_for_pca <- t(abs_scaled_data)

counts_for_pca <- counts_for_pca[, apply(counts_for_pca, 2, var) != 0]
pca_result <- prcomp(counts_for_pca, scale. = FALSE)
data = data.frame(pca_result$x, Sample = str_sub(rownames(counts_for_pca), 3, 3))

# Define a custom color palette with 3 distinct colors repeated 2 times
custom_colors <- rep(c("red", "blue", "#00cc00"), 2)

pca_plot1 <- ggplot(data, aes(x = PC1, y = PC2, col = Sample)) +
  geom_point(size = 5) +
  labs(title = "PCA LR only centered and scaled data") +
  scale_color_manual(values = custom_colors)  # Set custom color palette

# Save as png
png("PCA_LRonly_center_scaled.png", width = 8, height = 6)
print(pca_plot1)
dev.off()

# Save as PDF
pdf("PCA_LRonly_center_scaled.pdf", width = 8, height = 6)
print(pca_plot1)
dev.off()

############
# SR support
############

# Load data
flair.quantify.counts_SRsupport <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_SRsupport/flair.quantify.counts.tsv", header=TRUE)
isoform_sizes <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_SRsupport/flair.collapse.isoforms.sizes", header=TRUE, col.names=c("ID", "Length"))
# Convert flair.quantify.counts_SRsupport to a data frame if it's not already
flair.quantify.counts_SRsupport <- as.data.frame(flair.quantify.counts_SRsupport)

# Perform the inner join to add a column containing isoform sizes
data_SRsupport_with_length <- inner_join(flair.quantify.counts_SRsupport, isoform_sizes, by = "ID")

head(data_SRsupport_with_length)

# Select the numeric columns from column 2 to 7
numeric_data <- data_SRsupport_with_length[, 2:7]

# Scale the numeric data
scaled_data <- scale(numeric_data, center = TRUE, scale = TRUE)

# Check for values equal to 0 in scaled data
zero_values <- apply(scaled_data == 0, 2, all)
# zero_values will be a logical vector indicating whether each column contains only zeros
# Print columns with all values equal to 0
colnames(scaled_data)[zero_values]

# Check for NA values in scaled data
na_values <- apply(is.na(scaled_data), 2, any)
# na_values will be a logical vector indicating whether each column contains any NA values
# Print columns with NA values
colnames(scaled_data)[na_values]

# Rename the scaled columns and make absolute value
colnames(scaled_data) <- paste0(colnames(numeric_data), "_scaled")
abs_scaled_data <- abs(scaled_data)


head(abs_scaled_data)

counts_for_pca <- t(abs_scaled_data)

counts_for_pca <- counts_for_pca[, apply(counts_for_pca, 2, var) != 0]
pca_result <- prcomp(counts_for_pca, scale. = FALSE)
data = data.frame(pca_result$x, Sample = str_sub(rownames(counts_for_pca), 3, 3))

# Define a custom color palette with 3 distinct colors repeated 2 times
custom_colors <- rep(c("red", "blue", "#00cc00"), 2)

pca_plot2 <- ggplot(data, aes(x = PC1, y = PC2, col = Sample)) +
  geom_point(size = 5) +
  labs(title = "PCA SR support centered and scaled data") +
  scale_color_manual(values = custom_colors)  # Set custom color palette


# Save as png
png("PCA_SRsupported_center_scaled.png", width = 8, height = 6)
print(pca_plot2)
dev.off()

# Save as PDF
pdf("PCA_SRsupported_center_scaled.pdf", width = 8, height = 6)
print(pca_plot2)
dev.off()
