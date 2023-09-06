#path to save results
setwd("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_all")

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(gridExtra)

##Only TvsT CvsC and NvsN ROW Counts

#LRonly 
flair.quantify.counts_LRonly <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_LRonly/counts_matrix.tsv", header=TRUE)

# Create a list of column names to iterate through
column_names_LRonly <- colnames(flair.quantify.counts_LRonly)[-1]  # Exclude the 'ID' column

# Add a small constant to the data to avoid taking the log of zero
pseudocount <- 1e-10 + 1

# Loop through column combinations and create scatter plots
# Filter out non-numeric values and handle NA values
numeric_data <- flair.quantify.counts_LRonly %>%
  select(-ID) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  filter(if_any(everything(), is.numeric))

# Create a data frame with the specific column combinations you want to compare
comparisons <- data.frame(col1 = c("H0CRep1", "H0TRep1", "H0NRep1"),
                          col2 = c("H0CRep2", "H0TRep2", "H0NRep2"))

# Create a list to store the plots
plot_list <- list()

for (i in 1:nrow(comparisons)) {
  col1 <- comparisons$col1[i]
  col2 <- comparisons$col2[i]
  
  # Apply pseudocount to the data
  transformed_data <- numeric_data
  transformed_data[[col1]] <- transformed_data[[col1]] + pseudocount
  transformed_data[[col2]] <- transformed_data[[col2]] + pseudocount
  
  # Apply log2 transformation to the data
  transformed_data <- transformed_data %>%
    mutate(across(all_of(c(col1, col2)), ~ log2(. + pseudocount)))
  
  # Calculate the linear regression model
  regression_model <- lm(as.formula(paste(col2, "~", col1)), data = transformed_data)
  
  # Calculate the R-squared value
  r_squared <- summary(regression_model)$r.squared
  
  # Create a scatter plot with regression line
  plot <- ggplot(transformed_data, aes(x = .data[[col1]], y = .data[[col2]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste("Regression Plot: ", col1, " vs ", col2, "\nR-squared: ", round(r_squared, 4)),
         x = col1, y = col2) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 20),
      plot.title = element_text(size = 22),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
  # Add the plot to the list
  plot_list[[paste(col1, "_vs_", col2, ".png", sep = "")]] <- plot
}

# Combine all plots into a grid arrangement
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 1))

# Save as PDF
ggsave("combined_regression_plots_raw_LRonly.pdf", combined_plot, width = 12, height = 20)

# Save as png
ggsave("combined_regression_plots_raw_LRonly.png", combined_plot, width = 12, height = 20)

##############################
##### CLEAN R space before to RUN next steps
#############################
rm(list = ls())


##FILTERING
###########

flair.quantify.counts_LRonly <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_LRonly/counts_matrix.tsv", header=TRUE)

# Create a list of column names to iterate through
column_names_LRonly <- colnames(flair.quantify.counts_LRonly)[-1]  # Exclude the 'ID' column

# Add a small constant to the data to avoid taking the log of zero
pseudocount <- 1e-10 + 1

# Loop through column combinations and create scatter plots
# Filter out non-numeric values and handle NA values
numeric_data <- flair.quantify.counts_LRonly %>%
  select(-ID) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  filter(if_any(everything(), is.numeric))

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

# Convert the matrix to a data frame
abs_scaled_data_df <- as.data.frame(abs_scaled_data)

# Apply the data manipulation operations
numeric_data <- abs_scaled_data_df %>%
  select(H0CRep1_scaled, H0TRep1_scaled, H0NRep1_scaled, H0CRep2_scaled, H0TRep2_scaled, H0NRep2_scaled) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  filter(if_any(everything(), is.numeric))


# Create a data frame with the specific column combinations you want to compare
comparisons <- data.frame(col1 = c("H0CRep1_scaled", "H0TRep1_scaled", "H0NRep1_scaled"),
                          col2 = c("H0CRep2_scaled", "H0TRep2_scaled", "H0NRep2_scaled"))

# Create a list to store the plots
plot_list <- list()

for (i in 1:nrow(comparisons)) {
  col1 <- comparisons$col1[i]
  col2 <- comparisons$col2[i]
  
  # Apply pseudocount to the data
  transformed_data <- numeric_data
  transformed_data[[col1]] <- transformed_data[[col1]] + pseudocount
  transformed_data[[col2]] <- transformed_data[[col2]] + pseudocount
  
  # Apply log2 transformation to the data
  transformed_data <- transformed_data %>%
    mutate(across(all_of(c(col1, col2)), ~ log2(. + pseudocount)))
  
  # Calculate the linear regression model
  regression_model <- lm(as.formula(paste(col2, "~", col1)), data = transformed_data)
  
  # Calculate the R-squared value
  r_squared <- summary(regression_model)$r.squared
  
  # Create a scatter plot with regression line
  plot <- ggplot(transformed_data, aes(x = .data[[col1]], y = .data[[col2]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste("Regression Plot: ", col1, " vs ", col2, "\nR-squared: ", round(r_squared, 4)),
         x = col1, y = col2) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 20),
      plot.title = element_text(size = 22),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
  # Add the plot to the list
  plot_list[[paste(col1, "_vs_", col2, ".png", sep = "")]] <- plot
}

# Combine all plots into a grid arrangement
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 1))

# Save as PDF
ggsave("combined_regression_plots_scaled_LRonly.pdf", combined_plot, width = 12, height = 20)

# Save as png
ggsave("combined_regression_plots_scaled_LRonly.png", combined_plot, width = 12, height = 20)


## MERGE replicates

# Select the columns for each replicate group
h0c_columns <- c("H0CRep1_scaled", "H0CRep2_scaled")
h0t_columns <- c("H0TRep1_scaled", "H0TRep2_scaled")
h0n_columns <- c("H0NRep1_scaled", "H0NRep2_scaled")

# Calculate the mean for each replicate group
h0c_mean <- rowMeans(abs_scaled_data[, h0c_columns])
h0t_mean <- rowMeans(abs_scaled_data[, h0t_columns])
h0n_mean <- rowMeans(abs_scaled_data[, h0n_columns])

# Create a new data frame with the mean values and appropriate column names
scaled_mergedREP_data <- data.frame(H0C = h0c_mean, H0T = h0t_mean, H0N = h0n_mean)

# Display the merged data
head(scaled_mergedREP_data)

# Add the scaled columns to the filtered_data using cbind
filtered_scaled_data <- cbind(abs_scaled_data, scaled_mergedREP_data)


# correlation plots Merged vs Rep1
# Add a small constant to the data to avoid taking the log of zero
pseudocount <- 1e-10 + 1

# Loop through column combinations and create scatter plots
# Filter out non-numeric values and handle NA values
numeric_data <- filtered_scaled_data %>%
  select(c("H0CRep1_scaled", "H0TRep1_scaled", "H0NRep1_scaled", "H0C", "H0T", "H0N")) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  filter(if_any(everything(), is.numeric))


# Create a data frame with the specific column combinations you want to compare
comparisons <- data.frame(col1 = c("H0CRep1_scaled", "H0TRep1_scaled", "H0NRep1_scaled"),
                          col2 = c("H0C", "H0T", "H0N"))

# Create a list to store the plots
plot_list <- list()

for (i in 1:nrow(comparisons)) {
  col1 <- comparisons$col1[i]
  col2 <- comparisons$col2[i]
  
  # Apply pseudocount to the data
  transformed_data <- numeric_data
  transformed_data[[col1]] <- transformed_data[[col1]] + pseudocount
  transformed_data[[col2]] <- transformed_data[[col2]] + pseudocount
  
  # Apply log2 transformation to the data
  transformed_data <- transformed_data %>%
    mutate(across(all_of(c(col1, col2)), ~ log2(. + pseudocount)))
  
  # Calculate the linear regression model
  regression_model <- lm(as.formula(paste(col2, "~", col1)), data = transformed_data)
  
  # Calculate the R-squared value
  r_squared <- summary(regression_model)$r.squared
  
  # Create a scatter plot with regression line
  plot <- ggplot(transformed_data, aes(x = .data[[col1]], y = .data[[col2]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste("Regression Plot: ", col1, " vs ", col2, "\nR-squared: ", round(r_squared, 4)),
         x = col1, y = col2) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 20),
      plot.title = element_text(size = 22),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
  # Add the plot to the list
  plot_list[[paste(col1, "_vs_", col2, ".png", sep = "")]] <- plot
}

# Combine all plots into a grid arrangement
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 1))

# Save as PDF
ggsave("combined_regression_plots_scaled_MergedVSRep1_LRonly.pdf", combined_plot, width = 12, height = 20)

# Save as png
ggsave("combined_regression_plots_scaled_MergedVSRep1_LRonly.png", combined_plot, width = 12, height = 20)

# correlation plots Merged vs Rep2
# Add a small constant to the data to avoid taking the log of zero
pseudocount <- 1e-10 + 1

# Loop through column combinations and create scatter plots
# Filter out non-numeric values and handle NA values
numeric_data <- filtered_scaled_data %>%
  select(c("H0CRep2_scaled", "H0TRep2_scaled", "H0NRep2_scaled", "H0C", "H0T", "H0N")) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  filter(if_any(everything(), is.numeric))


# Create a data frame with the specific column combinations you want to compare
comparisons <- data.frame(col1 = c("H0CRep2_scaled", "H0TRep2_scaled", "H0NRep2_scaled"),
                          col2 = c("H0C", "H0T", "H0N"))

# Create a list to store the plots
plot_list <- list()

for (i in 1:nrow(comparisons)) {
  col1 <- comparisons$col1[i]
  col2 <- comparisons$col2[i]
  
  # Apply pseudocount to the data
  transformed_data <- numeric_data
  transformed_data[[col1]] <- transformed_data[[col1]] + pseudocount
  transformed_data[[col2]] <- transformed_data[[col2]] + pseudocount
  
  # Apply log2 transformation to the data
  transformed_data <- transformed_data %>%
    mutate(across(all_of(c(col1, col2)), ~ log2(. + pseudocount)))
  
  # Calculate the linear regression model
  regression_model <- lm(as.formula(paste(col2, "~", col1)), data = transformed_data)
  
  # Calculate the R-squared value
  r_squared <- summary(regression_model)$r.squared
  
  # Create a scatter plot with regression line
  plot <- ggplot(transformed_data, aes(x = .data[[col1]], y = .data[[col2]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste("Regression Plot: ", col1, " vs ", col2, "\nR-squared: ", round(r_squared, 4)),
         x = col1, y = col2) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 20),
      plot.title = element_text(size = 22),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
  # Add the plot to the list
  plot_list[[paste(col1, "_vs_", col2, ".png", sep = "")]] <- plot
}

# Combine all plots into a grid arrangement
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 1))


# Save as PDF
ggsave("combined_regression_plots_scaled_MergedVSRep2_LRonly.pdf", combined_plot, width = 12, height = 20)

# Save as png
ggsave("combined_regression_plots_scaled_MergedVSRep2_LRonly.png", combined_plot, width = 12, height = 20)




##############################
##### CLEAN R space before to RUN next steps
#############################
rm(list = ls())


#######################
##SRsupport

flair.quantify.counts_SRsupport <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_SRsupport/flair.quantify.counts.tsv", header = TRUE)

column_names_SRsupport <- colnames(flair.quantify.counts_SRsupport)[-1]  # Exclude the 'ID' column
# Create a list of column names to iterate through

# Add a small constant to the data to avoid taking the log of zero
pseudocount <- 1e-10 + 1

# Loop through column combinations and create scatter plots
# Filter out non-numeric values and handle NA values
numeric_data <- flair.quantify.counts_SRsupport %>%
  select(-ID) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  filter(if_any(everything(), is.numeric))

# Create a data frame with the specific column combinations you want to compare
comparisons <- data.frame(col1 = c("H0CRep1", "H0TRep1", "H0NRep1"),
                          col2 = c("H0CRep2", "H0TRep2", "H0NRep2"))

# Create a list to store the plots
plot_list <- list()

for (i in 1:nrow(comparisons)) {
  col1 <- comparisons$col1[i]
  col2 <- comparisons$col2[i]
  
  # Apply pseudocount to the data
  transformed_data <- numeric_data
  transformed_data[[col1]] <- transformed_data[[col1]] + pseudocount
  transformed_data[[col2]] <- transformed_data[[col2]] + pseudocount
  
  # Apply log2 transformation to the data
  transformed_data <- transformed_data %>%
    mutate(across(all_of(c(col1, col2)), ~ log2(. + pseudocount)))
  
  # Calculate the linear regression model
  regression_model <- lm(as.formula(paste(col2, "~", col1)), data = transformed_data)
  
  # Calculate the R-squared value
  r_squared <- summary(regression_model)$r.squared
  
  # Create a scatter plot with regression line
  plot <- ggplot(transformed_data, aes(x = .data[[col1]], y = .data[[col2]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste("Regression Plot: ", col1, " vs ", col2, "\nR-squared: ", round(r_squared, 4)),
         x = col1, y = col2) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 20),
      plot.title = element_text(size = 22),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
  # Add the plot to the list
  plot_list[[paste(col1, "_vs_", col2, ".png", sep = "")]] <- plot
}

# Combine all plots into a grid arrangement
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 1))

# Save as PDF
ggsave("combined_regression_plots_SRsupport.pdf", combined_plot, width = 12, height = 20)

# Save as png
ggsave("combined_regression_plots_SRsupport.png", combined_plot, width = 12, height = 20)


##############################
##### CLEAN R space before to RUN next steps
#############################
rm(list = ls())


##FILTERING
###########

flair.quantify.counts_SRsupport <- read.delim("~/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/MyMac/Documents/TFM/RESULTS/FLAIR_SRsupport/flair.quantify.counts.tsv", header = TRUE)

column_names_SRsupport <- colnames(flair.quantify.counts_SRsupport)[-1]  # Exclude the 'ID' column

# Add a small constant to the data to avoid taking the log of zero
pseudocount <- 1e-10 + 1

# Loop through column combinations and create scatter plots
# Filter out non-numeric values and handle NA values
numeric_data <- flair.quantify.counts_SRsupport %>%
  select(-ID) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  filter(if_any(everything(), is.numeric))

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

# Convert the matrix to a data frame
abs_scaled_data_df <- as.data.frame(abs_scaled_data)

# Apply the data manipulation operations
numeric_data <- abs_scaled_data_df %>%
  select(H0CRep1_scaled, H0TRep1_scaled, H0NRep1_scaled, H0CRep2_scaled, H0TRep2_scaled, H0NRep2_scaled) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  filter(if_any(everything(), is.numeric))


# Create a data frame with the specific column combinations you want to compare
comparisons <- data.frame(col1 = c("H0CRep1_scaled", "H0TRep1_scaled", "H0NRep1_scaled"),
                          col2 = c("H0CRep2_scaled", "H0TRep2_scaled", "H0NRep2_scaled"))

# Create a list to store the plots
plot_list <- list()

for (i in 1:nrow(comparisons)) {
  col1 <- comparisons$col1[i]
  col2 <- comparisons$col2[i]
  
  # Apply pseudocount to the data
  transformed_data <- numeric_data
  transformed_data[[col1]] <- transformed_data[[col1]] + pseudocount
  transformed_data[[col2]] <- transformed_data[[col2]] + pseudocount
  
  # Apply log2 transformation to the data
  transformed_data <- transformed_data %>%
    mutate(across(all_of(c(col1, col2)), ~ log2(. + pseudocount)))
  
  # Calculate the linear regression model
  regression_model <- lm(as.formula(paste(col2, "~", col1)), data = transformed_data)
  
  # Calculate the R-squared value
  r_squared <- summary(regression_model)$r.squared
  
  # Create a scatter plot with regression line
  plot <- ggplot(transformed_data, aes(x = .data[[col1]], y = .data[[col2]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste("Regression Plot: ", col1, " vs ", col2, "\nR-squared: ", round(r_squared, 4)),
         x = col1, y = col2) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 20),
      plot.title = element_text(size = 22),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
  # Add the plot to the list
  plot_list[[paste(col1, "_vs_", col2, ".png", sep = "")]] <- plot
}

# Combine all plots into a grid arrangement
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 1))

# Save as PDF
ggsave("combined_regression_plots_scaled_SRsupport.pdf", combined_plot, width = 12, height = 20)

# Save as png
ggsave("combined_regression_plots_scaled_SRsupport.png", combined_plot, width = 12, height = 20)

## MERGE replicates

# Select the columns for each replicate group
h0c_columns <- c("H0CRep1_scaled", "H0CRep2_scaled")
h0t_columns <- c("H0TRep1_scaled", "H0TRep2_scaled")
h0n_columns <- c("H0NRep1_scaled", "H0NRep2_scaled")

# Calculate the mean for each replicate group
h0c_mean <- rowMeans(abs_scaled_data[, h0c_columns])
h0t_mean <- rowMeans(abs_scaled_data[, h0t_columns])
h0n_mean <- rowMeans(abs_scaled_data[, h0n_columns])

# Create a new data frame with the mean values and appropriate column names
scaled_mergedREP_data <- data.frame(H0C = h0c_mean, H0T = h0t_mean, H0N = h0n_mean)

# Display the merged data
head(scaled_mergedREP_data)

# Add the scaled columns to the filtered_data using cbind
filtered_scaled_data <- cbind(abs_scaled_data, scaled_mergedREP_data)


# correlation plots Merged vs Rep1
# Add a small constant to the data to avoid taking the log of zero
pseudocount <- 1e-10 + 1

# Loop through column combinations and create scatter plots
# Filter out non-numeric values and handle NA values
numeric_data <- filtered_scaled_data %>%
  select(c("H0CRep1_scaled", "H0TRep1_scaled", "H0NRep1_scaled", "H0C", "H0T", "H0N")) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  filter(if_any(everything(), is.numeric))


# Create a data frame with the specific column combinations you want to compare
comparisons <- data.frame(col1 = c("H0CRep1_scaled", "H0TRep1_scaled", "H0NRep1_scaled"),
                          col2 = c("H0C", "H0T", "H0N"))

# Create a list to store the plots
plot_list <- list()

for (i in 1:nrow(comparisons)) {
  col1 <- comparisons$col1[i]
  col2 <- comparisons$col2[i]
  
  # Apply pseudocount to the data
  transformed_data <- numeric_data
  transformed_data[[col1]] <- transformed_data[[col1]] + pseudocount
  transformed_data[[col2]] <- transformed_data[[col2]] + pseudocount
  
  # Apply log2 transformation to the data
  transformed_data <- transformed_data %>%
    mutate(across(all_of(c(col1, col2)), ~ log2(. + pseudocount)))
  
  # Calculate the linear regression model
  regression_model <- lm(as.formula(paste(col2, "~", col1)), data = transformed_data)
  
  # Calculate the R-squared value
  r_squared <- summary(regression_model)$r.squared
  
  # Create a scatter plot with regression line
  plot <- ggplot(transformed_data, aes(x = .data[[col1]], y = .data[[col2]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste("Regression Plot: ", col1, " vs ", col2, "\nR-squared: ", round(r_squared, 4)),
         x = col1, y = col2) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 20),
      plot.title = element_text(size = 22),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
  # Add the plot to the list
  plot_list[[paste(col1, "_vs_", col2, ".png", sep = "")]] <- plot
}

# Combine all plots into a grid arrangement
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 1))

# Save as PDF
ggsave("combined_regression_plots_scaled_MergedVSRep1_SRsupport.pdf", combined_plot, width = 12, height = 20)

# Save as png
ggsave("combined_regression_plots_scaled_MergedVSRep1_SRsupport.png", combined_plot, width = 12, height = 20)

# correlation plots Merged vs Rep2
# Add a small constant to the data to avoid taking the log of zero
pseudocount <- 1e-10 + 1

# Loop through column combinations and create scatter plots
# Filter out non-numeric values and handle NA values
numeric_data <- filtered_scaled_data %>%
  select(c("H0CRep2_scaled", "H0TRep2_scaled", "H0NRep2_scaled", "H0C", "H0T", "H0N")) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  filter(if_any(everything(), is.numeric))


# Create a data frame with the specific column combinations you want to compare
comparisons <- data.frame(col1 = c("H0CRep2_scaled", "H0TRep2_scaled", "H0NRep2_scaled"),
                          col2 = c("H0C", "H0T", "H0N"))

# Create a list to store the plots
plot_list <- list()

for (i in 1:nrow(comparisons)) {
  col1 <- comparisons$col1[i]
  col2 <- comparisons$col2[i]
  
  # Apply pseudocount to the data
  transformed_data <- numeric_data
  transformed_data[[col1]] <- transformed_data[[col1]] + pseudocount
  transformed_data[[col2]] <- transformed_data[[col2]] + pseudocount
  
  # Apply log2 transformation to the data
  transformed_data <- transformed_data %>%
    mutate(across(all_of(c(col1, col2)), ~ log2(. + pseudocount)))
  
  # Calculate the linear regression model
  regression_model <- lm(as.formula(paste(col2, "~", col1)), data = transformed_data)
  
  # Calculate the R-squared value
  r_squared <- summary(regression_model)$r.squared
  
  # Create a scatter plot with regression line
  plot <- ggplot(transformed_data, aes(x = .data[[col1]], y = .data[[col2]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste("Regression Plot: ", col1, " vs ", col2, "\nR-squared: ", round(r_squared, 4)),
         x = col1, y = col2) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 20),
      plot.title = element_text(size = 22),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
  
  # Add the plot to the list
  plot_list[[paste(col1, "_vs_", col2, ".png", sep = "")]] <- plot
}

# Combine all plots into a grid arrangement
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 1))


# Save as PDF
ggsave("combined_regression_plots_scaled_MergedVSRep2_SRsupport.pdf", combined_plot, width = 12, height = 20)

# Save as png
ggsave("combined_regression_plots_scaled_MergedVSRep2_SRsupport.png", combined_plot, width = 12, height = 20)
