One Million MEPs (OMM): A Database of Human Muscle Responses to Brain Stimulation
Project Overview
The One Million MEPs (OMM) project aims to build a comprehensive database of human motor-evoked potentials (MEPs) recorded via Transcranial Magnetic Stimulation (TMS). TMS is a non-invasive brain stimulation technique used to study the brain's role in motor control by generating electromagnetic pulses that stimulate neurons in the motor cortex. These stimulations elicit small muscle twitches, termed motor-evoked potentials (MEPs), which are recorded through surface electrodes.

This project collects, processes, and analyzes a large-scale dataset of MEPs, totaling around 1 million MEP responses. The goal is to standardize the MEP data into a common format and to provide tools for statistical analysis and data visualization to the neuroscience research community.

The One Million MEPs (OMM) database contributes to the NIBS-BIDS (Non-invasive Brain Stimulation Brain Imaging Data Structure) standard, enabling easier data sharing and reproducibility in research.

Project Goals
Collating MEP Data: Process around 1 million MEP responses from various experiments into an annotated database.
Standardization: Convert the MEP data into a standard JSON format compatible with the NIBS-BIDS project for open sharing.
Statistical Analysis: Provide R scripts for the statistical analysis of MEP data, including t-tests, ANOVA, and other methods to study MEP response patterns.
Visualization: Develop graphical representations of MEP distributions and comparisons across different experimental conditions.
Data Export: Enable exporting the processed MEP data to JSON for easy sharing and standardization.

# preprocess_data.R
# Import libraries
library(dplyr)
library(tidyr)

# Load raw data
load_mep_data <- function(file_path) {
  data <- read.csv(file_path)
  return(data)
}

# Preprocess and clean the data
preprocess_mep_data <- function(data) {
  # Handle missing values (replace missing MEP responses with NA)
  data_clean <- data %>% drop_na()

  # Normalize the MEP amplitude data (assuming 'amplitude' is a column)
  data_clean <- data_clean %>%
    mutate(amplitude_normalized = (amplitude - min(amplitude)) / (max(amplitude) - min(amplitude)))

  # Filter out noisy or erroneous data (e.g., amplitude below a noise threshold)
  data_clean <- data_clean %>%
    filter(amplitude_normalized > 0.1)  # Example threshold
  
  return(data_clean)
}

# Save the cleaned data to CSV
save_clean_data <- function(data, output_path) {
  write.csv(data, file = output_path, row.names = FALSE)
}

# Main function to load, preprocess, and save the cleaned data
main_preprocessing <- function() {
  # Load raw data
  raw_data <- load_mep_data("../data/raw_meps.csv")
  
  # Preprocess the data
  cleaned_data <- preprocess_mep_data(raw_data)
  
  # Save the cleaned data
  save_clean_data(cleaned_data, "../data/cleaned_meps.csv")
  
  print("Data preprocessing completed!")
}

# Run preprocessing
main_preprocessing()
# analysis.R
# Import libraries
library(dplyr)

# Perform t-test and ANOVA
analyze_mep_data <- function(data) {
  # T-test: comparing two conditions (e.g., high vs low stimulation intensity)
  t_test_result <- t.test(amplitude_normalized ~ condition, data = data)
  print("T-Test Result:")
  print(t_test_result)
  
  # ANOVA: comparing multiple conditions
  aov_result <- aov(amplitude_normalized ~ condition, data = data)
  print("ANOVA Result:")
  print(summary(aov_result))
  
  return(list(t_test = t_test_result, aov = aov_result))
}

# Save statistical analysis results to file
save_analysis_results <- function(t_test_result, aov_result, output_file) {
  write(paste("T-Test Result:\n", capture.output(t_test_result), "\n\nANOVA Result:\n", capture.output(summary(aov_result))),
        file = output_file)
}

# Main function to perform analysis
main_analysis <- function() {
  # Load cleaned data
  cleaned_data <- read.csv("../data/cleaned_meps.csv")
  
  # Perform statistical analysis
  results <- analyze_mep_data(cleaned_data)
  
  # Save results to a file
  save_analysis_results(results$t_test, results$aov, "../data/analysis_results.txt")
  
  print("Statistical analysis completed!")
}

# Run analysis
main_analysis()
# visualization.R
# Import libraries
library(ggplot2)
library(dplyr)

# Plot MEP amplitude distribution
plot_mep_distribution <- function(data) {
  plot <- ggplot(data, aes(x = amplitude_normalized, fill = condition)) +
    geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
    labs(title = "Distribution of MEP Amplitudes", x = "Normalized MEP Amplitude", y = "Frequency") +
    theme_minimal()
  
  return(plot)
}

# Plot MEP amplitude comparison between conditions
plot_group_comparison <- function(data) {
  plot <- ggplot(data, aes(x = condition, y = amplitude_normalized, fill = condition)) +
    geom_boxplot() +
    labs(title = "Comparison of MEP Amplitudes Between Conditions", x = "Condition", y = "Normalized MEP Amplitude") +
    theme_minimal()
  
  return(plot)
}

# Save plots to files
save_plot <- function(plot, output_path) {
  ggsave(output_path, plot = plot, width = 8, height = 6)
}

# Main function to generate and save visualizations
main_visualization <- function() {
  # Load cleaned data
  cleaned_data <- read.csv("../data/cleaned_meps.csv")
  
  # Generate and save MEP amplitude distribution plot
  dist_plot <- plot_mep_distribution(cleaned_data)
  save_plot(dist_plot, "../figures/mep_distribution.png")
  
  # Generate and save MEP amplitude comparison plot
  comparison_plot <- plot_group_comparison(cleaned_data)
  save_plot(comparison_plot, "../figures/mep_group_comparison.png")
  
  print("Data visualizations created and saved!")
}

# Run visualization
main_visualization()
# export_json.R
# Import libraries
library(jsonlite)
library(dplyr)

# Convert MEP data to JSON format
export_mep_to_json <- function(data, output_path) {
  json_data <- toJSON(data, pretty = TRUE, auto_unbox = TRUE)
  
  # Save JSON to file
  write(json_data, file = output_path)
}

# Main function to load cleaned data and export it to JSON
main_export_json <- function() {
  # Load cleaned data
  cleaned_data <- read.csv("../data/cleaned_meps.csv")
  
  # Export data to JSON format
  export_mep_to_json(cleaned_data, "../data/mep_data.json")
  
  print("MEP data exported to JSON format!")
}

# Run JSON export
main_export_json()
