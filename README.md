Profiling the Peripheral Blood Immunophenotype in Pancreatic Adenocarcinoma

Pancreatic ductal adenocarcinoma (PDAC) is one of the deadliest cancers, with a dismal 5-year survival rate of just 9%, primarily due to late-stage diagnoses and limited treatment options. Despite substantial progress in oncology, there is still a critical need for new insights into PDAC's pathogenesis and identification of novel therapeutic targets.

This project aims to profile the immune signatures found in the peripheral blood of PDAC patients and compare them to healthy individuals and those with benign pancreatic diseases. Such an analysis could uncover immune biomarkers for early diagnosis, prognostic stratification, and personalized treatment approaches.

While previous research has focused on immune profiles in PDAC, a comprehensive systemic analysis of immune cell composition, particularly from peripheral blood and in relation to clinicopathological features, is missing in the current literature.

Immune phenotyping of peripheral blood from PDAC patients could reveal distinct immune signatures that are associated with disease progression, prognosis, and potential therapeutic interventions.

The study will leverage Cytometry by Time of Flight (CyTOF) data to provide a single-cell resolution analysis of 30 immune cell marker proteins in peripheral blood samples.

Cohorts:
PDAC patients (n=67)
Post-metastatic PDAC (n=31)
Benign pancreatic disease (n=32)
Healthy donors (n=32)
For each sample, detailed clinicopathological metadata is available, along with additional bulk exome and transcriptome data for the surgery-eligible PDAC patients.

Project Aims
Immune Phenotype Profiling: Use CyTOF data to comprehensively profile the immune cell composition of peripheral blood samples in PDAC patients and compare these profiles to healthy individuals and those with benign disease.

Biomarker Discovery: Identify immune signatures that could serve as potential biomarkers for early PDAC diagnosis and prognostic stratification.

Therapeutic Insights: Explore immune phenotypes that may inform personalized therapeutic interventions based on immune system responses in PDAC patients.

This project will utilize R as the primary computational environment. Key steps will follow established CyTOF analysis workflows, including data preprocessing, clustering, and visualization. Below is a summary of the workflow:

R Packages:
flowCore: For reading and handling flow cytometry data.
CATALYST: For preprocessing, normalization, and data transformation.
cytoTree: For clustering, unsupervised analysis, and marker enrichment modeling.
ggplot2/tSNE/UMAP: For dimensionality reduction and visualization of complex data.
survival: For survival analysis and Cox proportional hazard modeling.

Data Preprocessing: Data will be preprocessed using CyTOF-specific R packages. This will involve quality checks, batch correction, normalization, and transformation.

Clustering and Annotation: Unsupervised clustering will be performed to identify distinct immune cell types in the samples. Clusters will be annotated based on enriched immune markers.

Dimensionality Reduction: Visualization of high-dimensional data will be achieved using tSNE and UMAP to map immune cell types and their associations with disease status.

Statistical Analysis: Statistical methods will be employed to evaluate correlations between immune cell compositions, clinicopathological data, and patient survival using multivariate models, including Cox proportional hazards.

This project will apply advanced bioinformatics techniques to CyTOF data for identifying potential biomarkers and immune signatures relevant to PDAC. Through systematic profiling of immune phenotypes in patient blood samples, this study seeks to inform new avenues for early diagnosis, prognosis, and treatment in pancreatic cancer.


fcs_path <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/DATA/"

library(flowCore)
library(ggplot2)
library(dplyr)

# Read the FCS file
fcs_data <- read.FCS(fcs_file_path)

# Display a summary of the FCS data
summary(fcs_data)

# Define the meta data path
meta_path <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/META"

# Load metadata files
immune_panel <- read_csv(file.path(meta_path, "Immune_panel.csv"))
filebatch <- read_csv(file.path(meta_path, "filebatch.csv"))
sample_meta <- read_csv(file.path(meta_path, "sample_meta.csv"))

# List all CSV files in the DEBARCODING_KEYS subfolder
debarcoding_keys_path <- file.path(meta_path, "DEBARCODING_KEYS")
debarcoding_key_files <- list.files(debarcoding_keys_path, pattern = "*.csv", full.names = TRUE)

# Load and preview each debarcoding key file separately
debarcoding_keys <- lapply(debarcoding_key_files, function(file) {
  key <- read_csv(file)
  print(paste("Preview of", file))
  print(head(key))
  return(key)
})

# Preview other metadata files
print("Preview of immune_panel")
head(immune_panel)

print("Preview of filebatch")
head(filebatch)

print("Preview of sample_meta")
head(sample_meta)

# Set the directory where cleaned FCS files will be saved
clean_dir <- file.path(getwd(), "Cleaned")
if (!dir.exists(clean_dir)) dir.create(clean_dir)

# List FCS files
fcs_files <- list.files("/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/DATA/", pattern = "\\.fcs$", full.names = TRUE)

# Loop through each file
for (file in fcs_files) {
  
  # Read FCS file
  ff <- flowCore::read.FCS(filename = file, transformation = FALSE)

  # Clean flow rate using flowAI package
  ff_cl <- clean_flow_rate(flow_frame = ff, 
                           out_dir = clean_dir, 
                           to_plot = FALSE,
                           data_type = "MC",
                           a = 0.005,
                           pen_FS = 10000, 
                           remove = "all") 

  # Clean signal using flowCut package
  ff_cl_s <- clean_signal(flow_frame = ff_cl,
                          to_plot = "All",
                          out_dir = clean_dir,
                          Segment = 1000,
                          arcsine_transform = TRUE,
                          data_type = "MC",
                          non_used_bead_ch = "140",
                          MaxPercCut = 0.2) 

  # Write cleaned FCS files
  flowCore::write.FCS(ff_cl_s,
                      file = file.path(clean_dir, gsub("_Normalized","_cleaned", basename(file)))) 
}

# Define the source and destination directories
source_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/"
dest_dir <- file.path(source_dir, "Cleaned")

# Create the destination directory if it doesn't exist
if (!dir.exists(dest_dir)) {
  dir.create(dest_dir, recursive = TRUE)
}

# List all cleaned FCS files in the source directory
files_to_move <- list.files(source_dir, pattern = "_cleaned\\.fcs$", full.names = TRUE)

# Move files to the destination directory
for (file in files_to_move) {
  file.rename(file, file.path(dest_dir, basename(file)))
}

```


```{r}
file_list <- list.files("/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/CLEANED/")
```


```{r}
# Load necessary libraries
library(flowCore)
library(data.table)

# Define the file paths
cleaned_fcs_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/Cleaned/"
debarcoding_key_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/META/DEBARCODING_KEYS/"

# Function to load and visualize a cleaned .fcs file
load_and_visualize_fcs <- function(fcs_file) {
  # Load the .fcs file
  fcs_data <- read.FCS(fcs_file)
  
  # Print summary of the .fcs file
  print(summary(fcs_data))
  
  # Visualize the first few rows of the .fcs data
  exprs_data <- exprs(fcs_data)
  print(head(exprs_data))
}

# Function to load and visualize a debarcoding key
load_and_visualize_debarcoding_key <- function(debarcoding_key_file) {
  # Load the debarcoding key
  debarcoding_key <- fread(debarcoding_key_file)
  
  # Print summary of the debarcoding key
  print(summary(debarcoding_key))
  
  # Visualize the first few rows of the debarcoding key
  print(head(debarcoding_key))
}

# Example usage for the first file
# File name without extensions
base_name <- "T01_E17_M46_HD266_01"

# Construct full file paths
fcs_file <- paste0(cleaned_fcs_dir, base_name, "_cleaned.fcs")
debarcoding_key_file <- paste0(debarcoding_key_dir, base_name, ".csv")

# Load and visualize the cleaned .fcs file
load_and_visualize_fcs(fcs_file)

# Load and visualize the debarcoding key
load_and_visualize_debarcoding_key(debarcoding_key_file)

# Load necessary libraries
library(flowCore)
library(CATALYST)

# Define directories
cleaned_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/CLEANED/"
meta_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/META/"
debarcoding_key_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/META/DEBARCODING_KEYS/"
debarcoded_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/DEBARCODED/"

# Create debarcoded directory if it doesn't exist
if (!dir.exists(debarcoded_dir)) dir.create(debarcoded_dir)

# List all .fcs files in the cleaned directory
fcs_files <- list.files(cleaned_dir, pattern = "\\.fcs$", full.names = TRUE)

# Read file batch ID
file_batch_id <- read.csv(file.path(meta_dir, "filebatch.csv"), row.names = 1)

# Define the debarcode_files function
debarcode_files <- function(fcs_files, 
                            file_batch_id,
                            debarcoding_key_dir,
                            out_dir = getwd(), 
                            min_threshold = TRUE, 
                            threshold = 0.18, 
                            to_plot = TRUE, 
                            barcodes_used = NULL, 
                            less_than_th = FALSE){
  
  # Check for duplicated file names
  if(anyDuplicated(fcs_files) != 0){
    stop("Names of FCS files are duplicated")
  }
  
  # Iterate over each FCS file
  for (file in fcs_files){
    print(paste0("   ", Sys.time()))
    print(paste0("   Debarcoding ", file))
    
    # Read the FCS file
    ff <- flowCore::read.FCS(file, transformation = FALSE)
    file_id <- which(file == fcs_files)
    batch_id <- file_batch_id[file_id, 1]
    
    # Get corresponding debarcoding key file
    key_file <- file.path(debarcoding_key_dir, gsub("_cleaned.fcs", ".csv", basename(file)))
    
    # Check if the key file exists
    if (!file.exists(key_file)) {
      warning(paste0("Debarcoding key file not found for ", basename(file)))
      next
    }
    
    # Read debarcoding key
debarcode_key <- read.csv(key_file, row.names = 1)

# Remove the "X" prefix if present
colnames(debarcode_key) <- sub("^X", "", colnames(debarcode_key))

# Check if column names are numeric
if (!all(grepl("^\\d+$", colnames(debarcode_key)))) {
  warning("Column names of debarcoding key should be numeric")
  print("Non-numeric column names:")
  print(colnames(debarcode_key)[!grepl("^\\d+$", colnames(debarcode_key))])
  next
}

    
    # Create output directory if it doesn't exist
    if(!dir.exists(out_dir)) dir.create(out_dir)
    
    # Select barcodes to use for debarcoding
    s_key <- debarcode_key
    
    # Prepare data for debarcoding
    dat <- CATALYST::prepData(ff)
    dat <- CATALYST::assignPrelim(dat, bc_key = s_key)
    
    # Estimate cutoffs for barcoding
    dat <- CATALYST::estCutoffs(dat)
    
    # Check and adjust cutoffs if necessary
    less_than_th <- c()
    if (min_threshold == TRUE){
      # Handle NA cutoffs by setting them to a default value (e.g., 1)
      metadata(dat)$sep_cutoffs[is.na(metadata(dat)$sep_cutoffs)] <- 1
      
      if(any(metadata(dat)$sep_cutoffs < threshold)){
        warning(paste0("Cutoff lower than ", threshold, " has been detected for ", basename(file), 
                      ", cutoff will be set to ", threshold))
        less_than_th <- c(less_than_th, basename(file))
      }
      
      id <- metadata(dat)$sep_cutoffs < threshold
      metadata(dat)$sep_cutoffs[id] <- threshold 
    } else {
      if(any(metadata(dat)$sep_cutoffs < threshold)){
        warning(paste0("Cutoff lower than ", threshold, " detected for ", basename(file))) 
        less_than_th <- c(less_than_th, basename(file))
      }
    }
    
    # Plot yields if requested
    if (to_plot == TRUE){
      p <- CATALYST::plotYields(dat, which = rownames(s_key))
      
      pdf(file.path(out_dir, paste(gsub(".fcs", "_yields.pdf", basename(file)))))
      for (name in names(p)){
        print(p[[name]])
      }
      dev.off()
    }
    
    # Apply the cutoffs to the data
    dat <- CATALYST::applyCutoffs(dat)
    
    # Plot debarcoding quality if requested
    if (to_plot == TRUE){
      p <- CATALYST::plotEvents(dat, n = 500)
      
      pdf(file.path(out_dir, paste(gsub(".fcs", "_debarcode_quality.pdf", 
                                        basename(file)))))
      for (name in names(p)){
        print(p[[name]])
      }
      dev.off()
    }
    
    # Remove events with bc_id = 0 (unassigned barcodes)
    dat <- dat[, dat$bc_id != 0]
    fs <- CATALYST::sce2fcs(dat, split_by = "bc_id")
  
    # Create a subdirectory for the batch
    tmp_dir <- file.path(out_dir, as.character(batch_id))
    if(!dir.exists(tmp_dir)) dir.create(tmp_dir)
    
    # Generate filenames for the debarcoded FCS files
    file_name <- gsub("_cleaned.fcs|.fcs", "", basename(file))
    
    flowCore::write.flowSet(fs, outdir = tmp_dir, 
                            filename = paste0(rownames(fs@phenoData), "_", file_name, 
                                              "_debarcoded.fcs")) 
  }
  
  # Save the names of files with lower debarcoding threshold if requested
  if(less_than_th == TRUE){
    saveRDS(less_than_th, file.path(out_dir, "files_with_lower_debarcoding_threshold.RDS"))
  }
}

# Apply debarcoding to all .fcs files
debarcode_files(fcs_files = fcs_files, 
                file_batch_id = file_batch_id, 
                debarcoding_key_dir = debarcoding_key_dir,
                out_dir = debarcoded_dir, 
                min_threshold = TRUE, 
                threshold = 0.18, 
                to_plot = TRUE, 
                barcodes_used = NULL, 
                less_than_th = FALSE)

# Load required packages
library(flowCore)
library(flowWorkspace)
library(stringr)
library(CytoNorm)
library(flowDensity)

# Define gating functions (assuming they are sourced from functions.R)
source("/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/SCRIPTS/functions.R")

# Define directories
debarcoded_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/DEBARCODED/"
gate_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/GATED/"
plot_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/GATED_PLOTS/"

# Create the gate and plot directories if they don't exist
if (!dir.exists(gate_dir)) dir.create(gate_dir)
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# List subdirectories containing modified .fcs files
subdirs <- list.dirs(debarcoded_dir, recursive = FALSE)
subdirs <- subdirs[grep("_Normalized\\.fcs$", subdirs)]

# Function to gate live cells and save density plot
gate_live_cells_and_plot <- function(flow_frame, file_name, viability_channel, plot_dir,
                                     tinypeak_removal_viability = 0.8, alpha_viability = 0.1,
                                     tinypeak_removal_Iridium = 0.8, alpha_Iridium = 0.05,
                                     arcsine_transform = TRUE, density = TRUE, man_th = 1.8, ...) {
  
  ff <- flow_frame
  
  # Perform arcsine transformation if needed
  if (arcsine_transform) {
    ff_t <- flowCore::transform(ff, transformList(colnames(ff)[grep("Di", colnames(ff))], CytoNorm::cytofTransform))
  } else {
    ff_t <- ff
  }
  
  # Perform density-based gating if specified
  if (density) {
    selection <- matrix(TRUE, nrow = nrow(ff), ncol = 1, dimnames = list(NULL, c("live")))
    
    tr <- list()
    for (m in c("Ir191Di", viability_channel)) {
      if (m == viability_channel) {
        upper <- TRUE
        alpha <- alpha_viability
        tr[[m]] <- flowDensity::deGate(ff_t, m, tinypeak.removal = tinypeak_removal_viability, 
                                       upper = upper, use.upper = TRUE, alpha = alpha, verbose = FALSE, count.lim = 3)
      } else {
        alpha <- alpha_Iridium
        tr[[m]] <- c(flowDensity::deGate(ff_t, m, tinypeak.removal = tinypeak_removal_Iridium, 
                                         upper = FALSE, use.upper = TRUE, alpha = alpha, verbose = FALSE, count.lim = 3), 
                     flowDensity::deGate(ff_t, m, tinypeak.removal = tinypeak_removal_Iridium, 
                                         upper = TRUE, use.upper = TRUE, alpha = alpha, verbose = FALSE, count.lim = 3)) 
      }
    }
    
    for (m in c(viability_channel, "Ir191Di")) {
      if (m == viability_channel) {
        selection[ff_t@exprs[, m] > tr[[m]][1], "live"] <- FALSE 
      } else {
        selection[ff_t@exprs[, m] < tr[[m]][1], "live"] <- FALSE
        selection[ff_t@exprs[, m] > tr[[m]][2], "live"] <- FALSE  
      }
    }
    
    percentage <- (sum(selection) / length(selection)) * 100
    flowDensity::plotDens(ff_t, c(viability_channel, "Ir191Di"), 
                          main = paste0(file_name, " (", format(round(percentage, 2), nsmall = 2), "%)"),
                          xlim = c(0, 8), ylim = c(0, 8), ...)
    
    abline(h = tr[["Ir191Di"]])
    abline(v = tr[[viability_channel]])
    
    points(ff_t@exprs[!selection[,"live"], c(viability_channel, "Ir191Di")], pch = ".") 
    
    ff <- ff[selection[,"live"], ]
    
  } else {
    # Perform gating by manual threshold
    th <- man_th
    cells_to_keep <- ff_t@exprs[, viability_channel] < th
    percentage <- (sum(cells_to_keep) / length(cells_to_keep)) * 100
    
    flowDensity::plotDens(ff_t, c(viability_channel, "Ir191Di"), 
                          main = paste0(file_name, " (", format(round(percentage, 2), nsmall = 2), "%)"),
                          xlim = c(0, 8), ylim = c(0, 8), ...)
    
    abline(v = th)
    
    points(ff_t@exprs[!cells_to_keep, c(viability_channel, "Ir191Di")], pch = ".") 
    
    ff <- ff[cells_to_keep, ]
  }
  
  # Save the density plot
  plot_file <- file.path(plot_dir, paste0(file_name, "_density_plot.png"))
  print(paste("Saving density plot:", plot_file))
  dev.copy(png, plot_file)
  dev.off()
  
  return(ff) # returns gated flow frame
}

# Process each subdirectory and gate files
for (subdir in subdirs) {
  files <- list.files(path = subdir, pattern = ".fcs$", full.names = TRUE)
  for (file in files) {
    print(paste("Processing file:", file))
    uID <- unlist(str_extract_all(basename(file), "(?<=_).*?(?=_)"))[1]
    
    # Read the .fcs file
    ff <- read.FCS(filename = file, transformation = FALSE)
    
    # Check if the object is of class flowFrame
    if (class(ff) != "flowFrame") {
      warning("Skipping file ", file, " as it is not of class 'flowFrame'.")
      next
    }
    
    # Print the original $TOT value
    print(paste("Original $TOT:", ff@description['$TOT']))
    print(paste("Number of cells before gating:", nrow(ff@exprs)))
    
    # Gate intact cells
    print("Gating intact cells...")
    ff <- gate_intact_cells(flow_frame = ff, file_name = basename(file), arcsine_transform = TRUE)
    print(paste("Number of cells after gating intact cells:", nrow(ff@exprs)))
    
    # Gate singlet cells
    print("Gating singlet cells...")
    ff <- gate_singlet_cells(flow_frame = ff, channels = "Event_length", file_name = basename(file), arcsine_transform = TRUE)
    print(paste("Number of cells after gating singlet cells:", nrow(ff@exprs)))
    
    # Gate immune cells
    print("Gating immune cells...")
    ff <- gate_out_non_immune(flow_frame = ff, bead_channel = "Y89", arcsine_transform = TRUE, density = FALSE, man_th = 1)
    print(paste("Number of cells after gating immune cells:", nrow(ff@exprs)))
    
    # Gate live cells and save density plot
    print("Gating live cells and saving density plot...")
    ff <- gate_live_cells_and_plot(flow_frame = ff, file_name = basename(file), viability_channel = "Y89Di",
                                   plot_dir = plot_dir, arcsine_transform = TRUE, density = TRUE, man_th = 1)
    print(paste("Number of cells after gating live cells:", nrow(ff@exprs)))
    
    # Get the number of cells after all gating steps
    n_cells <- nrow(ff@exprs)
    print(paste("Number of cells after all gating steps:", n_cells))
    
    # Update the $TOT keyword in the FCS header
    ff@description['$TOT'] <- as.character(n_cells)
    
    # Print the updated $TOT value
    print(paste("Updated $TOT:", ff@description['$TOT']))
    
    # Write gated .fcs file to the gate directory
    write.FCS(ff, file.path(gate_dir, gsub(".fcs", "_gated.fcs", basename(file))))
  }
}

# Load necessary library
library(flowCore)

# Define directories
original_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/DEBARCODED/"
gated_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/GATED/"

# List all debarcoded FCS files in subfolders
list_debarcoded_files <- function(dir) {
  # List all subdirectories
  subdirs <- list.dirs(dir, full.names = TRUE, recursive = TRUE)
  # List all FCS files in subdirectories
  fcs_files <- unlist(lapply(subdirs, function(x) list.files(x, pattern = "_debarcoded\\.fcs$", full.names = TRUE)))
  return(fcs_files)
}

# Get list of original (debarcoded) and gated FCS files
original_files <- list_debarcoded_files(original_dir)
gated_files <- list.files(gated_dir, pattern = "_gated\\.fcs$", full.names = TRUE)

# Extract base names for matching
original_basenames <- gsub("_debarcoded\\.fcs$", "", basename(original_files))
gated_basenames <- gsub("_debarcoded_gated\\.fcs$", "", basename(gated_files))

# Initialize vectors to store cell counts
original_counts <- numeric(length(original_files))
gated_counts <- numeric(length(gated_files))

# Calculate cell counts for original files
for (i in seq_along(original_files)) {
  ff_original <- read.FCS(original_files[i], transformation = FALSE)
  original_counts[i] <- nrow(ff_original@exprs)
}

# Calculate cell counts for gated files
for (i in seq_along(gated_files)) {
  ff_gated <- read.FCS(gated_files[i], transformation = FALSE)
  gated_counts[i] <- nrow(ff_gated@exprs)
}

# Create a named vector for easier matching
names(original_counts) <- gsub("_debarcoded\\.fcs$", "", basename(original_files))
names(gated_counts) <- gsub("_debarcoded_gated\\.fcs$", "", basename(gated_files))

# Match the counts
matched_original_counts <- original_counts[names(original_counts) %in% names(gated_counts)]
matched_gated_counts <- gated_counts[names(gated_counts) %in% names(matched_original_counts)]

# Calculate the difference in cell counts
cell_count_difference <- matched_original_counts - matched_gated_counts

# Print the results
cat("Original cell counts:\n")
print(matched_original_counts)
cat("Gated cell counts:\n")
print(matched_gated_counts)
cat("Difference in cell counts:\n")
print(cell_count_difference)

# Load necessary libraries
library(flowCore)
library(ggplot2)
library(dplyr)

# Define the path to the directory containing the FCS files
base_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/DEBARCODED/"

# List all FCS files in the directory and its subdirectories
fcs_files <- list.files(base_dir, pattern = "debarcoded.fcs$", recursive = TRUE, full.names = TRUE)

# Initialize a data frame to store the results
cell_counts <- data.frame(Sample_ID = character(), Input_cells = numeric(), stringsAsFactors = FALSE)

# Loop through each FCS file and read the cell count
for (file in fcs_files) {
  # Read the FCS file
  fcs <- read.FCS(file)
  
  # Get the number of cells (events) in the file
  n_cells <- nrow(exprs(fcs))
  
  # Extract the Sample ID from the file name
  sample_id <- sub("_.*", "", basename(file))
  
  # Add the results to the data frame
  cell_counts <- rbind(cell_counts, data.frame(Sample_ID = sample_id, Input_cells = n_cells))
}

# Calculate the median cell count
median_cells <- median(cell_counts$Input_cells)

# Identify the number of cells above the median
above_median_count <- sum(cell_counts$Input_cells > median_cells)

# Calculate the threshold for the majority of samples
majority_threshold <- quantile(cell_counts$Input_cells, probs = 0.5)  # Adjust the quantile as needed, here it represents the median

# Identify the number of samples with cell counts above the threshold
above_threshold_count <- sum(cell_counts$Input_cells > majority_threshold)

# Identify the minimum cell count
min_cells <- min(cell_counts$Input_cells)

# Plot the cell count distribution
p <- ggplot(cell_counts, aes(x = Sample_ID, y = Input_cells)) +
  geom_bar(stat = "identity", fill = "darkgray") +
  geom_hline(aes(yintercept = median_cells), color = "red", linetype = "dashed", size = 1) +
  geom_hline(aes(yintercept = min_cells), color = "pink", linetype = "dashed", size = 1) +
  annotate("text", x = 0.5, y = min_cells, label = min_cells, color = "pink", size = 4, hjust = 0, vjust = -1) +
  annotate("text", x = 0.5, y = median_cells, label = majority_threshold, color = "red", size = 4, hjust = 0, vjust = -1) +
  labs(title = "Cell count distribution", x = "Sample ID", y = "Input cells") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Print the plot
print(p)

# Save the plot
ggsave("cell_count_distribution.png", p, width = 15, height = 7)

library(ggplot2)

# Convert the signal distribution matrix to a data frame
signal_df <- as.data.frame(signal_distribution)

# Add an index column to the data frame
signal_df$Index <- seq_len(nrow(signal_df))

# Reshape the data frame from wide to long format for easier plotting
library(tidyr)
signal_long <- pivot_longer(signal_df, cols = -Index, names_to = "Channel", values_to = "Intensity")

# Extract marker names from the channel names
signal_long$Marker <- gsub("_.*", "", signal_long$Channel)

# Log-transform the intensity values
signal_long$log_intensity <- log(signal_long$Intensity + 0.1)

# Plot colorful density plots for each marker
marker_plots <- ggplot(signal_long, aes(x = log_intensity, fill = Channel)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Marker, scales = "free") +
  scale_fill_viridis_d(option = "C") +  # Choose a rich color palette
  theme_minimal() +
  labs(title = "Signal Intensity Distribution for Each Marker", x = "Log(Intensity + 0.1)", y = "Marker") +
  theme(legend.position = "none")  # Remove legend to declutter the plot

# Print the plots
print(marker_plots)

# Save the plot as a PNG file
ggsave("marker_density_plots.png", marker_plots, width = 10, height = 6)

```

```{r}
# Define source and destination directories
source_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/DEBARCODED/PLOTS"
dest_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/"

# List all files in the source directory
files_to_move <- list.files(source_dir, full.names = TRUE)

# Move each file to the destination directory
for (file_path in files_to_move) {
  # Extract the file name
  file_name <- basename(file_path)
  
  # Define the new file path in the destination directory
  new_file_path <- file.path(dest_dir, file_name)
  
  # Move the file
  file.rename(file_path, new_file_path)
}

# Optionally, remove the now empty PLOTS directory
unlink(source_dir, recursive = TRUE)

# Load necessary libraries
library(SingleCellExperiment)
library(scater)
library(ggplot2)

# Define the file path for the input subsampled SCE object
input_file_path <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/subsampled_sce_1000_cells_per_sample.rds"

# Read the subsampled SCE object
sce <- readRDS(input_file_path)

# List of immune markers of interest
markers_of_interest <- c(
  "CD45", "CCR6", "CD123", "CD19", "CD4", "CD8a", "CD11c", "CD16",
  "CD45RO", "CD45RA", "CD161", "CCR4", "CD25", "CD27", "CD57",
  "CXCR3", "CXCR5", "CD28", "CD38", "CD56", "TCRgd", "CRTH2",
  "CCR7", "CD14", "CD3", "CD20", "CD66b", "HLA-DR", "IgD", "CD127"
)

# Subset the counts data for these markers
counts_subset <- counts(sce)[markers_of_interest, ]

# Set seed for reproducibility
set.seed(123)

# Perform PCA using the counts assay
sce <- runPCA(sce, exprs_values = "counts")

# Create a PCA plot focusing on Condition
pca_data <- reducedDim(sce, "PCA")
pca_df <- as.data.frame(pca_data)
pca_df$Condition <- colData(sce)$Condition

# Plot the PCA results with adjustments
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 1, alpha = 0.7) +  # Adjust point size and transparency
  labs(title = "PCA Plot of Immune Markers by Condition",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal() +
  theme(legend.position = "right", 
        legend.text = element_text(size = 10), 
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 12))

# Display the plot
print(p)

# Save the plot with increased size
output_plot_path <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/pca_condition_plot.pdf"
ggsave(output_plot_path, plot = p, width = 15, height = 12)  # Increased dimensions

message("PCA plot by Condition has been saved to ", output_plot_path)

# Load necessary libraries
library(SingleCellExperiment)
library(pheatmap)
library(dplyr)

# Define the file path for the input subsampled SCE object
input_file_path <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/subsampled_sce_1000_cells_per_sample.rds"

# Read the subsampled SCE object
sce <- readRDS(input_file_path)

# Access marker names from rowData
marker_names <- rownames(rowData(sce))

# Ensure marker_names is not empty
if (length(marker_names) == 0) {
  stop("No marker names found in rowData.")
}

# Define the markers of interest
markers_of_interest <- c("CD45", "CCR6", "CD123", "CD19", "CD4", "CD8a", "CD11c", "CD16", 
                         "CD45RO", "CD45RA", "CD161", "CCR4", "CD25", "CD27", "CD57", 
                         "CXCR3", "CXCR5", "CD28", "CD38", "CD56", "TCRgd", "CRTH2", 
                         "CCR7", "CD14", "CD3", "CD20", "CD66b", "HLA-DR", "IgD", "CD127")

# Check if all markers of interest are present in marker_names
missing_markers <- setdiff(markers_of_interest, marker_names)
if (length(missing_markers) > 0) {
  warning("The following markers are not found in the dataset: ", paste(missing_markers, collapse = ", "))
  # Proceed with the available markers
  markers_of_interest <- intersect(markers_of_interest, marker_names)
}

# Plot marker ranking based on non-redundancy scores
# Ensure you have the necessary function plotNRS loaded or defined
plotNRS(sce, features = markers_of_interest, color_by = "Condition")

# Adjust margins if necessary
par(mar = c(5, 8, 4, 2) + 0.1)

# Save the plot
png("marker_ranking_plot.png", width = 10, height = 8, units = "in", res = 300)  # Adjust width, height, and resolution as needed
plotNRS(sce, features = markers_of_interest, color_by = "Condition")
dev.off()  # Close the PNG device

# Load necessary libraries
library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(umap)

# Define the file path for the input subsampled SCE object
input_file_path <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/subsampled_sce_1000_cells_per_sample.rds"

# Read the subsampled SCE object
sce <- readRDS(input_file_path)

# Set seed for reproducibility
set.seed(123)

# Prepare data for UMAP (using counts directly)
counts_data <- counts(sce)

# Run UMAP
umap_result <- umap(t(counts_data))  # Transpose if necessary, as UMAP expects features in rows

# Create a UMAP data frame
umap_df <- as.data.frame(umap_result$layout)
umap_df$Condition <- colData(sce)$Condition

# Define color palette
condition_colors <- scale_color_manual(values = RColorBrewer::brewer.pal(n = length(unique(umap_df$Condition)), name = "Set1"))

# UMAP plot
p_umap <- ggplot(umap_df, aes(x = V1, y = V2, color = Condition)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "UMAP Plot of Immune Markers by Condition",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2") +
  condition_colors +
  theme_minimal() +
  theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 10), 
        plot.title = element_text(size = 16),
        axis.title = element_text(size = 12))

# Display the UMAP plot
print(p_umap)

# Save the plot
output_plot_path <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/umap_condition_plot.pdf"
ggsave(output_plot_path, plot = p_umap, width = 15, height = 12)

message("UMAP plot has been saved.")

```

```{r}
# Define color palette
condition_colors <- RColorBrewer::brewer.pal(n = length(unique(umap_df$Condition)), name = "Set1")

# Generate UMAP plots for each condition
conditions <- unique(umap_df$Condition)

for (condition in conditions) {
    # Subset data for the current condition
    subset_data <- umap_df[umap_df$Condition == condition, ]
    
    # Get color for the current condition
    condition_color <- condition_colors[which(unique(umap_df$Condition) == condition)]
    
    # Create UMAP plot for the current condition
    p_umap_condition <- ggplot(subset_data, aes(x = V1, y = V2)) +
        geom_point(size = 2, alpha = 0.7, color = condition_color) +  # Use the specific color for the condition
        labs(title = paste("UMAP Plot for Condition:", condition),
             x = "UMAP Dimension 1",
             y = "UMAP Dimension 2") +
        theme_minimal() +
        theme(plot.title = element_text(size = 16),
              axis.title = element_text(size = 12))

    # Display the UMAP plot
    print(p_umap_condition)

    # Save the plot for the current condition
    output_plot_path <- paste0("/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/umap_", gsub(" ", "_", condition), "_plot.pdf")
    ggsave(output_plot_path, plot = p_umap_condition, width = 15, height = 12)
    
    message(paste("UMAP plot for condition", condition, "has been saved to", output_plot_path))
}

# Create flowFrame object
flow_frame <- flowFrame(expr_data_transformed)

# Read input
flowSOM_res <- ReadInput(flow_frame)
flowSOM_res <- BuildSOM(flowSOM_res, colsToUse = markers)
flowSOM_res <- BuildMST(flowSOM_res)

# Meta-clustering using metaClustering_consensus
k <- 35
metaClustering_res <- metaClustering_consensus(flowSOM_res$map$codes, k = k)
flowSOM_res$metaClustering <- metaClustering_res

# Basic ggplot2 plot without advanced customization
png(filename = "FlowSOM_Clusters.png", width = 4000, height = 4000, res = 300)
ggplot(plot_data, aes(x = X, y = Y, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "FlowSOM Clusters", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal()
dev.off()

# PlotStars function with correct background values
backgroundValues <- flowSOM_res$metaClustering

# Basic ggplot2 plot without advanced customization
png(filename = "FlowSOM_Clusters.png", width = 4000, height = 4000, res = 300)
ggplot(plot_data, aes(x = X, y = Y, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "FlowSOM Clusters", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal()
dev.off()

# Load necessary libraries
library(SingleCellExperiment)
library(ggplot2)
library(umap)
library(dplyr)
library(flowCore)

# Load your data
sce <- readRDS("/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/updated_sce_with_modified_sample_ids.rds")

# Load clustering results
clustering_results <- readRDS("/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/flowSOM_clustering_results.rds")

# Extract expression data and metadata
expr_data <- assay(sce)
metadata <- colData(sce)

# Select markers
desired_markers <- c("CD45", "CD123", "CD19", "CD4", "CD8a", "CD11c", "CD16", "CD45RO", "CD45RA",
                     "CD161", "CCR4", "CD25", "CD27", "CD57", "CXCR5", "CD28", "CD38", "CD56",
                     "TCRgd", "CRTH2", "CCR7", "CD14", "CD3", "CD20", "CD66b", "HLA-DR", "IgD", "CD127")

# Filter markers to include only those present in the data
available_markers <- intersect(desired_markers, rownames(expr_data))

# Transpose and subset expression data
expr_data_subset <- expr_data[available_markers, ]

# Subset 1000 cells per sample
set.seed(123)
sample_ids <- unique(metadata$sample_id)
subsetted_indices <- unlist(lapply(sample_ids, function(sample_id) {
  sample_cells <- which(metadata$sample_id == sample_id)
  if(length(sample_cells) > 1000) {
    return(sample(sample_cells, 1000))
  } else {
    return(sample_cells)
  }
}))
expr_data_subset <- expr_data_subset[, subsetted_indices]
metadata_subset <- metadata[subsetted_indices, ]

# Transpose to ensure samples are rows and markers are columns
expr_data_subset <- t(expr_data_subset)

# Ensure it's a numeric matrix and set column names
expr_data_subset <- as.matrix(expr_data_subset)
colnames(expr_data_subset) <- available_markers

# Transform data
expr_data_transformed_subset <- asinh(expr_data_subset / 5)

# Extract clustering labels for the subsetted cells
cell_clusters_subset <- clustering_results$metaClustering[subsetted_indices]

# Perform UMAP
set.seed(123)
umap_res <- umap::umap(expr_data_transformed_subset, n_neighbors = 30, min_dist = 0.3)
umap_data <- data.frame(umap_res$layout)
colnames(umap_data) <- c("UMAP1", "UMAP2")
umap_data$Sample <- metadata_subset$sample_id

# Determine global min and max for color scaling
global_min <- 0
global_max <- 5.644849

# Create directory for saving plots
output_dir <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/UMAP_marker_expression"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Combine UMAP coordinates with marker expressions and plot
for (marker in available_markers) {
  umap_data[[marker]] <- expr_data_transformed_subset[, marker]
  
  # Plot UMAP with marker expression
  umap_plot <- ggplot(umap_data, aes_string(x = "UMAP1", y = "UMAP2", color = marker)) +
    geom_point(size = 1, alpha = 0.6) +
    scale_color_viridis_c(limits = c(global_min, global_max), na.value = "grey") +  # Apply consistent color scale
    labs(title = paste("UMAP Plot of", marker, "Expression"), x = "UMAP1", y = "UMAP2", color = marker) +
    theme_minimal() +
    guides(color = guide_colorbar(title = marker))
  
  # Save UMAP plot for the marker
  file_name <- file.path(output_dir, paste0("UMAP_Plot_", marker, ".png"))
  png(filename = file_name, width = 4000, height = 4000, res = 300)
  print(umap_plot)
  dev.off()
}


# Load necessary libraries
library(ComplexHeatmap)
library(circlize) # for color functions

# Define color scale
color_scale <- colorRamp2(c(min(heatmap_data), median(heatmap_data), max(heatmap_data)),
                          c("blue", "white", "red"))

# Compute a statistic for each column for the top annotation
# Example: Sum of values for each marker
col_stats <- colSums(heatmap_data)

# Create top annotation
top_annotation <- HeatmapAnnotation(
  barplot = anno_barplot(
    col_stats,
    height = unit(2, "cm"),
    width = unit(1, "cm"),
    ylim = c(0, max(col_stats))
  )
)

# Create Heatmap with the correct top annotation
heatmap <- Heatmap(heatmap_data,
        name = "Median Marker Intensity",
        col = color_scale,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        row_title = "Metaclusters",
        column_title = "Markers",
        row_names_side = "left",
        column_names_side = "top",
        top_annotation = top_annotation,
        show_row_names = TRUE,
        show_column_names = TRUE
)

# Additional customization for row sizes
row_annot <- rowAnnotation(
  size = anno_barplot(rowSums(heatmap_data)),
  show_annotation_name = FALSE,
  annotation_height = unit(4, "cm")
)

# Add relative sizes of clusters on the side of the heatmap
heatmap_with_size <- Heatmap(heatmap_data,
        name = "Median Marker Intensity",
        col = color_scale,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        row_title = "Metaclusters",
        column_title = "Markers",
        row_names_side = "left",
        column_names_side = "top",
        right_annotation = row_annot
)

# Draw heatmap
draw(heatmap_with_size)

# Load necessary libraries
library(SingleCellExperiment)
library(ggplot2)
library(tidyr)
library(dplyr)

# Load the SingleCellExperiment object
sce <- readRDS("/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/subsampled_sce_1000_cells_per_sample.rds")

# Check the row names of the counts assay to verify markers
marker_names <- rownames(assay(sce, "counts"))
print(marker_names)  # Check if the markers are here

# Define the markers based on actual row names (if they match)
markers <- intersect(c("CD45", "CCR6", "CD123", "CD19", "CD4", "CD8a", "CD11c", 
                       "CD16", "CD45RO", "CD45RA", "CD161", "CCR4", "CD25", 
                       "CD27", "CD57", "CXCR3", "CXCR5", "CD28", "CD38", 
                       "CD56", "TCRgd", "CRTH2", "CCR7", "CD14", "CD3", 
                       "CD20", "CD66b", "HLA-DR", "IgD", "CD127"), marker_names)

# Combine assay data with colData
counts_data <- as.data.frame(t(assay(sce, "counts")))  # Transpose to have cells as rows
combined_data <- cbind(as.data.frame(colData(sce)), counts_data)  # Convert colData to a data frame

# Convert to long format for plotting
freq_data_long <- combined_data %>%
  pivot_longer(cols = all_of(markers), names_to = "thr_major_type", values_to = "frequency")

# Frequency bar plot
freq_plot <- ggplot(freq_data_long, aes(x = sample_id, fill = thr_major_type)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Frequency")

# Display the frequency plot
print(freq_plot)

# Stacked bar plot by Condition
condition_plot <- ggplot(combined_data, aes(x = Condition, fill = thr_major_type)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(y = "Frequency", x = "Condition")

# Display the condition plot
print(condition_plot)

# Violin plots for various cell types across conditions
violin_plots <- combined_data %>%
  pivot_longer(cols = all_of(markers), names_to = "thr_major_type", values_to = "count") %>%
  ggplot(aes(x = Condition, y = count, fill = Condition)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ thr_major_type, scales = "free_y") +
  theme_minimal() +
  labs(y = "Count", x = "Condition")

# Display the violin plots
print(violin_plots)

 Load required libraries
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)

# Path to the RDS file
file_path <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/updated_sce_with_metacluster_annotations.rds"

# Load the SingleCellExperiment object
sce_updated <- readRDS(file_path)

# Verify the object
print(sce_updated)

# Convert colData to a data frame for easier manipulation with dplyr
col_data_df <- as.data.frame(colData(sce_updated))

# Check for the presence of the 'Condition' column
if (!"Condition" %in% colnames(col_data_df)) {
  stop("The 'Condition' column is missing from the data.")
}

# Calculate the proportions
proportions_df <- col_data_df %>%
  group_by(Sample_ID, Major_lineage, Condition) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(Sample_ID) %>%
  mutate(total_cells = sum(cell_count),
         proportion = cell_count / total_cells) %>%
  ungroup()

# Print the proportions data frame
print(proportions_df)

# Perform Kruskal-Wallis test for each Major_lineage
test_results <- proportions_df %>%
  group_by(Major_lineage) %>%
  do({
    df_subset <- .
    if (nrow(df_subset) > 1) {  # Check if there is more than one row to test
      test_result <- kruskal.test(proportion ~ Condition, data = df_subset)
      tibble(p_value = test_result$p.value)
    } else {
      tibble(p_value = NA_real_)
    }
  }) %>%
  ungroup()

# Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
test_results <- test_results %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH"))

# Print the test results
print(test_results)

# Create boxplots with consistent y-axis scaling across all facets
plot <- ggplot(proportions_df, aes(x = Condition, y = proportion, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~ Major_lineage, scales = "fixed") +  # Use fixed scales for uniform y-axis
  theme_minimal(base_size = 18) +  # Larger base size for better legibility
  labs(title = "Proportion of Major Lineages Across Sample Groups",
       x = "Condition",
       y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.position = "bottom",
        plot.title = element_text(size = 18, face = "bold"))

# Save the plot with higher resolution and larger size
output_plot_path = "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/major_lineage_proportions_plot.png"
ggsave(output_plot_path, plot, width = 20, height = 15, dpi = 300)  # Increased dimensions for better legibility

# Save the test results to a CSV file
output_test_results <- "/rds/projects/m/mossp-pancreatic-cancer/TASSELLA/major_lineage_test_results.csv"
write.csv(test_results, output_test_results, row.names = FALSE)

# Confirm the file is saved
cat("Test results have been saved to", output_test_results, "\n")

# Load necessary libraries
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)

# Assuming sce_updated is your SingleCellExperiment object and it includes Sample_ID, Major_lineage, Condition

# Extract expression data and combine with cell metadata
exprs_data <- as.data.frame(t(assay(sce_updated)))
col_data <- colData(sce_updated)

# Ensure the column names are correctly assigned
exprs_data$Sample_ID <- col_data$Sample_ID
exprs_data$Major_lineage <- col_data$Major_lineage
exprs_data$Condition <- col_data$Condition

# Check the columns to ensure they exist
print(head(exprs_data))

# Calculate median intensity for each marker within each Major_lineage and Sample_ID
median_intensity_df <- exprs_data %>%
  pivot_longer(cols = -c(Sample_ID, Major_lineage, Condition), names_to = "marker", values_to = "expression") %>%
  group_by(Sample_ID, Major_lineage, Condition, marker) %>%
  summarise(median_intensity = median(expression, na.rm = TRUE)) %>%
  ungroup()

# Perform Kruskal-Wallis test for each marker within each Major_lineage
kw_test_results <- median_intensity_df %>%
  group_by(Major_lineage, marker) %>%
  summarise(p_value = kruskal.test(median_intensity ~ Condition, data = .)$p.value,
            .groups = 'drop')

# Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
kw_test_results <- kw_test_results %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH"))

# Print the test results
print(kw_test_results)

# Create boxplots for marker median intensities
ggplot(median_intensity_df, aes(x = Condition, y = median_intensity, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~ Major_lineage + marker, scales = "free_y") +
  theme_minimal() +
  labs(title = "Marker Median Intensities Across Sample Groups",
       x = "Condition",
       y = "Median Intensity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

```
```{r}
# Load necessary libraries
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)

# Assuming sce_updated is your SingleCellExperiment object and it includes Sample_ID, Major_lineage, Condition

# Extract expression data and combine with cell metadata
exprs_data <- as.data.frame(t(assay(sce_updated)))
col_data <- colData(sce_updated)

# Ensure the column names are correctly assigned
exprs_data$Sample_ID <- col_data$Sample_ID
exprs_data$Major_lineage <- col_data$Major_lineage
exprs_data$Condition <- col_data$Condition

# Check the columns to ensure they exist
print(head(exprs_data))

# Calculate median intensity for each marker within each Major_lineage and Sample_ID
median_intensity_df <- exprs_data %>%
  pivot_longer(cols = -c(Sample_ID, Major_lineage, Condition), names_to = "marker", values_to = "expression") %>%
  group_by(Sample_ID, Major_lineage, Condition, marker) %>%
  summarise(median_intensity = median(expression, na.rm = TRUE), .groups = 'drop')

# Perform Kruskal-Wallis test for each marker within each Major_lineage
kw_test_results <- median_intensity_df %>%
  group_by(Major_lineage, marker) %>%
  summarise(p_value = kruskal.test(median_intensity ~ Condition, data = .)$p.value,
            .groups = 'drop')

# Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
kw_test_results <- kw_test_results %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH"))

# Print the test results
print(kw_test_results)

# Define the number of pages based on the number of unique markers and Major_lineage combinations
n_markers <- length(unique(median_intensity_df$marker))
n_major_lineages <- length(unique(median_intensity_df$Major_lineage))
plots_per_page <- 15  # Number of facets per page

# Create boxplots for marker median intensities with pagination
n_pages <- ceiling(n_markers * n_major_lineages / plots_per_page)

# Save each page of the plot
for (i in 1:n_pages) {
  plot <- ggplot(median_intensity_df, aes(x = Condition, y = median_intensity, fill = Condition)) +
    geom_boxplot() +
    facet_wrap_paginate(~ Major_lineage + marker, scales = "free_y", ncol = 5, nrow = 3, page = i) +
    theme_minimal() +
    labs(title = paste("Marker Median Intensities Across Sample Groups - Page", i),
         x = "Condition",
         y = "Median Intensity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the plot to a file
  ggsave(paste0("marker_median_intensities_boxplot_page_", i, ".png"), plot, width = 20, height = 15, units = "in", dpi = 300)
}

```
```{r}
# Create a large boxplot for marker median intensities
plot <- ggplot(median_intensity_df, aes(x = Condition, y = median_intensity, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~ Major_lineage + marker, scales = "free_y") +
  theme_minimal(base_size = 16) +
  labs(title = "Marker Median Intensities Across Sample Groups",
       x = "Condition",
       y = "Median Intensity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5))

# Save the plot with large dimensions
ggsave("marker_median_intensities_boxplot_large.png", plot, width = 40, height = 30, units = "in", dpi = 300)
```

