# Load necessary library for data manipulation and visualization
library(tidyverse)

# Load the provided dataset which contains 'Y' (gene expression matrix) and 'patient_data'
load("assess_data_resit-23.Rdata")

# ------------- Part 1: PCA (Principal Component Analysis) -------------
# Apply log transformation to the gene expression matrix 'Y', adding 1 to avoid log(0) issues
# Then transpose the matrix because PCA in R works on rows as observations, but here the samples are columns
pca_x  <- t(log(Y + 1))

# Perform PCA (Principal Component Analysis) on the log-transformed and transposed matrix
# Centering and scaling the data
pca_res1 <- prcomp(pca_x, scale = TRUE, center = TRUE)

# Create a data frame that contains the PCA results ('pca_res1$x' holds the principal components)
# Add metadata for 'tissue' (tumour/normal) and 'patient' from 'patient_data'
pca_df1 <- data.frame(
  pca_res1$x,  # Principal component scores for each sample
  tissue = patient_data$tissue,  # Tissue type: tumour or normal
  patient = patient_data$patient  # Patient identifier
)

# Plot the PCA results using ggplot
# X-axis represents the first principal component (PC1), and Y-axis represents the second principal component (PC2)
# Color the points by tissue type (tumour or normal), so we can visually identify any differences or clusters
ggplot(pca_df1, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point(size = 3) +  # Plot the PCA points with a size of 3
  theme_bw() +  # Use a clean black and white theme for the plot
  labs(x = "PC1", y = "PC2") +  # Label the axes with PC1 and PC2
  theme(legend.position = "bottom")  # Move the legend to the bottom of the plot

# The PCA plot helps visually inspect whether tumour and normal tissues separate in the first two PCs,
# indicating variance in gene expression between them or identifying outliers that may indicate problems in the data.

# ------------- Part 2: Differential Expression Analysis -------------
# Load the MASS library which contains functions for statistical modeling, such as GLM (Generalized Linear Models)
library(MASS)

# Select an arbitrary gene index, here we are analyzing the gene at row 20 in the expression matrix Y
idx <- 20

# Define a subset of 20 samples (first 20 columns)
c_cl <- 1:20

# Extract the tissue type (tumour/normal) for the first 20 samples
x <- patient_data$tissue[c_cl]

# Extract the patient information for the first 20 samples
z <- patient_data$patient[c_cl]

# Create a temporary data frame for regression analysis, including:
# - 'y': Gene expression values for the selected gene (idx) across the first 20 samples
# - 'x': Tissue type (tumour/normal)
# - 'z': Patient identifiers to account for patient-specific variability
# - 'lib_size': Library size, i.e., the total expression (sum of all gene counts) for each sample
tmp <- data.frame(y = Y[idx, c_cl], x = x, z = z, lib_size = colSums(Y[, c_cl]))

# Fit a generalized linear model (GLM) with a Poisson distribution (suitable for count data)
# We are modeling the gene expression 'y' as a function of tissue type ('x'), patient ('z'), and library size ('lib_size')
# This helps identify if tissue type (tumour vs. normal) has a significant effect on gene expression
out <- glm(y ~ x + z + lib_size, data = tmp, family = "poisson")

# Extract the p-value for the tissue type variable ('x')
# This tells us if there is a significant difference in expression of the selected gene between tumour and normal tissue

# Load necessary libraries
library(tidyverse)

# Load the dataset
load("assess_data_resit-23.Rdata")

# Part 1: Perform PCA on the full dataset
# Transpose the gene expression matrix and apply log transformation
pca_x  <- t(log(Y + 1))

# Perform PCA with centering and scaling the data
pca_res <- prcomp(pca_x, scale = TRUE, center = TRUE)

# Create a dataframe to store PCA results with tissue and patient metadata
pca_df <- data.frame(
  pca_res$x,                      # Principal component scores
  tissue = patient_data$tissue,    # Tissue type (tumour/normal)
  patient = patient_data$patient   # Patient ID
)

# Part 2: Plot the first two principal components
ggplot(pca_df, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = "PC1", y = "PC2", title = "PCA of Gene Expression Data") +
  theme(legend.position = "bottom")

# Part 3: Investigate problematic samples visually
# Look for samples that appear far from the main cluster in the plot
# They may represent technical issues (outliers) or noise in the data.

# Identify problematic samples manually based on the scatter plot
# suppose patient 'P5' appears to be an outlier.
# Let's filter out patient 'P5' from further analysis

# Identify the problematic sample(s) - Replace 'P5' with actual problematic samples from the plot
problematic_samples <- c("P5")

# Part 4: Remove problematic samples
cleaned_pca_df <- pca_df[!pca_df$patient %in% problematic_samples, ]
cleaned_Y <- Y[, !patient_data$patient %in% problematic_samples]

# Part 5: Re-plot PCA without the problematic samples
ggplot(cleaned_pca_df, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = "PC1", y = "PC2", title = "PCA After Removing Problematic Samples") +
  theme(legend.position = "bottom")

p_val <- summary(out)$coefficients[2, 4]

# 'p_val' provides the p-value for the tissue type coefficient, which helps determine whether this gene is significantly
# differentially expressed between tumour and normal tissue after accounting for patient-specific variation and library size.

PCA: If tumour and normal tissues form separate clusters, it suggests a clear distinction in gene expression profiles between the two tissue types. Any outliers or overlapping points may indicate batch effects or noise in the data.
Differential Expression: The p-value from the GLM model gives statistical evidence on whether the selected gene is differentially expressed between tumour and normal tissue. A low p-value (typically < 0.05) indicates that the gene’s expression is significantly different between the two conditions.

Why are they problematic?
Problematic samples in PCA are usually those that do not cluster well with the rest of the data. These samples may:
Be technical outliers (due to issues during data collection or sequencing).
Have poor data quality (low library size or insufficient read depth).
Have unexpected biology (e.g., contamination or mislabeling of tumour/normal samples).
Removing these samples can help avoid misleading conclusions in downstream analyses.

After removing the problematic samples, you should re-run your analyses (such as PCA and differential expression) to ensure that your data is clean and that the results are reliable.

The glm() function in the MASS package will be used for Poisson regression. The key steps are:

Use the gene expression data (Y) and metadata (patient_data) to perform a Poisson regression.
For each gene, fit a Poisson model with tissue type (tumour or normal) as the main predictor.
Extract the p-values from the regression model and plot the log10(p-value) for each gene.

# Load necessary libraries
library(tidyverse)
library(MASS)


# Initialize a vector to store the p-values for each gene
p_values <- numeric(nrow(cleaned_Y))

# Loop over each gene and fit a Poisson regression model
for (i in 1:nrow(cleaned_Y)) {
  # Subset the gene expression data for the current gene
  y <- cleaned_Y[i, ]
  
  # Create a temporary dataframe with tissue type and patient information
  tmp <- data.frame(
    y = y, 
    tissue = cleaned_patient_data$tissue,  # Tumour or Normal
    patient = cleaned_patient_data$patient,  # Patient ID
    lib_size = colSums(cleaned_Y)  # Library size (sum of counts per sample)
  )
  
  # Fit a Poisson regression model
  model <- glm(y ~ tissue + lib_size, data = tmp, family = "poisson")
  
  # Extract the p-value for the tissue variable (indicating tumour vs. normal)
  p_values[i] <- summary(model)$coefficients[2, 4]
}

# Calculate the -log10(p-value) for each gene
log10_p_values <- -log10(p_values)

# Plot the log10 p-values
plot(log10_p_values, pch = 20, col = "blue",
     main = "Log10(p-values) for Differential Expression",
     xlab = "Gene Index", ylab = "-log10(p-value)")

# Highlight significant genes (p < 0.05)
significant_genes <- which(p_values < 0.05)
points(significant_genes, log10_p_values[significant_genes], col = "red", pch = 20)

Poisson regression: We fit a Poisson regression model for each gene using the tissue type (tumour or normal) and library size as predictors. The p-value for the tissue variable is extracted to test if the gene is differentially expressed between tumour and normal samples.
Plot: The log10(p-value) for each gene is plotted. Red dots highlight genes with p-values < 0.05, indicating significant differential expression.


# Initialize a vector to store p-values from the model without the tissue covariate
p_values_no_tissue <- numeric(nrow(cleaned_Y))

# Loop over each gene and fit the Poisson model WITHOUT tissue type as a covariate
for (i in 1:nrow(cleaned_Y)) {
  y <- cleaned_Y[i, ]
  
  tmp <- data.frame(
    y = y,
    patient = cleaned_patient_data$patient,
    lib_size = colSums(cleaned_Y)
  )
  
  model_no_tissue <- glm(y ~ lib_size, data = tmp, family = "poisson")
  
  # Extract the p-value for the model without tissue
  p_values_no_tissue[i] <- summary(model_no_tissue)$coefficients[2, 4]
}

# Calculate the -log10(p-values) without tissue
log10_p_values_no_tissue <- -log10(p_values_no_tissue)

# Plot both p-value sets together for comparison
plot(log10_p_values, pch = 20, col = "blue", 
     main = "Comparison of Log10(p-values) With and Without Tissue Covariate",
     xlab = "Gene Index", ylab = "-log10(p-value)", ylim = range(log10_p_values, log10_p_values_no_tissue))

# Add the p-values from the model without the tissue covariate
points(log10_p_values_no_tissue, pch = 20, col = "green")

# Add a legend
legend("topright", legend = c("With Tissue", "Without Tissue"), col = c("blue", "green"), pch = 20)

The comparison plot shows the p-values from both models. Genes with larger differences between the two models will indicate that the tissue covariate had a significant impact on the gene's expression levels.

# Significant genes with tissue covariate
sig_with_tissue <- sum(p_values < 0.05)

# Significant genes without tissue covariate
sig_without_tissue <- sum(p_values_no_tissue < 0.05)

# Print out the number of significant genes
print(paste("Number of significant genes with tissue covariate: ", sig_with_tissue))
print(paste("Number of significant genes without tissue covariate: ", sig_without_tissue))

The inclusion of the tissue covariate is likely to have a large effect on the p-values, as the tissue type (tumour vs. normal) is expected to be a key driver of differential gene expression.

The library size (adjusting for sequencing depth) has a smaller impact compared to tissue type. However, it is important for normalizing the data to account for technical differences.

Based on the results, tissue type should have the largest effect on gene expression. This is because we expect large differences in gene expression between tumour and normal samples, and the p-values should be smaller (more significant) when tissue type is included.
More genes are likely to be identified as differentially expressed between tumour and normal samples.
The model ignores a major biological difference, leading to fewer significant results and less reliable differential expression analysis.

