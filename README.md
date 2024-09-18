This project aims to provide hands-on practice with RNA-Seq data analysis using the DESeq2 R package. DESeq2 is a powerful tool to identify differentially expressed genes (DEGs) between different conditions (in this case, samples that are sensitive or resistant to radiotherapy). Additionally, we generate visualizations such as heatmaps, PCA plots, and other graphical representations to help interpret the results.

The dataset used in this project is from the Pasilla RNA-Seq experiment on Drosophila melanogaster. The data was collected to investigate the effect of RNA interference on mRNA encoding RNA-binding proteins.

Objectives
Identify Differentially Expressed Genes (DEGs): Compare gene expression between samples sensitive and resistant to radiotherapy using DESeq2.
Normalization: Normalize read counts to account for differences in sequencing depth.
Visualizations:
Heatmaps of top DEGs.
Principal Component Analysis (PCA) plot for sample clustering.
Correlation matrix for evaluating sample similarities.


This analysis will require the following R packages:

DESeq2: For differential gene expression analysis.
ggplot2: For data visualization.
pheatmap: To generate heatmaps.
RColorBrewer: For color palettes in heatmaps.


You can install these packages using the following commands:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "pheatmap", "RColorBrewer"))
install.packages("ggplot2")

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Load the data
load("pasilla_data.RData") 

# Check if the data is loaded correctly
head(countData)
head(colData)

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

# Prefilter: Remove rows with 0 or 1 read
dds <- dds[rowSums(counts(dds)) > 1, ]

# Normalize the counts
dds <- DESeq(dds)

# Get the results of DESeq2 analysis
res <- results(dds)

# Order results by adjusted p-value
resOrdered <- res[order(res$padj), ]

# View the top differentially expressed genes
topGenes <- head(resOrdered, 15)
print("Top 15 Differentially Expressed Genes:")
print(topGenes)

# Save the top genes to a CSV file
write.csv(topGenes, "Top_15_Differentially_Expressed_Genes.csv")

# Log transform for visualization
rld <- rlog(dds, blind = FALSE)

# Generate a heatmap of the top 50 differentially expressed genes
top50genes <- head(order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE), 50)
mat <- assay(rld)[top50genes, ]

# Plot the heatmap
pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = TRUE, show_colnames = TRUE,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Top 50 Differentially Expressed Genes Heatmap")

# Save the heatmap to a file
pdf("Top_50_DEG_Heatmap.pdf")
pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = TRUE, show_colnames = TRUE,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Top 50 Differentially Expressed Genes Heatmap")
dev.off()

# PCA plot
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Generate PCA plot
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

# Save PCA plot to PDF
ggsave("PCA_plot.pdf")

# Generate correlation matrix heatmap
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colData$condition, colnames(countData), sep = "-")
colnames(sampleDistMatrix) <- NULL

# Plot the correlation heatmap
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample Correlation Matrix")

# Save correlation heatmap to PDF
pdf("Sample_Correlation_Heatmap.pdf")
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample Correlation Matrix")
dev.off()

Expected Outputs
Top 15 DEGs: The top 15 differentially expressed genes will be saved in a CSV file (Top_15_Differentially_Expressed_Genes.csv).
Heatmap: A PDF file containing the heatmap of the top 50 differentially expressed genes (Top_50_DEG_Heatmap.pdf).
PCA Plot: A PDF file showing the PCA plot for sample clustering (PCA_plot.pdf).
Sample Correlation Heatmap: A PDF file containing the correlation matrix of the samples (Sample_Correlation_Heatmap.pdf).

Important Considerations
Data Preprocessing: Ensure the RNA-Seq count matrix and metadata are correctly formatted before running the analysis.
Saving Images: You can adapt the code to save images in other formats (e.g., PNG, TIFF) by modifying the ggsave() function.
DESeq2 Normalization: The normalization performed by DESeq2 accounts for sequencing depth and ensures that read counts are comparable across samples.
