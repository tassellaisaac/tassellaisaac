This project implements a Differentially Expressed Genes (DEG) and KEGG Pathway Enrichment Analysis between patient samples that are sensitive and resistant to radiotherapy. The RData file provided contains the following datasets:

Expression matrix: Gene expression data for patients.
Sample information: Details about whether the samples are sensitive or resistant to radiotherapy.
Gene information: Gene annotations.

The project aims to:

Identify differentially expressed genes (DEGs) between the sensitive and resistant samples and display the top 15 genes.
Perform KEGG pathway enrichment analysis to identify upregulated and downregulated pathways in sensitive vs resistant samples.
Present the results of the DEGs and pathways identified, highlighting the most important genes and pathways.

You will need the following R packages for this analysis:

limma: For differential gene expression analysis.
edgeR: For RNA-Seq data normalization and statistical analysis.
clusterProfiler: For KEGG pathway enrichment analysis.
org.Hs.eg.db: For gene annotation (if working with human data).
ggplot2: For visualization.
biomaRt: To annotate the genes.

You can install the required packages using the following R commands:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "edgeR", "clusterProfiler", "org.Hs.eg.db", "ggplot2", "biomaRt"))

# Load necessary libraries
library(limma)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(biomaRt)

# Load the RData file
load("expression_data.RData")


# expression_matrix: Gene expression data
# sample_info: Information on whether samples are sensitive or resistant to radiotherapy
# gene_info: Gene annotations

# Prepare the DGEList object for analysis
dge <- DGEList(counts = expression_matrix)

# Normalization
dge <- calcNormFactors(dge)

# Create the design matrix based on sample information
design <- model.matrix(~ sample_info$status)
dge <- estimateDisp(dge, design)

# Fit the linear model using limma-voom
v <- voom(dge, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Get top differentially expressed genes (DEGs)
deg_results <- topTable(fit, coef = 2, n = Inf, sort.by = "P")
top_15_genes <- head(deg_results, 15)

# Display the top 15 differentially expressed genes
print("Top 15 Differentially Expressed Genes:")
print(top_15_genes)

# Save the top 15 DEGs to a CSV file
write.csv(top_15_genes, "top_15_degs.csv")

# Perform KEGG Pathway Enrichment Analysis
# Convert gene symbols to Entrez IDs
deg_entrez <- bitr(rownames(deg_results), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform enrichment analysis using clusterProfiler
kegg_enrichment <- enrichKEGG(gene = deg_entrez$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)

# Display the enriched KEGG pathways
print("KEGG Pathways Enriched in DEGs:")
print(kegg_enrichment)

# Save KEGG enrichment results to a CSV file
write.csv(kegg_enrichment, "kegg_enrichment_results.csv")

# Visualize top KEGG pathways
barplot(kegg_enrichment, showCategory = 10, title = "Top 10 Enriched KEGG Pathways")

# Identify important genes from the analysis
# Filter DEGs by logFC and p-value to identify significant up/down-regulated genes
significant_genes <- subset(deg_results, abs(logFC) > 1 & adj.P.Val < 0.05)
print("Significant Up/Down-regulated Genes:")
print(significant_genes)

# Save the significant genes to a CSV file
write.csv(significant_genes, "significant_genes.csv") 

Results and Interpretation
1. Top Differentially Expressed Genes (DEGs):
The top 15 DEGs between sensitive and resistant radiotherapy samples will be displayed. These genes are critical in understanding the molecular basis of sensitivity or resistance to radiotherapy.

2. KEGG Pathway Enrichment:
The KEGG pathway enrichment analysis identifies pathways that are significantly upregulated or downregulated in either the sensitive or resistant samples. This helps in understanding the broader biological processes involved.

Important Notes
Ensure the data is properly formatted and loaded correctly from the RData file.
The analysis assumes that gene symbols are used in the expression matrix and converts them to Entrez IDs for KEGG analysis.
Modify the filtering criteria (logFC, p-value) based on your specific requirements when identifying significant genes.

Further Extensions
The script can be extended to include GO term enrichment analysis using enrichGO() from the clusterProfiler package.
Consider visualizing the expression levels of the top DEGs using heatmaps or volcano plots.

