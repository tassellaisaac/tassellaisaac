This project demonstrates how to perform differential gene expression analysis using the DESeq2 R package on RNA-Seq data. The dataset used for this analysis comes from Drosophila melanogaster (Pasilla dataset), where samples are either treated or untreated. The goal of the project is to identify differentially expressed genes (DEGs) between these two conditions.

Objectives
Differential Expression Analysis:
Identify differentially expressed genes between treated and untreated conditions using DESeq2.
Perform normalization of the RNA-seq data to account for different sequencing depths.
Generate Visualizations:
Heatmap of highly expressed genes and differentially expressed genes.
PCA plot for visualization of sample grouping.
Correlation matrix of the samples.

The dataset includes:

Count Data: RNA-Seq gene expression data for treated and untreated Drosophila melanogaster samples.
Samples are provided in a count matrix format, where:
Rows represent genes.
Columns represent different samples.
Sample Information: Condition metadata (treated/untreated) and sequencing library type (paired-end/single-end) for the samples.


This analysis requires the following R packages:

DESeq2: For differential expression analysis.
gplots: For heatmap visualization.
RColorBrewer: For setting color palettes in plots.

You can install the required packages by running the following commands in R:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
install.packages("gplots")
install.packages("RColorBrewer")


Load the necessary R packages for analysis:
library("DESeq2")
library("gplots")
library('RColorBrewer')


Read the count data file containing gene expression data into R and assign it to an object (count.table):
count.table <- read.table("pasilla_gene_counts.tsv", header=TRUE, row.names=1, quote="", comment.char="")
head(count.table)

Define the experimental condition (treated vs untreated) and the type of sequencing library used (single-end vs paired-end):
cond.type <-  c( "untreated", "untreated", "untreated","untreated", "treated", "treated", "treated" )
lib.type  <-  c( "single-end", "single-end", "paired-end", "paired-end", "single-end", "paired-end", "paired-end" )


Filter the count data to include only paired-end samples:
isPaired <-  lib.type == "paired-end"
count.table <-  count.table[ , isPaired ]
head(count.table)


Create a DESeq2 dataset and run the differential expression analysis:
samples <- data.frame(row.names=c("untreated3", "untreated4", "treated2", "treated3"),
                      condition=as.factor(c(rep("untreated",2),rep("treated",2))))
pasCDS <- DESeqDataSetFromMatrix(countData = count.table, colData=samples, design=~condition)
pasCDS$condition <- factor(pasCDS$condition, levels=c("untreated", "treated"))
pasCDS <- DESeq(pasCDS)

# Extract the results
res <- results(pasCDS)
resOrdered <- res[order(res$padj), ]


Save the normalized counts and differentially expressed genes (DEGs) to CSV files:
write.csv(as.data.frame(resOrdered), file="condition_treated_results.csv")
write.csv(as.data.frame(counts(pasCDS, normalized=TRUE)), file="normalised_condition_treated.csv")


Filter for significant differentially expressed genes (adjusted p-value < 0.1):
resSig <- res[which(res$padj < 0.1), ]
head(resSig[order(resSig$log2FoldChange), ])  # Top down-regulated genes
tail(resSig[order(resSig$log2FoldChange), ])  # Top up-regulated genes

Heatmap of Top 30 Highly Expressed Genes:
select <- order(rowMeans(counts(pasCDS, normalized=TRUE)), decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(counts(pasCDS, normalized=TRUE)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale='none', 
          dendrogram='none', trace='none', margin=c(10,6))

select <- rownames(counts(pasCDS, normalized=TRUE)) %in% rownames(resOrdered[1:50,])
heatmap.2(counts(pasCDS, normalized=TRUE)[select,], col = hmcol, scale='row', trace='none', margin=c(10,6))

rld <- rlogTransformation(pasCDS, blind=TRUE)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
heatmap.2(mat, Rowv = as.dendrogram(hclust(distsRL)), symm = TRUE, trace='none', col = rev(hmcol), margin=c(13, 13))

plotPCA(rld, intgroup=c('condition'))

Expected Outputs
Differential Gene Expression Results: The top differentially expressed genes will be saved to condition_treated_results.csv.
Normalized Counts: Normalized counts across all samples will be saved to normalised_condition_treated.csv.
Heatmap: Heatmaps for the top 30 highly expressed genes and top 50 differentially expressed genes.
PCA Plot: A PCA plot that shows the clustering of samples based on their gene expression profiles.
Correlation Matrix: Heatmap showing the correlation between different samples.

Important Considerations
Filtering for Low Count Genes: It's essential to pre-filter the low-count genes to avoid noise in the analysis.
Normalization: The DESeq2 method performs normalization internally, which is crucial to ensure the counts are comparable between samples with different sequencing depths.
Fold Change: The log2 fold change values help in identifying genes that are significantly upregulated or downregulated between the conditions.
