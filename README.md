This project performs Gene Set Enrichment Analysis (GSEA) using differential gene expression data obtained from an RNA-seq experiment. The analysis is conducted using the RStudio Server on the BEAR Portal. We will use MSigDB gene sets to identify pathways enriched in the differentially expressed genes between two conditions. 

The following R libraries are required for the analysis:

biomaRt: For gene ID annotation.
fgsea: For GSEA analysis.
data.table: For efficient data manipulation.
BiocParallel: For parallel processing (optional).
tidyverse: For data cleaning and visualization.


# For gene ID annotation
library(biomaRt)

# For GSEA analysis
library(fgsea)
library(data.table)
library(BiocParallel)
library(tidyverse)

help(package="biomaRt")

# Load the differential expression result
result_diffExp <- read.table("gtex.brain.2Tiss.diff-exp.outcome.txt", sep="\t", header=TRUE, row.names=1)

# Check the first few rows
head(result_diffExp, 3)

# Check dimensions of the dataset
dim(result_diffExp)

# Find how many gene IDs were significantly differentially expressed
length(which(result_diffExp$adj.P.Val < 0.05))

# Create a BioMart object
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get gene annotations for differential expression dataset
myAnno <- getBM(attributes=c("ensembl_gene_id","chromosome_name","hgnc_symbol"), 
                filters="ensembl_gene_id", values=rownames(result_diffExp), mart=ensembl)

# Check the first few rows of the annotation
head(myAnno, 3)

# Remove duplicate entries and empty gene symbols
DupRow <- rownames(myAnno[which(duplicated(myAnno$ensembl_gene_id)),])
myAnno_clean <- myAnno[-as.integer(DupRow), ]
myAnno_clean <- myAnno_clean[which(!(myAnno_clean$hgnc_symbol == "")), ]

# Merge annotation data with differential expression results
result_diffExp.anno <- merge(result_diffExp, myAnno_clean, by.x=0, by.y="ensembl_gene_id")

# Save the annotated result for future reference
write.table(result_diffExp.anno, "Hypoth.vs.Hippocmps_limma_reslt.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Order the dataset by t-statistic
result_diffExp.anno.ord <- result_diffExp.anno[order(result_diffExp.anno$t), ]

# Create a ranked list using the t-statistic values
gene.rnk <- result_diffExp.anno.ord$t
names(gene.rnk) <- result_diffExp.anno.ord$hgnc_symbol

# Load MSigDB Hallmark pathways
MSig.Hallmark <- gmtPathways("h.all.v7.5.1.symbols.gmt")

# Run GSEA
gene.rnk.fgsea <- fgseaMultilevel(pathways=MSig.Hallmark, stats=gene.rnk, eps=0.0, minSize=15, maxSize=500)

# Filter for significant pathways
gene.rnk.fgsea.sig <- gene.rnk.fgsea %>% filter(padj < 0.001 & (NES >= 1.5 | NES <= -1.5))

# Save the GSEA output
fwrite(gene.rnk.fgsea, file="Hypoth.vs.Hippocmps.fGseaReslt.txt", sep="\t")

# Create enrichment plots for significant pathways
pdf("Hypoth.vs.Hippocmps.fGsea.EnrichmentPlt.pdf", width=3, height=3)
for(i in 1:nrow(gene.rnk.fgsea.sig)) {
  print(plotEnrichment(MSig.Hallmark[[gene.rnk.fgsea.sig$pathway[i]]], gene.rnk) + 
        labs(title=gene.rnk.fgsea.sig$pathway[i]) + 
        theme(plot.title = element_text(size=5), axis.text.x = element_text(size=3), 
              axis.text.y = element_text(size=3), axis.title.x = element_text(size=4), 
              axis.title.y = element_text(size=4)))
}
dev.off()

pdf("Hypoth.vs.Hippocmps.fGseaTablePlt.pdf", width=11, height=9)
plotGseaTable(MSig.Hallmark[gene.rnk.fgsea.sig$pathway], gene.rnk, gene.rnk.fgsea, gseaParam=0.5)
dev.off()

Output Files
Annotated Differential Expression Results: Hypoth.vs.Hippocmps_limma_reslt.txt
GSEA Results: Hypoth.vs.Hippocmps.fGseaReslt.txt
Enrichment Plots PDF: Hypoth.vs.Hippocmps.fGsea.EnrichmentPlt.pdf
GSEA Summary Table PDF: Hypoth.vs.Hippocmps.fGseaTablePlt.pdf
