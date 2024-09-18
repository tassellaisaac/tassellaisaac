This project provides a step-by-step guide to analyzing single-cell RNA sequencing (scRNA-seq) data using the Seurat package in R. The data  is derived from Patient N1469, and the goal is to perform quality control, normalization, dimensionality reduction, cell clustering, and marker gene identification. Additionally, cell type annotation is conducted using a reference bulk RNA-seq dataset.

You will need the following R packages for the analysis:

Seurat: For handling single-cell RNA-seq data.
SingleCellExperiment: For working with single-cell data structures.
scater: For quality control and visualization.
SC3: For single-cell clustering.
ggplot2: For visualizations.
vcd: For producing ternary plots.

# Install packages if not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("SingleCellExperiment", "scater", "SC3"))
install.packages(c("Seurat", "ggplot2", "vcd"))

The data used for this analysis includes:

Raw scRNA-seq count data for Patient N1469, which is loaded as a 10X dataset using Seurat’s Read10X() function.
Human mammary gland epithelium reference RNA-seq signature genes from a bulk RNA-seq study (GSE161892), used for cell type annotation.

Data Structure
Rows: Genes/features.
Columns: Individual cells (with barcodes as identifiers).

library(Seurat)
library(SingleCellExperiment)
library(scater)
library(SC3)
library(ggplot2)
library(vcd)  # For ternary plots

# Load the 10X data for Patient N1469
N1469.data <- Read10X(data.dir = "Data/N1469")
colnames(N1469.data) <- paste("N1469", colnames(N1469.data), sep="_")

# Create a Seurat object
N1469 <- CreateSeuratObject(counts=N1469.data, project="N1469", min.cells=3, min.features=200)

#Quality Control
N1469[["percent.mt"]] <- PercentageFeatureSet(N1469, pattern = "^MT-")
N1469 <- subset(N1469, subset = nFeature_RNA > 500 & percent.mt < 20)

#Normalization 
N1469 <- NormalizeData(N1469)
N1469 <- FindVariableFeatures(N1469, selection.method="vst", nfeatures=1500)

#Scaling and Dimensionality Reduction
N1469 <- ScaleData(N1469)
N1469 <- RunPCA(N1469, features=VariableFeatures(N1469))
DimPlot(N1469, reduction = "pca")

N1469 <- RunTSNE(N1469, dims=1:30, seed.use=2021)
DimPlot(N1469, reduction = "tsne")

#Cell Clustering
N1469 <- FindNeighbors(N1469, dims=1:30)
N1469 <- FindClusters(N1469, resolution=0.1)
DimPlot(N1469, label = TRUE)

#Marker Gene Identification
cluster1.markers <- FindMarkers(N1469, ident.1 = 1, min.pct = 0.25)
N1469.markers <- FindAllMarkers(N1469, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Feature and Dot Plotting
FeaturePlot(N1469, features = top3$gene, ncol=3)
DotPlot(N1469, features = top3$gene, dot.scale = 8) + RotatedAxis()


#Cell Type Annotation
load("Data/Human-PosSigGenes.RData")
HumanSig <- list(Basal=Basal, LP=LP, ML=ML, Str=Str)

SigScore <- list()
for(i in 1:length(HumanSig)){
  SigScore[[i]] <- colMeans(N1469@assays$RNA@data[HumanSig[[i]], ])
}
SigScores <- do.call("cbind", SigScore)
colnames(SigScores) <- c("Basal", "LP", "ML", "Stroma")
N1469@meta.data <- cbind(N1469@meta.data, SigScores)

VlnPlot(N1469, features = colnames(SigScores), ncol=2, pt.size=0.1)


#Ternary Plot
TN <- matrix(0L, ncol(N1469), 3L)
colnames(TN) <- c("LP", "ML", "Basal")
for(i in colnames(TN)){
  TN[, i] <- colSums(N1469@assays$RNA@counts[HumanSig[[i]], ] > 0L)
}
ternaryplot(TN, cex=0.2, pch=16, col=col.p[Idents(N1469)], grid=TRUE)

Output Files
scRNAseqTODAY.RData: Save the analysis at various stages to avoid losing progress.
scRNAseq_Complete.RData: Final output after all steps have been performed, including clustering and marker identification.

