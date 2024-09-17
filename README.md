# Multi-Omics Integration and Analysis for Environmental Study

## Overview
This project explores the impact of polar organic compounds on *Daphnia magna* through multi-omics analysis (RNA-Seq, metabolomics) and machine learning approaches. The analysis integrates data from environmental chemical concentrations and bioassays to draw biological insights.

## Datasets
The following datasets are used in this study:
- `water_chemicals.tsv`: Semi-quantified screening of 91 chemicals in 12 river water samples.
- `rna_norm_counts.csv`: RNA-Seq normalized read counts.
- `polar_pos_pqn_imputed.csv`: Metabolome data (positive mode) after imputation and normalization.
- `polar_neg_pqn_imputed.csv`: Metabolome data (negative mode) after imputation and normalization.
- `sample_sheet.csv`: Sample metadata with concentration levels.

## Key Steps
1. **Data Preprocessing**: Impute missing values, normalize RNA-Seq and metabolome data.
2. **Exploratory Data Analysis**: Perform PCA to visualize variance in data.
3. **Multi-Omics Integration**: Combine RNA-Seq and metabolome data for joint analysis.
4. **Machine Learning**: Train a Random Forest classifier to predict concentration levels.
5. **Network Biology**: Build co-expression networks to identify key gene-gene interactions.
6. **Pathway Enrichment**: Use external tools to perform pathway enrichment analysis based on gene expression data.

## Getting Started
### Prerequisites
Install the required packages by running:
```bash
pip install -r requirements.txt
