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

# data_preprocessing.py
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer

def load_and_impute_data(file_path):
    """Load dataset and impute missing values using k-NN."""
    data = pd.read_csv(file_path)
    imputer = KNNImputer(n_neighbors=5)
    imputed_data = pd.DataFrame(imputer.fit_transform(data), columns=data.columns)
    return imputed_data

def scale_data(data):
    """Scale the data using StandardScaler."""
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data)
    return scaled_data


    # Impute and scale metabolome data
    polar_pos_data = load_and_impute_data('../data/polar_pos_pqn_imputed.csv')
    polar_neg_data = load_and_impute_data('../data/polar_neg_pqn_imputed.csv')
    
    # Scaling
    scaled_polar_pos = scale_data(polar_pos_data)
    scaled_polar_neg = scale_data(polar_neg_data)
    
    # Save to file 
    pd.DataFrame(scaled_polar_pos).to_csv('../data/scaled_polar_pos.csv', index=False)
    pd.DataFrame(scaled_polar_neg).to_csv('../data/scaled_polar_neg.csv', index=False)

# exploratory_analysis.py
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def perform_pca(data, title, output_path):
    """Perform PCA and save the plot to the output directory."""
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(data)
    
    plt.figure(figsize=(8,6))
    plt.scatter(pca_result[:,0], pca_result[:,1], c='blue', alpha=0.5)
    plt.title(title)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)")
    plt.savefig(output_path)
    plt.show()


    # Load the RNA-Seq data
    rna_data = pd.read_csv('../data/rna_vst_counts.csv')

    # Scale the RNA data
    scaler = StandardScaler()
    scaled_rna = scaler.fit_transform(rna_data.T)
    
    # Perform PCA on RNA-Seq data
    perform_pca(scaled_rna, 'PCA of RNA-Seq Data', '../figures/pca_rna_seq.png')
    
    # Similarly, load and process metabolome data (positive mode)
    polar_pos_data = pd.read_csv('../data/scaled_polar_pos.csv')
    perform_pca(polar_pos_data, 'PCA of Metabolome Data (Positive Mode)', '../figures/pca_polar_pos.png')

# machine_learning.py
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

def train_random_forest(X_train, y_train):
    """Train a Random Forest classifier."""
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    rf.fit(X_train, y_train)
    return rf

def evaluate_model(rf, X_test, y_test):
    """Evaluate the model and return predictions."""
    y_pred = rf.predict(X_test)
    print(confusion_matrix(y_test, y_pred))
    print(classification_report(y_test, y_pred))
    return y_pred

def plot_feature_importance(rf, X_columns, output_path):
    """Plot and save feature importance."""
    importances = rf.feature_importances_
    indices = importances.argsort()[-10:]  # Top 10 features
    plt.figure(figsize=(10,6))
    plt.title("Top 10 Feature Importances")
    plt.barh(range(10), importances[indices], color='green', align='center')
    plt.yticks(range(10), [X_columns[i] for i in indices])
    plt.xlabel('Relative Importance')
    plt.savefig(output_path)
    plt.show()


    # Load data (RNA + Metabolomics)
    rna_data = pd.read_csv('../data/rna_vst_counts.csv').T
    polar_pos_data = pd.read_csv('../data/scaled_polar_pos.csv').T
    X = pd.concat([rna_data, polar_pos_data], axis=1)
    
    # Load labels
    y = pd.read_csv('../data/sample_sheet.csv')['concentration'].map({'Control': 0, '1x': 1, '10x': 2})

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Train model
    rf = train_random_forest(X_train, y_train)

    # Evaluate model
    evaluate_model(rf, X_test, y_test)

    # Plot feature importance
    plot_feature_importance(rf, X.columns, '../figures/feature_importance.png')

# network_biology.py
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def build_coexpression_network(data, threshold=0.9):
    """Build co-expression network from correlation matrix."""
    correlation_matrix = np.corrcoef(data.T)
    G = nx.Graph()
    for i in range(correlation_matrix.shape[0]):
        for j in range(i+1, correlation_matrix.shape[1]):
            if correlation_matrix[i, j] > threshold:
                G.add_edge(i, j)
    return G

def plot_network(G, output_path):
    """Plot and save the network graph."""
    plt.figure(figsize=(10, 10))
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=False, node_size=10, edge_color='b', alpha=0.3)
    plt.savefig(output_path)
    plt.show()


    # Load RNA-Seq data
    rna_data = pd.read_csv('../data/rna_vst_counts.csv')

    # Build co-expression network
    G = build_coexpression_network(rna_data)

    # Plot the network
    plot_network(G, '../figures/coexpression_network.png')

