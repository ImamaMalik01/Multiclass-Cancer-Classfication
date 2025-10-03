# Multiclass-Cancer-Classfication
Final Year Project: Classification of Colorectal, Skin, and Head and Neck cancers using ML/DL on gene expression (FPKM) data from GEO NCBI. Includes classification models, DGE analysis, heatmaps, and bioinformatics visualizations.
This project aims to classify Colorectal Cancer, Skin Cancer, and Head & Neck Cancer using Machine Learning (ML) and Deep Learning (DL) approaches applied to gene expression data (FPKM values) from the GEO NCBI database. The dataset contains 23,505 gene features per sample, obtained from high-throughput sequencing experiments.

The workflow includes:

Preprocessing gene expression data (FPKM normalized values).

Implementing ML models (XGBoosting, Random Forest) for binary and multiclass classification. A venn diagram for looking similar genes within the different cancer types. Using Shap value feature on Deep learning model.

Training a Feedforward Neural Network (FNN) for improved multiclass accuracy.

Performing Differential Gene Expression (DGE) analysis using R.

Visualizing results through heatmaps and Venn diagrams to identify significant upregulated and downregulated genes.

This project demonstrates how bioinformatics and AI can be combined to extract biological insights from sequencing data, identify potential biomarkers, and improve cancer classification accuracy.
