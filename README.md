# Predicting Diabetic Retinopathy from Peripheral Blood Gene Expression

**RNA-Seq analysis and machine learning pipeline for classifying diabetic retinopathy using leukocyte mRNA expression**

[![Dataset](https://img.shields.io/badge/GEO-GSE221521-orange)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221521)
![R](https://img.shields.io/badge/R-4.0+-blue)

## Overview

This project analyzes peripheral blood mRNA expression data from 193 participants to identify molecular signatures distinguishing diabetic patients with retinopathy (DR) from those without (DM). Using RNA-Seq data from leukocytes, we apply differential expression analysis, supervised learning, and dimension reduction to:

- Identify genes differentially expressed in diabetic retinopathy
- Build predictive models for DR classification (>70% accuracy)
- Discover biological pathways associated with retinopathy development
- Explore dimensionality reduction techniques for patient stratification

**Key Finding:** Gene expression from peripheral blood can predict diabetic retinopathy status with reasonable accuracy, suggesting potential biomarkers for non-invasive screening.

## Background

**Diabetic Retinopathy (DR)** is a complication of diabetes affecting the eyes, caused by damage to blood vessels in the retina due to prolonged high glucose levels. It's a leading cause of blindness in working-age adults.

**Research Question:** Can we identify molecular signatures in peripheral blood that distinguish diabetic patients with retinopathy from those without?

**Clinical Relevance:** Non-invasive blood-based biomarkers could enable earlier detection and monitoring of retinopathy risk without requiring specialized ophthalmologic examination.

## Dataset

**Source:** [GEO Accession GSE221521](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221521)  
**Reference:** Xiang et al. (2023). "Changes and Related Factors of Blood CCN1 Levels in Diabetic Patients." *Frontiers in Endocrinology*.

**Cohort:**
- **N = 193 participants** (95 female, 98 male)
- **Groups:** 
  - Control (non-diabetic): 50 participants
  - DM (Diabetes Mellitus without retinopathy): 54 participants
  - DR (Diabetic Retinopathy): 89 participants
- **Platform:** Illumina NovaSeq 6000
- **Features:** 16,697 genes (after filtering from 60,675 total)

## Methods

### 1. Preprocessing & Quality Control
- **Low count filtering:** Removed genes with <10 average counts per sample
- **Normalization:** Variance Stabilizing Transformation (VST) using DESeq2
- **Subset analysis:** Focused on DM vs. DR comparison for classification tasks

### 2. Differential Expression Analysis
- **Model:** Negative Binomial generalized linear model (DESeq2)
- **Comparisons:** DR vs. DM, DR vs. Control
- **Significance criteria:** 
  - Adjusted p-value < 0.05 (Benjamini-Hochberg FDR correction)
  - |log₂ Fold Change| > 1 for high-confidence genes
- **Functional analysis:** Gene Ontology enrichment using clusterProfiler

**Key Results:**
- **3,829 differentially expressed genes** (DR vs. DM)
  - 2,320 upregulated, 1,509 downregulated
  - Enriched pathways: RNA splicing, ribonucleoprotein complex biogenesis, ncRNA processing
- **5,110 differentially expressed genes** (DR vs. Control)

### 3. Supervised Learning

#### Lasso Logistic Regression
- **Feature selection:** L1 regularization for sparse model
- **Cross-validation:** 10-fold CV to optimize lambda parameter
- **Train/Test split:** 70/30
- **Performance:** 67.14% test accuracy, AUC = 0.70
- **Selected features:** 35 genes with non-zero coefficients
- **Top genes:** TFIP11, AL356750.1, FRG1BP

#### Random Forest
- **Configuration:** 500 trees, default mtry
- **Out-of-bag error:** 23.29%
- **Performance:** 75.71% test accuracy, AUC = 0.74
- **Variable importance:** MED4, PLXND1, PIEZO1, NCL, AC073957.3
- **Class-specific error:** Better DM prediction (10.26% error) than DR (38.24% error)

### 4. Dimension Reduction

Explored three techniques for visualizing high-dimensional gene expression:

- **PCA:** Linear method; first 2 PCs explain 26.42% variance, showed substantial overlap
- **t-SNE:** Non-linear manifold learning; better separation with DR trending toward lower first component
- **UMAP:** Non-linear topology preservation; revealed distinct DR cluster, suggesting molecular subtypes

## Key Findings

### Biological Insights
- Differentially expressed genes are enriched in **RNA processing pathways**, suggesting altered post-transcriptional regulation in DR
- Genes like **RN7SK, RPS9, CEACAM6, OLR1** are downregulated in DR
- Upregulated genes include **ATP5PDP4, TATDN2P2, ACTBP2**

### Predictive Models
- **Random Forest outperforms Lasso** (75.71% vs. 67.14% accuracy)
- Both models achieve AUC > 0.70, indicating moderate discriminative ability
- **Challenge:** Higher misclassification rate for DR (false negatives), suggesting molecular heterogeneity

### Patient Stratification
- **UMAP reveals a distinct DR subgroup** with similar expression profiles
- Some DR patients cluster with DM patients, suggesting:
  - Early-stage DR with minimal molecular changes
  - Potential DR molecular subtypes
  - Heterogeneity in disease progression

## Technologies

**Language:** R 4.0+

**Core Packages:**
- `DESeq2`: Differential expression analysis
- `glmnet`: Lasso regression with regularization
- `randomForest`: Random forest classification
- `caret`: Model training and evaluation
- `tidyverse`: Data manipulation (dplyr, ggplot2)
- `clusterProfiler`: Gene ontology enrichment
- `org.Hs.eg.db`: Human gene annotation
- `Rtsne`: t-SNE implementation
- `umap`: UMAP dimension reduction

**Analysis Techniques:**
- RNA-Seq normalization and quality control
- Differential gene expression (Negative Binomial GLM)
- Multiple testing correction (Benjamini-Hochberg FDR)
- Regularized logistic regression
- Ensemble learning (Random Forest)
- Dimensionality reduction (PCA, t-SNE, UMAP)

### Data Access

Download the dataset from GEO:
```r
# Using GEOquery
library(GEOquery)
gse <- getGEO("GSE221521", GSEMatrix = TRUE)
```

Or download directly from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221521

## Reproducibility

- **Random seed:** Set to 2024 for train/test splits and cross-validation
- **Software versions:** R 4.0+, DESeq2 1.30+, glmnet 4.1+

## Clinical Implications

This analysis suggests that:

1. **Peripheral blood transcriptomics** may serve as a non-invasive screening tool for DR
2. **RNA processing pathways** warrant further investigation as therapeutic targets
3. **Patient heterogeneity** in DR may require personalized monitoring strategies
4. **Top predictive genes** (MED4, PLXND1, PIEZO1) could be validated as biomarker candidates

**Limitations:**
- Cross-sectional design; cannot establish causality
- Limited sample size for rare subtypes
- External validation needed on independent cohorts
- Blood-based signals may reflect systemic diabetes effects, not retina-specific changes

## Future Directions

- **Validation:** Test models on independent DR cohorts
- **Longitudinal analysis:** Track gene expression changes during DR progression
- **Integration:** Combine with clinical variables (HbA1c, duration of diabetes, blood pressure)
- **Single-cell RNA-Seq:** Identify specific leukocyte populations driving differential expression
- **Functional validation:** Experimental studies of candidate genes in retinal cell models

## References

**Dataset:**  
Xiang, Z.-Y., Chen, S.-L., Qin, X.-R., Lin, S.-L., Xu, Y., Lu, L.-N., & Zou, H.-D. (2023). Changes and Related Factors of Blood CCN1 Levels in Diabetic Patients. *Frontiers in Endocrinology*, 14. https://doi.org/10.3389/fendo.2023.1131993

**Methods:**  
Love, M.I., Huber, W. & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. https://doi.org/10.1186/s13059-014-0550-8

Yu, G., Wang, L., Han, Y., & He, Q. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS: A Journal of Integrative Biology*, 16(5), 284-287. https://doi.org/10.1089/omi.2011.0118

## Author

**Jaden Thomas**  
M.S. Biostatistics Student, University of Minnesota  
B.S. Data Analytics, The Ohio State University

Analysis completed as part of PUBH 5280 - Genomic Data Analysis coursework at The Ohio State University (Fall 2024)

## License

Code is available for educational and research purposes. Please cite the original dataset publication (Xiang et al., 2023) if using this data.

---

*For questions about the analysis or code, please open an issue or contact me via [your email/LinkedIn].*
