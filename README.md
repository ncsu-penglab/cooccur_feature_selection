# README
Feature selection using co-occurrence correlation improves cell clustering and embedding in single cell RNAseq data

# Figure 1 
![alt text](https://github.com/ncsu-penglab/cooccur_feature_selection/blob/master/Results/figure1_1.jpg?raw=true) [Fig. 1. Outlier UMI counts inflate standardized variance in PBMC 3k data. (a) Non zero UMI count distributions of DPH6 (HVG; red) and FCGR1B (LVF, teal). Both genes were detected in 61 out of 2700 cells. (b) Distribution of standardized variances. Colored lines indicate the value of FCGRB1 (teal) and DHP6 (red). (c) Standardized variances of 2,000 HVGs and 2,000 randomly sampled LVGs before and after converting the maximum UMI count to the median non zero UMI count. Red line indicates the minimum standardized variance of the original 2,000 HVGs.]

# Figure 2
![alt text](https://github.com/ncsu-penglab/cooccur_feature_selection/blob/master/Results/figure2_1.jpg?raw=true)
* Fig. 2. UMAP reductions comparing the performance of VST- and co-occurrenceselected features. Reductions were generated using the intersection of VST and cooccurrence features (top, left), features exclusive to co-occurrence (top, right),
features exclusive to VST (bottom, left), or 2,000 randomly selected and previously
unused features (bottom, right).

# Figure 3
![alt text](https://github.com/ncsu-penglab/cooccur_feature_selection/blob/master/Results/figure3_1.jpg?raw=true)
* Fig. 3. Concordance between cell clusters generated using either the co-occurrence or VST derived features and known cell types in FACS-sorted PBMC datasets. Numbers represent the percentage of cells from a cell type assigned to the respective cell cluster. Cell clusters were generated using variable and co-occurrence features selected using the PBMC 3k dataset.

# Figure 4 
![alt text](https://github.com/ncsu-penglab/cooccur_feature_selection/blob/master/Results/figure4_1.jpg?raw=true)
* Fig. 4. Dimension reduction usingco-occurrence derived features improves cell clustering as measured by adjusted Rand Index and silhoutte coefficient in simulated sc-RNAseq data. (a) UMAP reductions of a simulated gene expression matrix with a de.prob equal to 0.3 and four cell groups. Reductions were generated using the cooccurrence (co-occur) (left) or VST (right) derived features. Cells are colored based on their ground-truth, simulated cell type. (b) Average difference in silhoutte coefficient and ARI between co-occurrence and VST derived features based cell reductions and clusterings. Dimension reduction and cell clustering were performed with Seurat’s UMAP and Louvain algorithms with the first 10 PCs and a resolution
of 1.0.

# Table I
| Features | Silhouette Width | ARI |
| :------: | :------: | :------: |
| VST ⋂ cooccur | 0.366±0.44 | 0.627±0.02 |
| cooccur ⊄ VST | 0.268±0.44 | 0.553±0.04 |
| VST ⊄ cooccur | -0.083±0.29 | 0.215±0.01 |
| SCT ⋂ cooccur | 0.392±0.45 | 0.668±0.03 |
| cooccur ⊄ SCT | 0.213±0.44 | 0.518±0.02 |
| SCT ⊄ cooccur | 0.02±0.32 | 0.333±0.07 |
| random | -0.141±0.24 | 0.186±0.03 |
* Table. I. Performance comparison of the co-occurrence (co-occur)
and VST or SCTransform (SCT) derived features in the PBMC 3k
dataset measured by silhoutte width and adjusted Rand Index (ARI).

# Table II
|             | Silhouette Width | Silhouette Width | ARI | ARI |
| :---------: | :--------------: |:--:|:---:|:--:|
| Dataset     | cooccur | VST | cooccur | VST|
| Zhengmix4eq | 0.728 ± 0.21 | 0.643 ± 0.3 | 0.96 | 0.691 |  
| Zhengmix4uneq | 0.669 ± 0.27 | 0.537 ± 0.41 | 0.962 | 0.833| 
| Zhengmix8eq | 0.453 ± 0.45 | 0.337 ± 0.5  | 0.653 | 0.522 |
* Table. II. Performance comparison of cell clustering using the co-occurrence (cooccur) vs. VST derived features
in the FACS PBMC datasets measured by silhoutte width and adjusted Rand Index (ARI).
