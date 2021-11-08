# README
Feature selection using co-occurrence correlation improves cell clustering and embedding in single cell RNAseq data

# Figure 1 
* Fig. 1. Outlier UMI counts inflate standardized variance in PBMC 3k data. (a) Non zero UMI count distributions of DPH6 (HVG; red) and FCGR1B (LVF, teal). Both genes were detected in 61 out of 2700 cells. (b) Distribution of standardized variances. Colored lines indicate the value of FCGRB1 (teal) and DHP6 (red). (c) Standardized variances of 2,000 HVGs and 2,000 randomly sampled LVGs before and after converting the maximum UMI count to the median non zero UMI count. Red line indicates the minimum standardized variance of the original 2,000 HVGs.

# Figure 2
* Fig. 2. UMAP reductions comparing the performance of VST- and co-occurrenceselected features. Reductions were generated using the intersection of VST and cooccurrence features (top, left), features exclusive to co-occurrence (top, right),
features exclusive to VST (bottom, left), or 2,000 randomly selected and previously
unused features (bottom, right).

# Figure 3
* Fig. 3. Concordance between cell clusters generated using either the co-occurrence or VST derived features and known cell types in FACS-sorted PBMC datasets. Numbers represent the percentage of cells from a cell type assigned to the respective cell cluster. Cell clusters were generated using variable and co-occurrence features selected using the PBMC 3k dataset.

# Figure 4 
* Fig. 4. Dimension reduction usingco-occurrence derived features improves cell clustering as measured by adjusted Rand Index and silhoutte coefficient in simulated sc-RNAseq data. (a) UMAP reductions of a simulated gene expression matrix with a de.prob equal to 0.3 and four cell groups. Reductions were generated using the cooccurrence (co-occur) (left) or VST (right) derived features. Cells are colored based on their ground-truth, simulated cell type. (b) Average difference in silhoutte coefficient and ARI between co-occurrence and VST derived features based cell reductions and clusterings. Dimension reduction and cell clustering were performed with Seurat’s UMAP and Louvain algorithms with the first 10 PCs and a resolution
of 1.0.

# Method
* Bullet list
  * Nested bullet
    * Sub-nested bullet etc
  * Bullet list item 2
# Example
