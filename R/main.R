#' Co-occurrence feature selection
#'
#' Search for seed genes in the gene expression data
#' Selecting features that form significant correlations with seed genes 
#'
# @param mtx raw UMI count gene expression matrix
# @param n_seeds number of seeds genes to use in feature selection (estimate would be # cell types)
# @param lower_bound b/w 0 and 1; search for seeds in genes detected in > lower_bound fraction of cells 
# @param uppe_bound b/w 0 and 1; search for seeds in genes detected in < upper_bound fraction of cells 
# @return dataframe consisting of significant correlations between genes
feature.selection <- function(
mtx, 
n_seeds=5, 
lower_bound=0.15, 
upper_bound=0.85
){
require(tidyverse)
require(philentropy)
## convert count data to binary 
mtx[mtx > 0] <- 1
## filter any undetected genes
mtx <- mtx[rowSums(mtx) > 1,]
seeds <- seed.find(mtx, n_seeds, lower_bound, upper_bound)
return(bind_rows(lapply(seeds, FUN=correlation.calc, mtx=mtx)))
}
