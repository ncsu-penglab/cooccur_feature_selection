
	feature.selection <- function(mtx, n_genes, lower_bound, upper_bound){
		require(tidyverse)
		require(philentropy)
		# mtx <- as.data.frame(mtx)
		mtx[mtx > 0] <- 1
		mtx <- mtx[rowSums(mtx) > 1,]
		seeds <- seed.find(mtx, n_genes, lower_bound, upper_bound)
		return(bind_rows(lapply(seeds, FUN=correlation.calc, mtx=mtx)))
	}
