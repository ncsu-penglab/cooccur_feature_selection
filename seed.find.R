
	seed.find <- function(mtx, n_genes, lower_bound, upper_bound){
		t1 <- apply(mtx, MARGIN=1, FUN=sum) %>% as.data.frame %>%  rownames_to_column %>% 
			  setNames(c('gene', 'nCells')) %>% mutate(frac=nCells/ncol(mtx)) %>% 
			  filter(frac > lower_bound) %>% filter(frac < upper_bound)
		## jaccard distance calcuation
		t2 <- philentropy::distance(mtx[t1$gene,], method = "jaccard", use.row.names=TRUE) %>% as.data.frame %>% 
			  rownames_to_column(var='g1') %>% gather(key='g2', value='J.dist', -g1) 
		t2 <- t2[order(t2$J.dist, decreasing=TRUE),]
		## initialize seed list with most distant genes 
		seeds <- unlist(t2[1,,drop=TRUE])[names(unlist(t2[1,,drop=TRUE])) %in% c('g1','g2')]
		for (i in 1:(n_genes - length(seeds))){
			t3 <- t2[grep(paste(seeds, collapse='|'), t2$g1),] %>% group_by(g2) %>% summarize(mean.dist=mean(J.dist))
			t3 <- t3[order(t3$mean.dist, decreasing=TRUE),]
			seeds <- c(seeds, as.character(t3[1,'g2']))
			
		}
		return(seeds)
	}
	