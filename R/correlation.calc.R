
	correlation.calc <- function(marker, mtx){
		## occurrences
		freq.df <- as.data.frame(rowSums(mtx)) %>% 
					rownames_to_column %>% 
					setNames(c('G2', 'G2.occur')) %>% 
					filter(G2 != marker)
		freq.df[,'G1'] <- marker
		freq.df[,'G1.occur'] <- sum(mtx[marker,])
		freq.df[,'cooccur'] <- apply(freq.df, MARGIN=1, FUN=obs.cooccur, mtx=mtx)
		freq.df[is.na(freq.df$cooccur),'cooccur'] <- 0 
		## calculate p value 
		freq.df$res <- apply(freq.df, MARGIN=1, FUN=new.pvalue, mtx=mtx)
		freq.df <- freq.df %>% separate(col='res', into=c('p_lt', 'p_gt'), sep='_') %>% 
					mutate(p_lt=as.numeric(p_lt)) %>% 
					mutate(p_gt=as.numeric(p_gt))
		## MTC 
		freq.df$p_lt_adj <- p.adjust(p=freq.df$p_lt, method='bonferroni', n=nrow(freq.df))
		freq.df$p_gt_adj <- p.adjust(p=freq.df$p_gt, method='bonferroni', n=nrow(freq.df))
		freq.df <- freq.df[freq.df$p_gt_adj < 0.05 | freq.df$p_lt_adj < 0.05,]
		return(freq.df %>% select(G1, G2, G1.occur, G2.occur, cooccur, p_lt, p_gt, p_lt_adj, p_gt_adj))
	}
	
	obs.cooccur <- function(x, mtx){
		return(table(colSums(mtx[x[c('G1', 'G2')],])==2)['TRUE'])
	}
	
	
