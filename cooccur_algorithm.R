	
###########################################################################################################
	# preliminary script for running co-occurrecne feature selection of sc-RNAseq gene expression matrix:
	#	Parameters:
	#		mtx: raw, unnormalized gene expression matrix 
	#		n_genes: number of seed genes used for feature selection (should be the estimated number of cells in population)
	#		lower_bound: to be considered as a seed, each gene must be detected in at least lower_bound percentage of cells 
	#		upper_bound: to be considered as a seed, each gene must be detected in no more than upper_bound percentage of cells 
###########################################################################################################

	main <- function(mtx, n_genes, lower_bound, upper_bound){
		require(tidyverse)
		require(philentropy)
		# mtx <- as.data.frame(mtx)
		mtx[mtx > 0] <- 1
		mtx <- mtx[rowSums(mtx) > 1,]
		seeds <- seed.find(mtx, n_genes, lower_bound, upper_bound)
		return(bind_rows(lapply(seeds, FUN=marker.fun, mtx=mtx)))
	}
	
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
	
	marker.fun <- function(marker, mtx){
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
	
	new.pvalue <- function(x, mtx){
		obs_cooccur <- as.numeric(x['cooccur'])
		feature.inc <- as.numeric(x['G2.occur'])
		marker.inc <- as.numeric(x['G1.occur'])
		max_inc <- max(feature.inc, marker.inc)
	    min_inc <- min(feature.inc, marker.inc)
		## number of possible cooccurences is the total number of cells 
		nsite <- ncol(mtx)
		all.probs <- phyper(q=0:min_inc, 
							m=min_inc, 
							n=(nsite - min_inc), 
							k=max_inc)
		
		psite <- as.numeric(nsite + 1)
		
		prob_share_site <- mapply(FUN=prob_share_site.fun, x=all.probs[2:length(all.probs)], xminus1=all.probs[1:length(all.probs)-1])
		prob_share_site <- c(all.probs[1], prob_share_site, rep(0, psite - length(prob_share_site) - 1))
		
		p_lt <- sum(prob_share_site[1:(obs_cooccur+1)])
		p_gt <- sum(prob_share_site[(obs_cooccur+1):length(prob_share_site)])
		res <- paste(p_lt, '_', p_gt)
		return(res)
	
	}
	
	prob_share_site.fun <- function(x, xminus1){
		res <- x - xminus1
		return(res)
	}
	

