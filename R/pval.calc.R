
	pvalue.calc <- function(x, mtx){
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
