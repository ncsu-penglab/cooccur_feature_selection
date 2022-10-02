#/peng_1/peng_lab/tools/stat/miniconda3/envs/EW/bin/python
################################################################################################
#	@Author: Evan Walsh
#	@File: /peng_1/peng_lab/scripts/cooccur_algorithm_python/cooccur/1_cooccur_algorithm.py
#	
#	Objective: script the co-occurrence feature selection algorithm in python
#	
################################################################################################

import numpy as np
import pandas as pd 
import scanpy as sc
import itertools as it
from scipy.stats import fisher_exact
from scipy.spatial import distance
from sklearn.metrics import jaccard_score
from statsmodels.stats.multitest import multipletests
import os 

# os.environ['R_HOME'] = '/opt/R/4.1.0/bin/R'
# from rpy2.robjects.packages import importr
# from rpy2.robjects.vectors import FloatVector
# stats = importr('stats')

## I/O
root = "/peng_1/peng_lab/results/ongoing/cooccur_algorithm_python"

## 1) filter low and high abundant features 
## 2) calculate a jaccard distance matrix to determine seed genes 
def cooccur(mtx, num_of_seeds, lower_bound, upper_bound):
	## convert matrix to binary 
	mtx[mtx.columns] = mtx[mtx.columns].apply(pd.to_numeric)
	mtx[:] = np.where(mtx < 1, 0, 1)
	
	## remove features that are not present 
	mtx = mtx.loc[:, mtx.sum(axis=0) > 0]
	
	## find seed genes 
	seeds = seed_find(mtx, num_of_seeds, lower_bound, upper_bound)
	
	df_l = []
	for seed in seeds:
		print(seed)
		df_l.append(marker_function(marker=seed, mtx=mtx))
	
	dat = pd.concat(df_l)
	return dat

def seed_find(mtx, num_of_seeds, lower_bound, upper_bound):
	## caclulate the frequency of each gene across all cells 
	ncells = len(mtx.columns)
	t1 = mtx.sum(axis=1).to_frame()
	t1 = t1.rename(columns={0: 'Freq'})
	t1['rate'] = t1['Freq'].div(ncells)
	
	## filter features to only keep those within frequency bounds 
	t2 = t1[(t1['rate'] > lower_bound) & (t1['rate'] < upper_bound)] 
	features2keep = t2.index.tolist()
	mtx = mtx.filter(items = features2keep, axis=0)
	
	## calculate jaccard distance between features 
	sim_df = pd.DataFrame(np.ones((len(features2keep), len(features2keep))), index=features2keep, columns=features2keep)
	for col_pair in it.combinations(features2keep, 2):
		sim_df.loc[col_pair] = sim_df.loc[tuple(reversed(col_pair))] = jaccard_score(mtx.T[col_pair[0]], mtx.T[col_pair[1]])
	
	sim_df.index.name = 'G1'
	sim_df.reset_index(inplace=True)
	sim_df = pd.melt(sim_df, id_vars='G1', var_name='G2')
	
	## get seed genes 
	seeds = sim_df[sim_df['value'] == min(sim_df['value'])].iloc[:,0].values
	for i in range(num_of_seeds-len(seeds)):
		tmp = sim_df[sim_df['G1'].isin(seeds) & -sim_df['G2'].isin(seeds)].groupby(by=['G2']).mean()
		seed_new = tmp[tmp['value'] == min(tmp['value'])].index
		seeds = np.append(seeds, seed_new)
	
	return seeds

## 3) measure correlation between a marker gene and other genes (rows)
def marker_function(marker, mtx):
	## calculate the number of occurrences of each gene 
	
	# mtx.sum(axis=1)
	freq_df = pd.DataFrame({"G2occur":mtx.sum(axis=1)})
	
	## convert rownames to a column
	freq_df.index.name = 'G2'
	freq_df.reset_index(inplace=True)
	
	## remove the G1 value from freq_df 
	G1_freq = freq_df[freq_df['G2']==marker].G2occur
	freq_df = freq_df[freq_df.G2 != marker]
	
	freq_df['G1'] = marker
	# freq_df['G1occur'] = G1_freq[0]
	freq_df['G1occur'] = G1_freq.values[0]
	
	## reorder columns 
	freq_df = freq_df[['G1','G2', 'G1occur','G2occur']]
	
	## caclulate the number of observed co-occurrence
	freq_df['cooccur'] = freq_df.apply(lambda row : observed_cooccur(row['G1'], row['G2'], mtx), axis = 1)
	
	mtx_t = mtx.T
	## calculate probability of co-occur less than what was observed 
	freq_df['p_lt'] = freq_df.apply(lambda row: contigency_tab(row['G1'], row['G2'], mtx_t=mtx_t, alt='less'), axis=1)
	
	## calculate probability of co-occur greater than what was observed 
	freq_df['p_gt'] = freq_df.apply(lambda row: contigency_tab(row['G1'], row['G2'], mtx_t=mtx_t, alt='greater'), axis=1)
	
	## calculate adjusted p-values 
	freq_df['p_lt_adj'] = multipletests(freq_df['p_lt'], method='bonferroni')[1].tolist()
	freq_df['p_gt_adj'] = multipletests(freq_df['p_gt'], method='bonferroni')[1].tolist()
	
	return freq_df

def observed_cooccur(G1, G2, mtx):
	s_obs = mtx.loc[[G1,G2]].sum(axis=0)
	s_obs_freq = s_obs.value_counts()
	
	if s_obs_freq.size < 3: 
		miss_vals = set([0,1,2]) - set(s_obs_freq.index)
		for i in miss_vals:
			# s_miss_val = pd.Series(0, index=str(i))
			s_miss_val = pd.Series(0, index=[i])
			s_obs_freq = s_obs_freq.append(s_miss_val)
	
	return s_obs_freq[2]

def contigency_tab(G1, G2, mtx_t, alt):
	## transform gene expression matrix 
	tmpo = pd.crosstab(index=mtx_t[G1], columns=mtx_t[G2])
	if tmpo.shape != (2,2): ## do not need to consider if no presences b/c those will be filtered 
		tmpo[0] = 0
		tmpo = tmpo[[0,1]]
	
	res = fisher_exact(table = tmpo, alternative=alt)
	return res[1]

## practice co-occurrence algorithm using the finches data from co-occur R package
finches_file = os.path.join(root, 'finches.csv')
finches = pd.read_csv(finches_file, dtype='a', header='infer', index_col=0)
mtx = finches

##
cooccur(mtx, 5, 0.15, 0.85)
