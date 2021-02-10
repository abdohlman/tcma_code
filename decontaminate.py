#!/usr/bin/env python

# USAGE: python decontaminate.py $project $assay

##############################################
# Created by Anders B. Dohlman               #
# Contact anders.dohlman@duke.edu            #
# Last updated 2-5-21                        #
# Publication 10.1016/j.chom.2020.12.001     #
##############################################

import os
import sys
import pandas as pd
import numpy as np


def compute_mixture_fractions(data,ix,level_clf='species', min_count=5):
	'''Recursively estimates mixtures of tissue-resident and contaminant reads
	for all taxa. This is done using by calculating the relative fractions of
	unambiguously-aligned sequencing reads from tissue-resident or contaminant
	species within the taxa's clade.'''

	data = data.fillna(0.0)
	data[ data == 1 ] = 0

	fname = './taxonomy/tax_table.all.ids.txt'
	taxtable = pd.read_csv(fname, sep='\t',index_col=0)
	taxtable = taxtable.loc[taxtable[level_clf].isin(data.index)]

	taxa = pd.read_csv('./taxonomy/taxa.all.txt',sep='\t',index_col=0)

	df = pd.DataFrame(columns=data.columns)

	for level in taxtable.columns:

		if level == level_clf: break

		otus = taxtable[level].unique()

		for otu in otus:

			if otu not in data.index: continue
			children = taxtable.loc[taxtable[level] == otu,level_clf].unique()
			if len(children) == 0: continue
			
			C_ix = data.loc[np.intersect1d(ix,children)].sum() # sum clade in index
			C_clade = data.loc[data.index.isin(children)].sum() # sum clade

			C_clade[ C_clade < min_count ] = np.nan 

			pct_ix = C_ix.divide(C_clade)

			df.loc[otu] = pct_ix

	return df


def extract_subpopulation(pct, data, ix):
	'''Use estimated fractions (pct) and multiply with the dataset to get mixture
	from a given set of classified taxa (ix).'''

	data_sub = pd.concat([pct_raw.multiply(data.loc[pct.index]), data.loc[ix]])
	data_sub[ data_sub == 0 ] = np.nan
	data_sub = data_sub.dropna(0,'all')

	return data_sub


def nan_median(v,n=5):
	'''Returns the median of a series ignoring NaNs, unless there are
	fewer than n non-null values.'''

	if len(v) - sum(v.isnull()) < n: return np.nan
	else: return np.nanmedian(v)


def impute_mixtures(pct):
	'''The mixtures (tissue-resident/contamination) are difficult to estimate when
	very few reads from a given clade were found in a given sample. Assuming that
	mixtures are similar for samples processed at on the same plate/center, we
	first try to impute using the plate average (if there are sufficient wells),
	then impute using the center average. Since it primarily affects low-abundance
	taxa, imputation does not have a significant affect on the full dataset.'''

	fname = './metadata/metadata.{}.all.txt'.format(PROJECT)
	meta = pd.read_csv(fname,index_col=0,sep='\t',low_memory=False)
	meta = meta.loc[pct.columns]

	frames = []

	for sample_type in ['blood','tissue']:

		if sample_type == 'tissue':
			samples = meta.index[meta['sample.Sample'].isin(['PT','STN'])]
		else:
			samples = meta.index[meta['sample.Sample']=='BDN']

		if len(samples) == 0: continue

		plates = meta.loc[samples,'aliquot.plate_id'].value_counts()
		plates = plates.index[ plates > 5 ] # Only grab plates with < 5 instances

		centers = meta.loc[samples,'aliquot.center_name'].unique()

		df_pct = pct[samples].dropna(0,'all')

		df_real = (~df_pct.isnull())

		print(len(plates),'plates',len(centers),'centers')

		# Impute by plate
		for plate in plates:

			batch_samples = np.intersect1d(df_pct.columns, meta.index[meta['aliquot.plate_id']==plate])
			df_pct[batch_samples] = df_pct[batch_samples].apply(
				lambda x: x.fillna(nan_median(x)),axis=1)
			# print("plate {}: {:.2f}% sparsity".format(plate,100*df_pct.isnull().sum().sum()/len(df_pct.unstack())))

		# Impute by center
		for center in centers:
    
			batch_samples = np.intersect1d(df_pct.columns, meta.index[meta['aliquot.center_name']==center])
			df_pct[batch_samples] = df_pct[batch_samples].apply(
				lambda x: x.fillna(nan_median(x[df_real.loc[x.name,batch_samples]])),axis=1)
			# print("center {}: {:.2f}% sparsity".format(center,100*df_pct.isnull().sum().sum()/len(df_pct.unstack())))

		df_pct = df_pct.apply(lambda x: x.fillna(nan_median(x[df_real.loc[x.name]],n=1)),axis=1)
		# print("Final sparsity: {:.2f}%".format(100*df_pct.isnull().sum().sum()/len(df_pct.unstack())))
		frames.append(df_pct)
		print(df_pct.shape)

	pct_imputed = pd.concat(frames,axis=1)

	return pct_imputed


def predict_contaminants(PROJECT,assay,domain,stat, max_fdr=0.05, max_prev=0.2):
	'''Classifies taxa as contamination (ctm) or tissue-resident (tss) using prevalence.'''

	print("Classifying tissue-resident taxa from contaminants")

	fname = './prevalence/prevalence.{}.{}.{}.{}.txt'.format(PROJECT,assay,domain,stat)
	df = pd.read_csv(fname,sep='\t',index_col=0)
	df = df.loc[df['type']==level_clf]

	ix_tss = df.index[(df['fisher.fdr'] < max_fdr) & (df['prev.neg.pct'] < max_prev)]
	ix_ctm = df.index[~df.index.isin(ix_tss)]

	return ix_ctm, ix_tss



if __name__ == "__main__":

	# Use these parameters to classify taxa and estimate mixtures
	assay_clf = 'WGS'
	stat_clf = 'unambiguous'
	level_clf = 'species'

	# Apply classification to the following PROJECT, assay, domain, statistic
	try:
		PROJECT, ASSAY = sys.argv[1:]
	except:
		print("ERROR: Missing assay and/or PROJECT")
		exit(1)

	domain = 'bacteria'
	stat = 'unambiguous'

	print("Decontaminating the following dataset:")
	print(PROJECT, ASSAY, stat)
	print("Classification will be performed at the {} level using the following dataset:".format(level_clf))
	print(PROJECT, assay_clf, stat_clf)

	# Load taxonomy data
	taxa = pd.read_csv('./taxonomy/taxa.all.txt',sep='\t',index_col=0)

	# Load data to use for classification
	fname = './results/{}/{}/{}.{}.txt'.format(PROJECT,assay_clf,domain,stat_clf)
	data_clf = pd.read_csv(fname,sep='\t',index_col=0)

	# Load data to classify
	fname = './results/{}/{}/{}.{}.txt'.format(PROJECT,ASSAY,domain,stat)
	data = pd.read_csv(fname,sep='\t',index_col=0)

	# Classify tissue-resident vs. contaminant taxa
	ix_ctm_clf, ix_tss_clf = predict_contaminants(
			PROJECT,
			assay_clf,
			domain,
			stat_clf
	)

	# Intersect indices with data to be classified
	ix_all = data.index[data.index.isin(taxa.index[taxa['type']==level_clf])]
	ix_tss = ix_tss_clf
	ix_ctm = np.setdiff1d(ix_all, ix_tss_clf)

	print("Identified {} tissue-resident {}-level taxa".format(len(ix_tss),level_clf))

	# Save mixtures and assocated fractions for both subpopulations
	for sub_pop in ['decontam','contam']: 

		print("Extracting subplopulation {}".format(sub_pop))

		# Get species classified as tissue-resident (decontam) or contaminant (contam)
		ix = {'decontam':ix_tss, 'contam':ix_ctm}[sub_pop]

		# Use raw mixtures to compare tissue-resident vs. contamination fractions.
		pct_raw = compute_mixture_fractions(data, ix, level_clf=level_clf)
		fname = './mixtures/{}.{}.{}.{}.{}.raw.txt'.format(sub_pop,PROJECT,ASSAY,domain,stat_clf)
		pct_raw.to_csv(fname,sep='\t')

		# Use raw mixtures to generate the subpopulation.
		data_sub = extract_subpopulation(pct_raw, data, ix)
		fname = './results/{}/{}/{}.{}.{}.raw.txt'.format(PROJECT, ASSAY, domain, stat, sub_pop)
		data_sub.to_csv(fname,sep='\t')

		# Use imputed mixtures to compare tissue-resident vs. contamination fractions.
		pct_imputed = impute_mixtures(pct_raw)
		fname = './mixtures/{}.{}.{}.{}.{}.imputed.txt'.format(sub_pop,PROJECT,assay_clf,domain,stat_clf)
		pct_imputed.to_csv(fname,sep='\t')

		# Use imputed mixtures to generate the subpopulation.
		data_sub = extract_subpopulation(pct_imputed, data, ix)	
		fname = './results/{}/{}/{}.{}.{}.txt'.format(PROJECT, ASSAY, domain, stat, sub_pop)
		data_sub.to_csv(fname,sep='\t')


