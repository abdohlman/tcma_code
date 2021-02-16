#!/bin/bash python

# USAGE: python collapse_levels.py $project $assay $statistic

##############################################
# Created by Anders B. Dohlman               #
# Contact anders.dohlman@duke.edu            #
# Last updated 2-16-21                       #
# Publication 10.1016/j.chom.2020.12.001     #
##############################################


import os
import sys
import pandas as pd
import numpy as np

from skbio.stats.composition import clr


def normalize_reads(df, meta, collapse_var=''):
	''''Normalize reads using geometric mean'''
	if collapse_var:
		df = df.transpose()
		df = np.log10(df.fillna(0.0) + 1)
		df = df.join(meta[collapse_var])
		df = df.groupby(collapse_var).mean()
		df = 10**df - 1
		df = df.transpose() 

	df = df.fillna(0.0)
	return df



def normalize_relabund(df, meta, collapse_var=''):
	'''Normalizes dataset to sample-wise relative abundance, such that
	the sum of all taxa of a given rank is equal to 1 for all samples.'''

	df = df.fillna(0.0)

	for taxon in taxtable.columns:

		otus = np.intersect1d(taxtable[taxon].unique(), df.index)
		if len(otus) == 0: continue

		df.loc[otus] = df.loc[otus].divide(df.loc[otus].sum())

	if collapse_var:
		df = collapse_relabund(df, meta, collapse_var=collapse_var)

	df = df.fillna(0.0)

	return df


def normalize_clr(df, meta, collapse_var='', min_prev=0.25):
	'''Applies the centered log-ratio transform to the dataset. This
	maps relative abundances to a continuous, normally-distributed
	distribution.'''

	keep_otus = [ otu for otu in df.index if sum(df.loc[otu] > 0) > min_prev*len(df.columns) ]

	df = df.loc[keep_otus]
	df = df.fillna(0.0) + 0.1 # Ensure no zeros

	if collapse_var:
		df = collapse_relabund(df, meta, collapse_var=collapse_var)
		df = normalize_relabund(df, meta)

	for taxon in taxtable.columns: 
		otus = np.intersect1d(taxtable[taxon].unique(), df.index)
		if len(otus) == 0: continue

		df.loc[otus] = df.loc[otus].apply(clr,1,result_type='broadcast') #broadcast=True)#result_type='broadcast')

	return df


def collapse_relabund(df,meta,collapse_var='sample.bcr_sample_barcode'):
	''''Normalize to relative abundance after summing by collapse_var.'''

	df = df.transpose().fillna(0.0)
	df = df.join(meta[collapse_var])
	df = df.groupby(collapse_var).sum()
	df = df.transpose()
	df = normalize_relabund(df, meta)
	return df


def normalize_rpm(df, meta, flagstats_var='flagstats.total'):
	'''Normalize to reads per million, calculated according to total sequencing reads
	as measured using samtools flagstats.'''

	read_counts = meta[flagstats_var]

	for center in meta['aliquot.center_name'].unique():
		ix_center = meta.index[meta['aliquot.center_name']==center]
		ix_nan = read_counts.index[read_counts.isna()]
		ix_impute = np.intersect1d(ix_center,ix_nan)   
		read_counts.loc[ix_impute] = read_counts.loc[ix_center].median()

	df = 1e6*df.divide(read_counts,axis=1)

	return df


if __name__ == "__main__":

	domain = 'bacteria'
	prenormalize_rpm = True # This was set to false for the manuscript

	try:
		PROJECT, ASSAY, STAT = sys.argv[1:]
	except:
		print("ERROR: Missing ASSAY and/or PROJECT")
		exit(1)

	taxa = pd.read_csv('./taxonomy/taxa.all.txt',index_col=0,sep='\t')

	fname = './taxonomy/tax_table.all.ids.txt'
	taxtable = pd.read_csv(fname, sep='\t',index_col=0)

	fname = './results/{}/{}/{}.{}.txt'.format(PROJECT, ASSAY, domain, STAT)
	DATA = pd.read_csv(fname,sep='\t',index_col=0)

	fname = './metadata/metadata.{}.all.txt'.format(PROJECT)
	META = pd.read_csv(fname,sep='\t',index_col=0)

	# Normalize within sample types
	for tissue in ['solid','blood']:

		if tissue == 'blood':
			ix_tissue = META.index[META['sample.Sample']=='BDN']
		elif tissue == 'solid':
			ix_tissue = META.index[META['sample.Sample'].isin(['PT','STN'])]

		samples = np.intersect1d(DATA.columns, ix_tissue)

		data = DATA[samples]
		meta = META.loc[samples]

		if prenormalize_rpm:
			print("Pre-normalizing to RPM...")
			data = normalize_rpm(data,meta)

		# Collapse to sample-, patient-, and file level
		for collapse_var in ['sample.bcr_sample_barcode','case.bcr_patient_barcode','']:

			print("Collapsing to {}-level:".format(collapse_var.split('.')[0] if collapse_var else 'file'))

			collapse_name = collapse_var.split('.')[0] if collapse_var else 'file'
			stat_name = 'rpm.' + STAT if calculate_rpm else STAT

			print("Normalizing to read average...")
			data_reads = normalize_reads(data, meta, collapse_var=collapse_var)
			fname = './results/{}/{}/{}.{}.{}.{}.reads.txt'.format(PROJECT, ASSAY, domain, stat_name, tissue, collapse_name)
			data_reads.to_csv(fname,sep='\t')

			print("Normalizing to relative abundance...")
			data_relabund = normalize_relabund(data, meta, collapse_var=collapse_var)
			fname = './results/{}/{}/{}.{}.{}.{}.relabund.txt'.format(PROJECT, ASSAY, domain, stat_name, tissue, collapse_name)
			data_relabund.to_csv(fname,sep='\t')

			print("Normalizing to centered-log ratio (CLR)...")
			data_clr = normalize_clr(data, meta, collapse_var=collapse_var)
			fname = './results/{}/{}/{}.{}.{}.{}.clr.txt'.format(PROJECT, ASSAY, domain, stat_name, tissue, collapse_name)
			data_clr.to_csv(fname,sep='\t')

