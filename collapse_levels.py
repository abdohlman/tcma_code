#!/bin/env/bash python

# USAGE: python collapse_levels.py -p $project -a $assay

##############################################
# Created by Anders B. Dohlman               #
# Contact anders.dohlman@duke.edu            #
# Last updated 11-29-21                      #
# Publication 10.1016/j.chom.2020.12.001     #
##############################################


import os, argparse
import pandas as pd
import numpy as np

from skbio.stats.composition import clr


def main():

	# Normalize within sample types
	for sample_type in ['tissue','blood']:

		if sample_type == 'blood':
			ix_tissue = META.index[META['sample.Sample']=='BDN']
		elif sample_type == 'tissue':
			ix_tissue = META.index[META['sample.Sample'].isin(['PT','STN'])]

		samples = np.intersect1d(DATA.columns, ix_tissue)

		# data = DATA[samples].fillna(0)
		# meta = META.loc[samples]

		# Collapse to sample-, patient-, and file level
		for collapse_var in ['sample.bcr_sample_barcode','case.bcr_patient_barcode','']: #,'']:

			data = DATA[samples].fillna(0.0)
			meta = META.loc[samples]

			print("\nCollapsing {} to {}-level:".format(sample_type,collapse_var.split('.')[0] if collapse_var else 'file'))

			collapse_name = collapse_var.split('.')[0] if collapse_var else 'file'
			# stat_name = 'rpm.' + STAT if prenormalize_rpm else STAT
			reads_tag = 'reads' if no_rpm else 'rpm'

			print("Normalizing to read average...")
			data_reads = collapse_reads(data, meta, collapse_var=collapse_var, use_rpm=no_rpm is False)
			fname = '{}/{}.{}.{}.{}.{}.txt'.format(OUTPUT_DIR, domain, STAT, sample_type, collapse_name, reads_tag)
			print("Writing to",fname)
			data_reads.to_csv(fname,sep='\t')
			print('\t',data_reads.shape)

			print(data_reads.head())

			print("Normalizing to relative abundance...")
			data_relabund = collapse_relabund(data, meta, collapse_var=collapse_var, use_rpm=no_rpm is False)
			fname = '{}/{}.{}.{}.{}.{}.relabund.txt'.format(OUTPUT_DIR, domain, STAT, sample_type, collapse_name, reads_tag)
			print("Writing to",fname)
			data_relabund.to_csv(fname,sep='\t')
			print('\t',data_relabund.shape)

			print(data_relabund.head())

			print("Normalizing to centered-log ratio (CLR)...")
			data_clr = collapse_clr(data, meta, collapse_var=collapse_var, use_rpm=no_rpm is False, min_prev=min_clr_prev)
			fname = '{}/{}.{}.{}.{}.{}.clr.txt'.format(OUTPUT_DIR, domain, STAT, sample_type, collapse_name, reads_tag)
			print("Writing to",fname)
			data_clr.to_csv(fname,sep='\t')
			print('\t{}'.format(data_clr.shape))

			print(data_clr.head())


def normalize_rpm(df, meta, collapse_var=''):
	'''	Normalizes by total reads. If collapsing, sum counts and 
	total reads across samples.
		V = (c_1 + ... + c_N)/(t_1 + ... + t_N)
		c_i = counts in sample i
		t_i = total reads in sample i
		N = number of items to collapse
	'''
	df = df.transpose()
	df = df.join(READS)
	if collapse_var:
		df = df.join(meta[collapse_var]).groupby(collapse_var).sum()
	df = 1e6*df.divide(df.PRIMARY_READS, axis=0).drop(columns='PRIMARY_READS')
	df = df.transpose()

	# print("Normalized rpm")
	# print(df.head())

	return df


def normalize_relabund(df):
	'''Normalizes samples to relative abundance with each taxonomic rank.'''

	ranked_otus = []
	for taxon in taxtable.columns:
		otus = np.intersect1d(taxtable[taxon].unique(), df.index)
		if len(otus) == 0: continue
		df.loc[otus] = df.loc[otus].divide(df.loc[otus].sum())
		ranked_otus += list(otus)

	df = df.loc[ranked_otus]

	return df


def normalize_clr(df):
	'''Normalizes samples to centered log-ratio (CLR) with each taxonomic rank.'''

	ranked_otus = []
	for taxon in taxtable.columns: 
		otus = np.intersect1d(taxtable[taxon].unique(), df.index)
		if len(otus) == 0: continue
		df.loc[otus] = df.loc[otus].apply(clr,1,result_type='broadcast') #broadcast=True)#result_type='broadcast')
		ranked_otus += list(otus)

	df = df.loc[ranked_otus]

	return df


def collapse_reads(df, meta, collapse_var='', use_rpm=True):
	''''Average read counts using geometric mean'''

	if use_rpm:
		df = normalize_rpm(df, meta, collapse_var)
	elif collapse_var:
		df = df.transpose()
		df = np.log10(df.fillna(0.0) + 1)
		df = df.join(meta[collapse_var])
		df = df.groupby(collapse_var).mean()
		df = 10**df - 1
		df = df.transpose() 

	df = df.fillna(0.0)

	return df


def collapse_relabund(df, meta, collapse_var='',use_rpm=True):
	'''Normalizes dataset to sample-wise relative abundance, such that
	the sum of all taxa of a given rank is equal to 1 for all samples.
		V = c_1/Sum_1 + ... + c_N/Sum_N
		c_i := counts in sample i
		Sum_i := sum of reads in sample i for a given rank
		N := number of items to collapse
	'''

	if use_rpm:
		df = normalize_rpm(df, meta, collapse_var)

	elif collapse_var:

		df = normalize_relabund(df)
		df = df.transpose()
		df = df.join(meta[collapse_var])
		df = df.groupby(collapse_var).sum()
		df = df.transpose()

	df = normalize_relabund(df)

	df = df.fillna(0.0)

	return df


def collapse_clr(df, meta, collapse_var='',use_rpm=True, min_prev=0.25):
	'''Collapse by relative abundance then calculate CLR.'''

	keep_otus = [ otu for otu in df.index if sum(df.loc[otu] > 0) > min_prev*len(df.columns) ]

	df = df.loc[keep_otus]
	df = df.fillna(0.0) + 0.01 # pseudocount

	if use_rpm:
		df = normalize_rpm(df, meta, collapse_var)
		df = normalize_relabund(df)

	elif collapse_var:
		df = collapse_relabund(df, meta, collapse_var, use_rpm)

	df = normalize_clr(df)

	return df


def get_total_reads():
	'''Get read statistics.'''

	fpath = '{}/read_statistics.txt'.format(RESULTS_DIR)
	df_reads = pd.read_table(fpath, index_col=0)
	read_counts = df_reads['PRIMARY_READS']

	assert sum(read_counts.isna()) == 0
	assert read_counts.min() > 0

	return read_counts


def create_parser():

	parser = argparse.ArgumentParser(
		description="""Parses the PathSeq output files from ./pathseq_out and
		collates unambiguously-aligned reads for each taxa into an n x m table,
		with n taxa (NCBI taxonomy ID) and m sequencing runs (UUIDs). This
		table is saved to ./results/$project/$assay.""")

	parser.add_argument('-p','--project', nargs='?', required=True,
		help="""TCGA sequencing project (eg. COAD)""")

	parser.add_argument('-a','--assay', nargs='?', required=True,
		help="""TCGA experimental strategy (eg. WGS)""")

	parser.add_argument('-s','--statistic', nargs='?',default='unambiguous.decontam',
		help="""Statistic to acquire from pathseq output. Options include
		"unambiguous.decontam" (default), "score", "reads".""")

	parser.add_argument('-d','--domain', nargs='?',default='bacteria',
		help="""Domain or kingdom of interest. Options include:
		"bacteria" (default), "archaea", "fungi", "viruses".""")

	parser.add_argument('--no-rpm',action='store_true',default=False,
		help="""Do not normalize to reads-per-million prior to collapsing.""")

	parser.add_argument('--min-clr-prevalence',default=0.25,
		help="""Remove OTUs with lesser than this fraction of prevalence prior
		to calculating CLR. The CLR transform performs best with lower sparsity.
		Default = 0.25.""")

	return parser


if __name__ == "__main__":

	domain = 'bacteria'
	prenormalize_rpm = True # This was set to false for the manuscript

	parser = create_parser()
	args = parser.parse_args()
	print(args)

	d = vars(args)
	PROJECT = d['project']
	ASSAY = d['assay']
	DOMAIN = d['domain']
	STAT = d['statistic']
	no_rpm = d['no_rpm']
	min_clr_prev = float(d['min_clr_prevalence']) # Use 0.05 for pan-cancer, 0.25 otherwise

	OUTPUT_DIR = './collapsed/{}/{}'.format(PROJECT, ASSAY)
	RESULTS_DIR = './results/{}/{}'.format(PROJECT, ASSAY)

	if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)

	taxa = pd.read_csv('./taxonomy/taxa.all.txt',index_col=0,sep='\t')

	fname = './taxonomy/tables/tax_table.all.ids.txt'
	taxtable = pd.read_csv(fname, sep='\t',index_col=0)

	fname = '{}/{}.{}.txt'.format(RESULTS_DIR, domain, STAT)
	DATA = pd.read_csv(fname,sep='\t',index_col=0)

	fname = './metadata/metadata.{}.file.txt'.format(PROJECT)
	META = pd.read_csv(fname,sep='\t',index_col=0, low_memory=False)

	READS = get_total_reads()

	main()
	print("Done.\n")
