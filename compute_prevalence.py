#!/usr/bin/env python

# USAGE: python compute_prevalence.py -p $project -a $assay

##############################################
# Created by Anders B. Dohlman               #
# Contact anders.dohlman@duke.edu            #
# Last updated 11-29-21                      #
# Publication 10.1016/j.chom.2020.12.001     #
##############################################

import os
import argparse
import pandas as pd
import numpy as np

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection


def main():

	# Calculate project tissue vs. project blood prevalence
	df = compute_prevalence_vs_blood()
	fname = './prevalence/prevalence.{}.{}.{}.{}.txt'.format(PROJECT,ASSAY,DOMAIN,STAT)
	if len(df) > 0: df.to_csv(fname,sep='\t')
	else: print("Skipping tissue v.s. blood prevalence comparison")

	# Calculate project tissue vs. control tissue prevalence
	df = compute_prevalence_vs_project(sample_type='tissue', compare_projects=compare_projects)
	fname = './prevalence/prevalence_{}.{}.{}.{}.{}.tissue.txt'.format(TISSUE_COMP,PROJECT,ASSAY,DOMAIN,STAT)
	if len(df) > 0: df.to_csv(fname,sep='\t')
	else: print("Skipping {} tissue v.s. {} tissue prevalence comparison".format(PROJECT, TISSUE_COMP))

	# Calculate project blood vs. control blood prevalence
	df = compute_prevalence_vs_project(sample_type='blood', compare_projects=compare_projects)
	fname = './prevalence/prevalence_{}.{}.{}.{}.{}.blood.txt'.format(TISSUE_COMP,PROJECT,ASSAY,DOMAIN,STAT)
	if len(df) > 0: df.to_csv(fname,sep='\t')
	else: print("Skipping {} blood v.s. {} blood prevalence comparison".format(PROJECT, TISSUE_COMP))


def compute_fisher(df, nobs):
	'''For each taxon, calculate fisher exact test statistics comparing the
	detection rates of two sample populations (eg. tissue and blood).'''

	for ix in df.index:
		count = np.array(df.loc[ix,['prev.pos','prev.neg']]) #+ 1
		STAT, pval = fisher_exact(np.array([ 1 + count, 1 + nobs - count]),alternative='greater')
		df.loc[ix,'fisher.STAT'] = np.log2(STAT) if STAT else 0.0
		df.loc[ix,'fisher.p'] = pval if pval else 1.0

	df['fisher.fdr'] = fdrcorrection(df['fisher.p'])[1]

	return df


def compute_prevalence_vs_blood(min_reads=2, downsample=False):
	'''Calculate prevalence of each taxon in tissue (pos) vs. blood (neg).'''

	print("\nComparing {} tissue vs. {} blood...".format(PROJECT,PROJECT))

	taxa = pd.read_csv('./taxonomy/taxa.all.txt',sep='\t',index_col=0)

	fname = './metadata/metadata.{}.file.txt'.format(PROJECT)
	if not os.path.exists(fname):
		print("ERROR: no metadata file for {} present".format(PROJECT))
		return pd.DataFrame()
	meta = pd.read_csv(fname,sep='\t',index_col=0,low_memory=False)

	fname = './results/{}/{}/{}.{}.txt'.format(PROJECT, ASSAY, DOMAIN, STAT)
	if not os.path.exists(fname):
		print("ERROR: no results file present.")
		return pd.DataFrame()
	data = pd.read_csv(fname,sep='\t',index_col=0)

	meta = meta.loc[data.columns]

	if 'BDN' not in meta['sample.Sample'].unique():
		print("ERROR: no blood samples present")
		return pd.DataFrame()

	samples_neg = meta.index[meta['sample.Sample']=='BDN']
	samples_pos = meta.index[meta['sample.Sample'].isin(['PT','STN'])]
	print('{} tissue (pos) vs. {} blood (neg) samples'.format(len(samples_pos),len(samples_neg)))

	df = pd.DataFrame(index=data.index)
	data = data.fillna(0.0)
	df['prev.neg'] = np.array(data.loc[df.index,samples_neg].fillna(0.0) >= min_reads).sum(1)
	df['prev.pos'] = np.array(data.loc[df.index,samples_pos].fillna(0.0) >= min_reads).sum(1)
	df['prev.neg.pct'] = df['prev.neg']/len(samples_neg)
	df['prev.pos.pct'] = df['prev.pos']/len(samples_pos)

	nobs = np.array([len(samples_pos),len(samples_neg)])
	df = compute_fisher(df,nobs)
	df = df.join(taxa)

	return df


def compute_prevalence_vs_project(sample_type, compare_projects=['GBM','LGG'], min_reads=2):
	'''Calculate prevalence of each taxon in tissue (pos) vs. control tissue (neg).'''

	print("\nComparing {} {} vs. {} {}...".format(PROJECT, sample_type, TISSUE_COMP, sample_type))

	fname = './metadata/metadata.{}.file.txt'.format(PROJECT)
	if not os.path.exists(fname):
		print("ERROR: no metadata file for {} present".format(PROJECT))
		exit(1)

	meta_proj = pd.read_csv(fname,sep='\t',index_col=0,low_memory=False)

	ctrl_frames = []
	for project_ctrl in compare_projects:
		fname = './metadata/metadata.{}.file.txt'.format(project_ctrl)
		if not os.path.exists(fname):
			print("ERROR: Comparison {} metadata is not present".format(project_ctrl))
			return pd.DataFrame()
		df = pd.read_csv(fname,sep='\t',index_col=0,low_memory=False)
		ctrl_frames.append(df)

	meta = pd.concat([meta_proj] + ctrl_frames)

	fname = './results/{}/{}/{}.{}.txt'.format(PROJECT, ASSAY, DOMAIN, STAT)
	data_proj = pd.read_csv(fname,sep='\t',index_col=0)

	ctrl_frames = []
	for project_ctrl in compare_projects:
		fname = './results/{}/{}/{}.{}.txt'.format(project_ctrl, ASSAY, DOMAIN, STAT)
		if not os.path.exists(fname):
			print("ERROR: Comparison {} data is not present".format(project_ctrl))
			return pd.DataFrame()
		df = pd.read_csv(fname,sep='\t',index_col=0)
		ctrl_frames.append(df)

	data = pd.concat([data_proj] + ctrl_frames,axis=1)
	meta = meta.loc[data.columns]

	sample_names = {'blood':['BDN'],'tissue':['PT','STN']}[sample_type]

	if len( set(sample_names) & set(meta['sample.Sample'].unique()) ) == 0:
		print("ERROR: no {} samples present".format(sample_type))
		return pd.DataFrame()

	meta = meta.loc[meta['sample.Sample'].isin(sample_names)]

	data = data[meta.index]
	meta = meta.loc[data.columns]

	samples_neg = meta.index[meta['case.acronym'].isin(project_ctrl.split('-'))].unique()
	samples_pos = meta.index[meta['case.acronym'].isin(PROJECT.split('-'))].unique()
	print('{} {} (pos) vs. {} {} (neg) samples'.format(
		len(samples_pos),sample_type,len(samples_neg),sample_type))

	df = pd.DataFrame(index=data.index)

	data = data.fillna(0.0)
	df['prev.neg'] = np.array(data.loc[df.index,samples_neg].fillna(0.0) >= min_reads).sum(1)
	df['prev.pos'] = np.array(data.loc[df.index,samples_pos].fillna(0.0) >= min_reads).sum(1)
	df['prev.neg.pct'] = df['prev.neg']/len(samples_neg)
	df['prev.pos.pct'] = df['prev.pos']/len(samples_pos)

	nobs = np.array([len(samples_pos),len(samples_neg)])
	df = compute_fisher(df,nobs)
	df = df.join(taxa)

	return df


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
	parser.add_argument('-s','--statistic', nargs='?',default='unambiguous',
		help="""Statistic to acquire from pathseq output. Options include
		"unambiguous" (default), "score", "reads".""")
	parser.add_argument('-d','--domain', nargs='?',default='bacteria',
		help="""Domain or kingdom of interest. Options include:
		"bacteria" (default), "archaea", "fungi", "viruses".""")
	parser.add_argument('--tissue-comparison', nargs='?',default='brain',
		help="""Tissue type to use as comparision. Options: brain (default), ovary.
		Requires data for GBM/LGG and OV TCGA projects, respectively.""")

	return parser


if __name__ == "__main__":

	parser = create_parser()
	args = parser.parse_args()
	print(args)

	d = vars(args)
	PROJECT = d['project']
	ASSAY = d['assay']
	DOMAIN = d['domain']
	STAT = d['statistic']
	TISSUE_COMP = d['tissue_comparison']

	compare_projects = {'brain':['GBM','LGG'],'ovary':['OV']}[TISSUE_COMP]

	taxa = pd.read_csv('./taxonomy/taxa.all.txt',sep='\t',index_col=0)

	print('Computing prevalence for {} {}'.format(PROJECT, ASSAY))
	print('Domain:',DOMAIN)
	print('Statistic:',STAT)
	print('Tissue comparison:',TISSUE_COMP,compare_projects)

	main()
	print("Done.\n")


