#!/usr/bin/env python

##############################################
# Created by Anders B. Dohlman               #
# Contact anders.dohlman@duke.edu            #
# Last updated 1-19-21                       #
# Publication 10.1016/j.chom.2020.12.001     #
##############################################

# USAGE: python compute_prevalence.py <PROJECT> <ASSAY>

import os
import sys
import pandas as pd
import numpy as np

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests, fdrcorrection


def compute_fisher(df, nobs):
	'''For each taxon, calculate fisher exact test statistics comparing the
	detection rates of two sample populations (eg. tissue and blood).'''

	for ix in df.index:
		count = np.array(df.loc[ix,['prev.pos','prev.neg']]) #+ 1
		stat, pval = fisher_exact(np.array([ 1 + count, 1 + nobs - count]),alternative='greater')
		df.loc[ix,'fisher.stat'] = np.log2(stat) if stat else 0.0
		df.loc[ix,'fisher.p'] = pval if pval else 1.0

	df['fisher.fdr'] = fdrcorrection(df['fisher.p'])[1]

	return df


def compute_prevalence_vs_blood(project, assay, domain, stat, min_reads=2, downsample=False):
	'''Calculate prevalence of each taxon in tissue (pos) vs. blood (neg).'''

	print("\nComparing {} tissue vs. blood".format(project))

	taxa = pd.read_csv('./taxonomy/taxa.all.txt',sep='\t',index_col=0)

	fname = './metadata/metadata.{}.all.txt'.format(project)
	if not os.path.exists(fname): return pd.DataFrame()
	meta = pd.read_csv(fname,sep='\t',index_col=0,low_memory=False)

	fname = './results/{}/{}/{}.{}.txt'.format(project, assay, domain, stat)
	if not os.path.exists(fname): return pd.DataFrame()
	data = pd.read_csv(fname,sep='\t',index_col=0)

	meta = meta.loc[data.columns]

	if 'BDN' not in meta['sample.Sample'].unique(): return pd.DataFrame()

	samples_neg = meta.index[meta['sample.Sample']=='BDN']
	samples_pos = meta.index[meta['sample.Sample'].isin(['PT','STN'])]
	print('{} tissue vs {} blood samples'.format(len(samples_pos),len(samples_neg)))

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


def compute_prevalence_vs_project(project, assay, domain, stat, sample_type, compare_projects=['GBM','LGG'], min_reads=2):
	'''Calculate prevalence of each taxon in tissue (pos) vs. control tissue (neg).'''

	print("\nComparing {} {} vs. {} {}".format(project, sample_type, compare_projects, sample_type))

	taxa = pd.read_csv('./taxonomy/taxa.all.txt',sep='\t',index_col=0)

	fname = './metadata/metadata.{}.all.txt'.format(project)
	if not os.path.exists(fname):
		print("ERROR: metadata file not found.")
		exit(1)

	meta_proj = pd.read_csv(fname,sep='\t',index_col=0,low_memory=False)

	ctrl_frames = []
	for project_ctrl in compare_projects:
		fname = './metadata/metadata.{}.all.txt'.format(project_ctrl)
		if not os.path.exists(fname): return pd.DataFrame()
		df = pd.read_csv(fname,sep='\t',index_col=0,low_memory=False)
		ctrl_frames.append(df)

	meta = pd.concat([meta_proj] + ctrl_frames)

	fname = './results/{}/{}/{}.{}.txt'.format(project, assay, domain, stat)
	data_proj = pd.read_csv(fname,sep='\t',index_col=0)

	ctrl_frames = []
	for project_ctrl in compare_projects:
		fname = './results/{}/{}/{}.{}.txt'.format(project_ctrl, assay, domain, stat)
		df = pd.read_csv(fname,sep='\t',index_col=0)
		ctrl_frames.append(df)

	data = pd.concat([data_proj] + ctrl_frames,axis=1)
	meta = meta.loc[data.columns]

	sample_names = {'blood':['BDN'],'tissue':['PT','STN']}[sample_type]
	meta = meta.loc[meta['sample.Sample'].isin(sample_names)]

	data = data[meta.index]
	meta = meta.loc[data.columns]

	samples_neg = meta.index[meta['case.acronym'].isin(project_ctrl.split('-'))].unique()
	samples_pos = meta.index[meta['case.acronym'].isin(project.split('-'))].unique()
	print('{} {} vs {} {} samples'.format(
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


if __name__ == "__main__":

	try:
		project, assay = sys.argv[1:]
	except:
		print("ERROR: Missing assay and/or project")
		exit(1)

	stat = 'unambiguous'
	domain = 'bacteria'

	# Select sterile tissue projects to compare with PROJECT
	use_projects_as_control = []

	print('Screening project {} {}'.format(project,assay))
	print('Domain:',domain)
	print('Statistic:',stat)

	compare_projects = ['GBM','LGG']

	project_ctrl = '-'.join(compare_projects) # These projects were concatenated before analysis
	project_ctrl_name = {'OV':'ovary','GBM-LGG':'brain'}[project_ctrl]

	print('\nUsing {} as control tissue'.format(project_ctrl_name))

	# Calculate project tissue vs. project blood prevalence
	df = compute_prevalence_vs_blood(project, assay, domain, stat)
	fname = './prevalence/prevalence.{}.{}.{}.{}.txt'.format(project,assay,domain,stat)
	df.to_csv(fname,sep='\t')

	# Calculate project tissue vs. control tissue prevalence
	df = compute_prevalence_vs_project(project, assay, domain, stat, sample_type='tissue', compare_projects=compare_projects, )
	fname = './prevalence/prevalence_{}.{}.{}.{}.{}.tissue.txt'.format(project_ctrl_name,project,assay,domain,stat)
	df.to_csv(fname,sep='\t')

	# Calculate project blood vs. control blood prevalence
	df = compute_prevalence_vs_project(project, assay, domain, stat, sample_type='blood', compare_projects=compare_projects)
	fname = './prevalence/prevalence_{}.{}.{}.{}.{}.blood.txt'.format(project_ctrl_name,project,assay,domain,stat)
	df.to_csv(fname,sep='\t')


