#!/usr/bin/env python

# USAGE: python decontaminate.py -p $project -a $assay

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

TAXONOMY = [
    'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'
]

def main():

	print("\nDecontaminating the following dataset:")
	print("Project = {}".format(PROJECT))
	print("Statistic = {}".format(STAT))
	print("Assay = {}".format(ASSAY))

	print("\nClassification to be performed with the following parameters:")
	print("Rank = {}".format(rank_clf))
	print("Project = {}".format(PROJECT))
	print("Statistic = {}".format(stat_clf))
	print("Assay = {}".format(assay_clf))

	# Ensure no missing taxa
	assert len(data_obs.index) == len(np.intersect1d(data_obs.index, taxa.index))

	if custom_list:
		print('Using custom list {}'.format(custom_list))
		with open(custom_list) as f:
			ix_tss_clf = [ int(s) for s in f.read().splitlines() ]
		f.close()
	else:
		ix_ctm_clf, ix_tss_clf = predict_contaminants(
			PROJECT,
			assay_clf,
			DOMAIN,
			stat_clf
		)

	ix_all = data_obs.index[data_obs.index.isin(taxa.index[taxa['type'] == rank_clf])]
	ix_tss = np.intersect1d(ix_all, ix_tss_clf)
	# ix_ctm = np.setdiff1d(ix_all, ix_tss_clf)
	ix_ctm = np.setdiff1d(ix_all, ix_tss)

	print("Classified {} of {} {}-level taxa as tissue-resident".format(
		len(ix_tss), len(ix_all), rank_clf))

	data_tss, data_ctm, pct_tss, pct_ctm = compute_mixture_fractions(
        data_obs, ix_tss, ix_ctm, rank_clf=rank_clf
    )

	results_path = './results/{}/{}/'.format(PROJECT, ASSAY)
	if not os.path.exists(results_path):
	    os.mkdir(results_path)

	fname = './results/{}/{}/{}.{}.{}.txt'.format(
	        PROJECT, ASSAY, DOMAIN, STAT, 'decontam')
	data_tss.to_csv(fname, sep='\t')

	fname = './results/{}/{}/{}.{}.{}.txt'.format(
	        PROJECT, ASSAY, DOMAIN, STAT, 'contam')
	data_ctm.to_csv(fname, sep='\t')

	fname = './mixtures/{}.{}.{}.{}.{}.txt'.format(
	        'decontam', PROJECT, ASSAY, DOMAIN, stat_clf)
	pct_tss.to_csv(fname, sep='\t')

	fname = './mixtures/{}.{}.{}.{}.{}.txt'.format(
	        'contam', PROJECT, ASSAY, DOMAIN, stat_clf)
	pct_ctm.to_csv(fname, sep='\t')


def get_children(clade):
	'''Returns a list of the nearest children for a given clade from
	the taxonomy table.'''

	children = clade.apply(lambda x: x.dropna()[-1], axis=1)
	children.name = 'nearest_child'

	children = children.unique()

	return children


def compute_mixture_fractions(data_obs, ix_tss, ix_ctm, rank_clf='species', min_count=5):
    '''Recursively estimates mixtures of tissue-resident and contaminant reads
    for all taxa. This is done using by calculating the relative fractions of
    unambiguously-aligned sequencing reads from tissue-resident or contaminant
    species within the taxa's clade.'''

    # Observed data   
    data_obs = data_obs.fillna(0.0)
    data_obs[data_obs == 1] = 0 # remove singletons

    # Instantiate mixture components & mixture fractions
    data_tss = pd.DataFrame(index=data_obs.index, columns=data_obs.columns)
    data_ctm = pd.DataFrame(index=data_obs.index, columns=data_obs.columns)
    pct_tss = pd.DataFrame(index=data_obs.index, columns=data_obs.columns)
    pct_ctm = pd.DataFrame(index=data_obs.index, columns=data_obs.columns)

    # Taxonomy ranks
    taxonomy_order = taxtable.columns[::-1]

    n_rank_clf = list(taxonomy_order).index(rank_clf)
    n_rank_kingdom = list(taxonomy_order).index('superkingdom')
    taxonomy_order = list(taxonomy_order[n_rank_clf:n_rank_kingdom + 1])

    otus_clf = taxtable[rank_clf].dropna().unique()

    assert len(otus_clf) == len(ix_tss) + len(ix_ctm)

    pct_tss.loc[ix_tss] = 1
    pct_tss.loc[ix_ctm] = 0

    pct_ctm.loc[ix_ctm] = 1
    pct_ctm.loc[ix_tss] = 0

    data_tss.loc[otus_clf] = pct_tss.loc[otus_clf]*data_obs.loc[otus_clf]
    data_ctm.loc[otus_clf] = pct_ctm.loc[otus_clf]*data_obs.loc[otus_clf]

    for n, rank in enumerate(taxonomy_order):

        if rank == rank_clf: #species 
            continue
            
        # print('[[ {} ]]'.format(rank.upper()))
        print('Calculating mixture at {}-level...'.format(rank))


        # List OTUs at this rank in data_obs
        otus = np.intersect1d(taxtable[rank].dropna().unique(), data_obs.index)

        # Determine rank below
        ranks_below = taxonomy_order[0:taxonomy_order.index(rank)]

        for otu in otus:

        	# Get nearest available children within clade
            clade = taxtable.loc[taxtable[rank] == otu, ranks_below]
            children = get_children(clade)
            children = np.intersect1d(data_obs.index, children)

            # Calculate sum of reads among nearest children
            children_sum = data_tss.loc[children].sum() + data_ctm.loc[children].sum()
            children_sum = children_sum.replace(0, np.nan) # avoid divide by zero error

            # Compute fractions
            pct_tss.loc[otu] = data_tss.loc[children].sum().divide(children_sum).fillna(0)
            pct_ctm.loc[otu] = data_ctm.loc[children].sum().divide(children_sum).fillna(0)

            # Compute mixtures
            data_tss.loc[otu] = (pct_tss.loc[otu]*data_obs.loc[otu]).fillna(0)
            data_ctm.loc[otu] = (pct_ctm.loc[otu]*data_obs.loc[otu]).fillna(0)

    return data_tss, data_ctm, pct_tss, pct_ctm


def nan_median(v, n=5):
	'''Returns the median of a series ignoring NaNs, unless there are
	fewer than n non-null values.'''

	if len(v) - sum(v.isnull()) < n: return np.nan
	else: return np.nanmedian(v)


def predict_contaminants(PROJECT,assay,domain,stat, max_fdr=0.05, max_prev=0.2):
	'''Classifies taxa as contamination (ctm) or tissue-resident (tss) using prevalence.'''

	print("\nClassifying tissue-resident taxa from contaminants...")

	fname = './prevalence/prevalence.{}.{}.{}.{}.txt'.format(PROJECT,assay,domain,stat)
	df = pd.read_csv(fname,sep='\t',index_col=0)
	df = df.loc[df['type']==rank_clf]

	ix_tss = df.index[(df['fisher.fdr'] < max_fdr) & (df['prev.neg.pct'] < max_prev)]
	ix_ctm = df.index[~df.index.isin(ix_tss)]

	return ix_ctm, ix_tss


def create_parser():

    parser = argparse.ArgumentParser(
        description="""Parses the PathSeq output files from ./pathseq_out and
        collates unambiguously-aligned reads for each taxa into an n x m table,
        with n taxa (NCBI taxonomy ID) and m sequencing runs (UUIDs). This
        table is saved to ./results/$project/$assay.""")

    parser.add_argument(
        '-p', '--project', nargs='?', required=True,
        help="""TCGA sequencing project (eg. COAD)""")
    parser.add_argument(
        '-a', '--assay', nargs='?', required=True,
        help="""TCGA experimental strategy (eg. WGS)""")
    parser.add_argument(
        '-s', '--statistic', nargs='?', default='unambiguous',
        help="""Statistic to acquire from pathseq output. Options include
        "unambiguous" (default), "score", "reads".""")
    parser.add_argument(
        '-d', '--domain', nargs='?', default='bacteria',
        help="""Domain or kingdom to classify. Options include:
        "bacteria" (default), "archaea", "fungi", "viruses".""")
    parser.add_argument(
        '--q-value', nargs='?', default=0.05,
        help="""Maximum FDR q-value at which to classify taxa
        as tissue-resident. Default: 0.05".""")
    parser.add_argument(
        '--max-prevalence-blood', nargs='?', default=0.2,
        help="""Maximum blood prevalence at which to classify
        taxa as tissue-resident. Default: 0.2 (20%%)""")
    parser.add_argument(
        '--rank-for-classification', nargs='?', default='species',
        help="""Taxonomic rank to be used for classification. Abundance
        for higher-order taxa will be adjusted recursively. No features
        below this rank will be preserved. Default: species.""")
    parser.add_argument(
        '--stat-for-classification', nargs='?', default='unambiguous',
        help="""Pathseq output statistic to be be used for classification.
        Default: unambiguous.""")
    parser.add_argument(
        '--assay-for-classification', nargs='?', default='WGS',
        help="""Assay to be be used for classification. Default: WGS.""")
    parser.add_argument(
        '--custom-tissue-resident-list', nargs='?', default='',
        help="""Use a custom list of tissue-resident taxa. Taxa not in this
        list will be classified as contamination. Provide a file path to a
        text file containing the list of taxa separated by newlines.""")

    return parser


if __name__ == "__main__":

    parser = create_parser()
    args = parser.parse_args()
    print(args)

    # Load arguments
    d = vars(args)
    PROJECT = d['project']
    ASSAY = d['assay']
    DOMAIN = d['domain']
    STAT = d['statistic']
    QVALUE = d['q_value']
    MAXPREV = d['max_prevalence_blood']
    rank_clf = d['rank_for_classification']
    stat_clf = d['stat_for_classification']
    assay_clf = d['assay_for_classification']
    custom_list = d['custom_tissue_resident_list']

    # Load taxonomy index
    fname = './taxonomy/taxa.all.txt'
    taxa = pd.read_table(fname, index_col=0)

    # Load data to use for classification
    fname = './results/{}/{}/{}.{}.txt'.format(
        PROJECT, assay_clf, DOMAIN, stat_clf)
    data_clf = pd.read_table(fname, index_col=0)

    # Load observed data to be decomposed
    fname = './results/{}/{}/{}.{}.txt'.format(PROJECT, ASSAY, DOMAIN, STAT)
    data_obs = pd.read_table(fname, index_col=0)

    # Load taxonomy table
    fname = './taxonomy/tables/tax_table.all.ids.txt' #.format(DOMAIN)
    taxtable = pd.read_table(fname, index_col=0)
    taxtable = taxtable.loc[taxtable[rank_clf].isin(data_obs.index)]

    main()
    print("Done.\n")

