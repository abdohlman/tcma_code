#!/usr/bin/env python

# USAGE: python collate_pathseq_output.py $project $assay $min_clipped_read_length

##############################################
# Created by Anders B. Dohlman               #
# Contact anders.dohlman@duke.edu            #
# Last updated 2-10-21                       #
# Publication 10.1016/j.chom.2020.12.001     #
##############################################

import os, sys, argparse
import pandas as pd

TAXONOMY = ['superkingdom',
			'kingdom',#'subkingdom',
			'phylum',#'subphylum',
			'class',#'subclass',
			'order',#'suborder',
			'family',#'subfamily',
			'genus',
			'species']

def parse_results(resultfiles):
	'''Iterates through pathseq output files and creates a matrix with
	dimensions (m taxa) x (n sequencing runs).'''

	missing_ct = 0 # Count missing files

	frames = []

	for fname in resultfiles:

		if not fname in os.listdir(OUT_DIR):
			missing_ct += 1
			continue

		filepath = os.path.join(OUT_DIR, fname)
		df = pd.read_table(filepath, index_col=0)
		frames.append(df[STAT].rename(fname2id[fname]))

	if missing_ct > 0:
		print('WARNING: missing pathseq output for {} of {} files'.format(
			missing_ct,len(resultfiles)))

	data = pd.concat(frames,axis=1)

	# Filter by domain
	ix_domain = taxa.index[taxa['kingdom'].isin([DOMAIN.title(),'root'])]
	data_domain = data.loc[data.index.isin(ix_domain)]

	print("Resulting data table:",data_domain.shape)
	fpath = os.path.join(RESULT_DIR,'{}.{}.txt'.format(DOMAIN.lower(),STAT))
	data_domain.to_csv(fpath,sep='\t')


def make_dicts(meta):
	'''Dictonaries for converting between file name and file uuid.'''
	ids = meta.index.values
	fnames = meta['filename'].values
	# fbases  = [ f.replace('.bam','') for f in fnames ]
	fnames = [ f.replace('.bam','.{}.scores.txt'.format(minClippedReadLength)) for f in fnames ]
	fname2id = dict(zip(fnames, ids))
	id2fname = dict(zip(ids, fnames))

	return fname2id, id2fname


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
	parser.add_argument('-m','--min_clipped_read_length', nargs='?', required=True,
		help="""Minimum clipped read length given to pathseq.sh""")
	parser.add_argument('-s','--statistic', nargs='?',default='unambiguous',
		help="""Statistic to acquire from pathseq output. Options include
		"unambiguous" (default), "score", "reads".""")
	parser.add_argument('-d','--domain', nargs='?',default='bacteria',
		help="""Domain or kingdom of interest. Options include:
		"bacteria" (default), "archaea", "fungi", "viruses".""")

	return parser


if __name__ == "__main__":

	parser = create_parser()
	args = parser.parse_args()
	print(args)

	d = vars(args)
	PROJECT = d['project']
	ASSAY = d['assay']
	minClippedReadLength = d['min_clipped_read_length']
	DOMAIN = d['domain']
	STAT = d['statistic']

	OUT_DIR = './pathseq_out/{}/{}'.format(PROJECT, ASSAY)
	RESULT_DIR = './results/{}/{}'.format(PROJECT, ASSAY)

	if not os.path.exists(RESULT_DIR): os.makedirs(RESULT_DIR)

	METAFILE = './manifests/gdc_manifest.{}.{}.txt'.format(ASSAY, PROJECT)
	TAXAFILE = './taxonomy/taxa.all.txt'

	meta = pd.read_table(METAFILE,index_col=0)
	taxa  = pd.read_table(TAXAFILE,index_col=0)
	fname2id, id2fbase = make_dicts(meta)

	resultfiles = list(fname2id.keys())

	print("Collecting {} data for TCGA project {}...".format(ASSAY, PROJECT))
	print("Domain: {}".format(DOMAIN))
	print("Statistic: {}".format(STAT))
	parse_results(resultfiles)
	print("Done.")


