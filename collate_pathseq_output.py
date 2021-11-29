#!/usr/bin/env python

# USAGE: python collate_pathseq_output.py $project $assay $min_clipped_read_length

##############################################
# Created by Anders B. Dohlman               #
# Contact anders.dohlman@duke.edu            #
# Last updated 11-29-21                      #
# Publication 10.1016/j.chom.2020.12.001     #
##############################################

import os
import argparse
import pandas as pd

def main():

	parse_results()
	parse_read_statistics()


def make_dicts(meta):
	'''Creates dictonaries for converting between file name and file uuid.'''
	ids = meta.index.values
	fnames = meta['filename'].values
	fbases = [ f.replace('.bam', '') for f in fnames ]
	fbases2id = dict(zip(fbases, ids))
	id2fbases = dict(zip(ids, fbases))

	return fbases2id, id2fbases


def parse_results():
	'''Iterates through pathseq output files and creates a matrix with
	dimensions (m taxa) x (n sequencing runs).'''

	missing_ct = 0 # Count missing files

	frames = []

	fbases = fbase2id.keys()

	for fbase in fbases:

		fname = '{}.{}.scores.txt'.format(fbase, MCRL)

		if not fname in os.listdir(OUT_DIR):
			missing_ct += 1
			continue

		filepath = os.path.join(OUT_DIR, fname)
		df = pd.read_table(filepath, index_col=0)
		frames.append(df[STAT].rename(fbase2id[fbase]))

	if missing_ct > 0:
		print('WARNING: missing pathseq output for {} of {} files'.format(
			missing_ct, len(fbases)))

	data = pd.concat(frames, axis=1)

	# Filter by domain
	ix_domain = taxa.index[taxa['kingdom'].isin([DOMAIN.title(),'root'])] # include no rank?
	data_domain = data.loc[data.index.isin(ix_domain)]

	print("OTUs table:",data_domain.shape)
	fpath = os.path.join(RESULT_DIR,'{}.{}.txt'.format(DOMAIN.lower(),STAT))
	data_domain.to_csv(fpath,sep='\t')
	

def parse_read_statistics():
	'''Iterates through pathseq output filter metrics and scores metrics
	and saves them to read_statistics.txt.'''

	missing_ct_fm = 0 # Count missing filter metrics
	missing_ct_sm = 0 # Count missing scores metrics

	frames = []

	fbases = fbase2id.keys()

	for fbase in fbases:

		fname = '{}.{}.filter_metrics.txt'.format(fbase, MCRL)
		filepath = os.path.join(OUT_DIR, fname)
		if not os.path.exists(filepath):
			missing_ct_fm += 1
			continue

		df_filt = pd.read_table(filepath, skiprows=6)
		df_filt.index = [fbase2id[fbase]]

		fname = '{}.{}.scores_metrics.txt'.format(fbase, MCRL)
		filepath = os.path.join(OUT_DIR, fname)
		if not os.path.exists(filepath):
			missing_ct_sm += 1
			continue

		df_score = pd.read_table(filepath, skiprows=6)
		df_score.index = [fbase2id[fbase]]

		df = pd.concat([df_filt, df_score], axis=1)

		frames.append(df)

	if missing_ct_fm > 0:
		print('WARNING: missing filter metrics for {} of {} files'.format(
			missing_ct_fm, len(fbases)))

	if missing_ct_sm > 0:
		print('WARNING: missing scores metrics for {} of {} files'.format(
			missing_ct_sm, len(fbases)))

	read_statistics = pd.concat(frames)
	print("Read statistics table:",read_statistics.shape)

	fname = 'read_statistics.txt'.format(PROJECT, ASSAY, MCRL)
	fpath = os.path.join(RESULT_DIR, fname)
	read_statistics.to_csv(fpath,sep='\t')


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
	parser.add_argument('-m','--min-clipped-read-length', nargs='?', required=True,
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
	MCRL = d['min_clipped_read_length']
	DOMAIN = d['domain']
	STAT = d['statistic']

	OUT_DIR = './pathseq_out/{}/{}'.format(PROJECT, ASSAY)
	RESULT_DIR = './results/{}/{}'.format(PROJECT, ASSAY)

	if not os.path.exists(RESULT_DIR): os.makedirs(RESULT_DIR)

	METAFILE = './manifests/gdc_manifest.{}.{}.txt'.format(ASSAY, PROJECT)
	TAXAFILE = './taxonomy/taxa.all.txt'

	meta = pd.read_table(METAFILE,index_col=0)
	taxa  = pd.read_table(TAXAFILE,index_col=0)

	fbase2id, id2fbase = make_dicts(meta)

	print("Collecting {} data for project {}...".format(ASSAY, PROJECT))
	print("Domain: {}".format(DOMAIN))
	print("Statistic: {}".format(STAT))

	main()
	print("Done.\n")


