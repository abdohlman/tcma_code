#!/usr/bin/env python

# USAGE: python collate_pathseq_output.py $project $assay $min_clipped_read_length

##############################################
# Created by Anders B. Dohlman               #
# Contact anders.dohlman@duke.edu            #
# Last updated 2-10-21                       #
# Publication 10.1016/j.chom.2020.12.001     #
##############################################

import os
import sys
import pandas as pd

DOMAINS = ['bacteria'] 
# DOMAINS = ['bacteria','archaea','viruses','archaea']

SCORES = ['unambiguous','reads','score']

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

	dataFrames = dict(zip(SCORES,([] for i in range(len(SCORES)))))

	missing_ct = 0 # Count missing files

	for fname in resultfiles:

		if not fname in os.listdir(OUT_DIR):
			missing_ct += 1
			continue

		filepath = os.path.join(OUT_DIR, fname)
		df = pd.read_table(filepath, index_col=0)

		for score in SCORES:
			dataFrames[score].append(df[score].rename(fname2id[fname]))

	if missing_ct > 0:
		print('WARNING: did not find pathseq output for {} of {} files'.format(
			missing_ct,len(resultfiles)))


	for domain in DOMAINS:
		for score in SCORES:
		
			data = pd.concat(dataFrames[score],axis=1)

			print("Joining {} {}".format(domain,score))
			ix_domain = taxa.index[taxa['kingdom'].isin([domain.title(),'root'])]
			# ix = np.intersect1d(data.index, )
			data_domain = data.loc[data.index.isin(ix_domain)]

			print("Resulting data table:",data_domain.shape)
			fpath = os.path.join(RESULT_DIR,'{}.{}.txt'.format(domain.lower(),score))
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


if __name__ == "__main__":

	try:
		PROJECT, ASSAY, minClippedReadLength = sys.argv[1:]
	except:
		print("ERROR: Missing assay and/or project")
		exit(1)

	print("Summarizing results for {} {}...".format(PROJECT, ASSAY))
	print("Min clipped read length = {}".format(minClippedReadLength))

	OUT_DIR = './pathseq_out/{}/{}'.format(PROJECT,ASSAY)
	RESULT_DIR = './results/{}/{}'.format(PROJECT,ASSAY)

	if not os.path.exists(RESULT_DIR): os.makedirs(RESULT_DIR)

	METAFILE = './manifests/gdc_manifest.{}.{}.txt'.format(ASSAY, PROJECT)
	TAXAFILE = './taxonomy/taxa.all.txt'

	meta = pd.read_table(METAFILE,index_col=0)
	taxa  = pd.read_table(TAXAFILE,index_col=0)
	fname2id, id2fbase = make_dicts(meta)

	outfiles = os.listdir(OUT_DIR)

	resultfiles = list(fname2id.keys())

	print("Checking {} files...".format(len(outfiles),len(meta)))
	parse_results(resultfiles)
	print("Done.")


