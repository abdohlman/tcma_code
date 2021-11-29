#!/usr/bin/env python

# USAGE: python combine_projects.py -a $assay -n TCMA -p COAD,READ,ESCA,STAD,HNSC 

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
	combine_data()
	combine_metadata()


def combine_metadata():
	'''Combines project metadata.'''

	fpath = './metadata/'

	for level in ['file','sample','case']:

		print("Combining {} metadata:".format(level)) #, projects_list)
		print("\n\tSamples x Features")

		meta_frames = []

		for project in projects_list:

			fname = 'metadata.{}.{}.txt'.format(project, level)
			meta = pd.read_table(fpath + fname, index_col=0, low_memory=False)

			print("{}\t{}".format(project, meta.shape))

			meta_frames.append(meta)

		meta_combined = pd.concat(meta_frames)

		print('--------------------')
		print("{}\t{}\n".format(combined_project_name, meta_combined.shape))

		fname = 'metadata.{}.{}.txt'.format(combined_project_name, level)
		meta_combined.to_csv(fpath + fname, sep='\t')


def combine_data():
	'''Combines project read counts and read statistics.'''


	data_frames = []
	read_stats = []

	print("Combining projects data:", projects_list)
	print("\n\tTaxa x Samples")

	for project in projects_list:

		fpath = './results/{}/{}/'.format(project, assay)
		fname = '{}.{}.txt'.format(domain, stat)

		df = pd.read_table(fpath + fname, index_col=0)
		data_frames.append(df)

		print("{}\t{}".format(project, df.shape))

		rs = pd.read_table(fpath + 'read_statistics.txt', index_col=0)
		read_stats.append(rs)

	df_combined = pd.concat(data_frames, axis=1)
	rs_combined = pd.concat(read_stats)

	print('--------------------')
	print("{}\t{}\n".format(combined_project_name, df_combined.shape))

	fpath = './results/{}/{}/'.format(combined_project_name, assay)
	if not os.path.exists(fpath): os.makedirs(fpath)

	fname = '{}.{}.txt'.format(domain, stat)
	df_combined.to_csv(fpath + fname, sep='\t')

	fname = 'read_statistics.txt'
	rs_combined.to_csv(fpath + fname, sep='\t')



def create_parser():

    parser = argparse.ArgumentParser(
        description="""Combines various projects.""")

    parser.add_argument(
	    '-a', '--assay', nargs='?', required=True,
		help="""TCGA experimental strategy (eg. WGS)""")

    parser.add_argument(
        '-n', '--combined-project-name', nargs='?', required=True,
        help="""Name for combined project""")

    parser.add_argument(
        '-l', '--projects-list', nargs='?', required=True,
        help="""TCGA sequencing projects to combine, comma separated
        (eg. COAD,READ,ESCA,STAD,HNSC)""")

    parser.add_argument(
        '-s', '--statistic', nargs='?', default='unambiguous.decontam',
        help="""Statistic to acquire from pathseq output. Options include
        "unambiguous.decontam" (default), "score", "reads".""")

    parser.add_argument(
        '-d', '--domain', nargs='?', default='bacteria',
        help="""Domain or kingdom to classify. Options include:
        "bacteria" (default), "archaea", "fungi", "viruses".""")

    return parser


if __name__ == "__main__":

	parser = create_parser()
	args = parser.parse_args()
	print(args)

	d = vars(args)

	combined_project_name = d['combined_project_name']
	projects_list = d['projects_list'].split(',')
	assay = d['assay']
	domain = d['domain']
	stat = d['statistic']

	main()
	print("Done.\n")
