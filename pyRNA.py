from threading import Thread
import argparse
import pysam
import pandas as pd
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from utils.utils import *
from utils.runj_module import *

def get_options():
	parser = argparse.ArgumentParser(description='What to do ?')
	parser.add_argument('--pathM', metavar="path_to_bam_of_analysis", type=str, 
						default=None, nargs='?')
	parser.add_argument('--pathC', metavar="path_to_control_bam", type=str, 
						default=None, nargs='?')
	parser.add_argument('-S', metavar='sample', nargs='?', default=None, type = str,
						help='sample name .bam')
	parser.add_argument('-V', metavar='sample_c', nargs='?', default=None, type = str,
						help='sample name .bam')
	parser.add_argument('-A', metavar='action', required=True, type=str)
	parser.add_argument('-L', metavar='window_len', nargs='?', default=100, type = int,
						help='Position input in int format')
	parser.add_argument('-O', metavar='output', nargs='?', default=None, type=str,
						help='output name') 
	parser.add_argument('-P', metavar='Position_input', nargs='?', default=None, type = int,
						help='Position input in int format')
	parser.add_argument('-E', metavar='Event', nargs='?', default=None, type = str,
						help='"AG", "AL", "DG", "DL"')
	parser.add_argument('-C', metavar='Chrom_input', nargs='?', default=None,type=str,
						help='Chrom input in the chrX format')   
	parser.add_argument('-G', metavar='gene', nargs='?', default=None, type = str,
						help='If known, add the spliceAI predicted geneName"')						
	parser.add_argument('-R', metavar='reference', required=True, type=str, 
						help='Reference of the genome (grch37 or grch38)')
	parser.add_argument('--ret', metavar="reticule", help='reticule on the print sequences module', default=False)
	args = parser.parse_args()
	return (args)

#######################################
#Transforms a Bam into Pandas DataFrame
#######################################

def load_bam(CHROM, POS, WINDOW, path_to_bam, SAMPLE): 
	'''
	CHROM : str
	POS : int (hg38 or hg19)
	WINDOW: int
	path_to_bam : str
	SAMPLE : str
	'''
#loads the bam
    try :
        samfile = pysam.AlignmentFile(path_to_bam+"/"+SAMPLE, "rb") 
    except :
        return (0,0,0) #ERROR CODE

 #gets read from the proper region
    list_reads = []
    for read in samfile.fetch(str(CHROM), POS-WINDOW, POS + WINDOW):
        list_reads = list_reads + [space(str(read))]
    loc = pd.DataFrame(list_reads)
    try :
        loc.columns = ['QNAME', 'FLAG', 'RNAME','pos_start', 'MAPQ', 'Cigar_String', 'RNEXT', 'PNEXT', 'TLEN', 'seq', 'QUAL', 'L']
    except :
        return (0,0, 0)  #ERROR CODE
    for p in range(len(loc)):
        loc.QUAL[p]=to_list(loc.QUAL[p][12:-2])

#helpers for zero-padding reads
    loc['temp'] = loc['Cigar_String'] + ['_']*len(loc) + loc['seq'] + ['&']*len(loc) + loc['QUAL'].map(str) 
    loc['true_seq'] = loc['temp'].map(true_seq)
    loc['has_inserts'] = loc['temp'].map(has_inserts) #corrections needed if there exist virals insertions
    del loc['temp']
    loc['loc_end'] = loc['pos_start'].map(int) + loc['true_seq'].map(len)

#Sequence alignment via zero-padding
    loc = align_seq(loc, POS-WINDOW-1, WINDOW)
    try :
        loc['coverage'] = loc['true_aligned_seq'].map(comptage)
        return (loc)
    except :
        return (1)
    return (loc)

  





args = get_options()

try:
	if (args.R=="grch37"):
		ref = pd.read_csv("grch37.txt", sep="\t")
	if (args.R=="grch38"):
		ref = pd.read_csv("grch38.txt", sep="\t")
except : 
	print("missing reference")

try : 

#######################################
# Functionality 1 : transfomation of a bam into a pandas Dataframe, no loss of information. Writes it as a tsv. 
#######################################

	if (args.A == 'reads_to_dataframe'): 
		try : 
			CHROM = args.C
		except : 
			print("missing chrom")
		try : 
			POS = args.P
		except : 
			print ("missing position")
		try :
			WINDOW = args.L
		except : 
			print ("missing window")
		try :
			path_to_bam = args.pathM
		except : 
			print ("missing path to bam")
		try :
			sample = args.S
		except : 
			print ("missing sample name")
		try : 
			EVENT = args.E
		except : 
			print("missing event argument")
		try : 
			gene = args.G
		except : 
			print("missing gene name argument")
		load_bam(CHROM, POS, WINDOW, path_to_bam, sample).to_csv(args.O+".tsv", sep="\t", index=False)


#######################################
# Functionality 2 : writes ONLY the reads of a region, aligned from POS-WINDOW to POS+WINDOW
#######################################
	elif (args.A == 'write_reads'): 
		try : 
			CHROM = args.C
		except : 
			print("missing chrom")
		try : 
			POS = args.P
		except : 
			print ("missing position")
		try :
			WINDOW = args.L
		except : 
			print ("missing window")
		try :
			path_to_bam = args.pathM
		except : 
			print ("missing path to bam")
		try :
			sample = args.S
		except : 
			print ("missing sample name")
		loc = load_bam(CHROM, POS, WINDOW, path_to_bam, sample)
		write_sequence(loc, WINDOW, sample, CHROM, POS, reticule=False)


	else : 

#######################################
# Functionality 3 : Splicing Analysis, computes 4 empirical RUNJ scores (cf ReadME) at a specific position
#######################################
		try : 
			CHROM = args.C
		except : 
			print("missing chrom")
		try : 
			POS = args.P
		except : 
			print ("missing position")
		try :
			WINDOW = args.L
		except : 
			print ("missing window")
		try :
			path_to_bam = args.pathM
		except : 
			print ("missing path to bam")
		try :
			sample = args.S
		except : 
			print ("missing sample name")
		try : 
			EVENT = args.E
		except : 
			print("missing event argument")
		try : 
			gene = args.G
		except : 
			print("missing gene name argument")
		ref = pd.read_csv(args.R, sep='\t')
		loc = load_bam(CHROM, POS, WINDOW, path_to_bam, sample)
		COV, Score_AG, Score_AL, Score_DG, Score_DL = check_splicing_events(loc, WINDOW, gene, ref, EVENT)
		print("*** Region coverage *** \n COV="+ str(COV) + "\n*** Splicing scores ***\n"+ "* AG = " + str(Score_AG) + "\n* AL = "+ str(Score_AL) + "\n* DG = " + str(Score_DG) + "\n* DL = "+ str(Score_DL) + "\n**********************")
except :
		ref = pd.read_csv(args.R, sep='\t')
		loc = load_bam(CHROM, POS, WINDOW, path_to_bam, sample)
		COV, Score_AG, Score_AL, Score_DG, Score_DL = check_splicing_events(loc, WINDOW, gene, ref, EVENT)
		print("*** Region coverage *** \n COV="+ str(COV) + "\n*** Splicing scores ***\n"+ "* AG = " + str(Score_AG) + "\n* AL = "+ str(Score_AL) + "\n* DG = " + str(Score_DG) + "\n* DL = "+ str(Score_DL) + "\n**********************")






