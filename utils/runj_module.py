from threading import Thread
import argparse
import pysam
from utils import *

import pandas as pd
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



def is_in_region(seq, window):
	'''
	seq : str
	window: int
	'''
	try :
		if ((seq[window-1]=="0") and (seq[window]=="0") and (seq[window+1]=="0")) :
			return False
		else :
			return True
	except :
		False

def get_strand(gene, ref):
	'''
	gene : str
	ref: pandas dataframe
	'''
	try :
		fo = ref[ref.NAME==str(gene)]
	except :
		fo = ref
	try :
		out = fo.STRAND.iloc[0]
	except :
		return ("?")
	if (len(fo)>1):
		return ("?")
	else :
		return (out)

# Creates the string of reads, separated with a \n
def sequences(loc, window, reticule):
	'''
	loc : pandas dataframe
	window: int
	reticule : bool
	'''
	out = ""
	for k in range (len(loc)):
		seq = loc.true_aligned_seq.iloc[k]
		if (is_in_region(seq, window)):
			if reticule:
				out = out +  seq[:window] + "||" + seq[window:window+1] + "||" +seq[window+1:] + "\n"
			else : 
				out = out +  seq + "\n"
	return (out)

# Writes down the reads as a text file
def write_sequence(loc, window, sample, chrom, var, reticule):
	'''
	loc : pandas dataframe
	window: int
	sample: str
	chrom : str
	var : position of variant(POS), int
	reticule : bool
	'''
	if isinstance(loc, tuple):
		print("sequence is empty")
	else : 
		reads = sequences(loc, window, reticule)
		r = open(str(sample) + "_" + str(chrom) + "_" + str(var)+'.txt', "w")
		r.write(reads)
		r.close()


def check_sequence_for_splicing_events(seq, window, strand, event):
	'''
	seq : str
	window: int
	strand: str
	event : str
	'''
	if (is_in_region(seq, window)):
		if (strand=="+"):
			if (event=="AG"):
				window -=1
				if  (is_null(seq[window]) and is_not_null(seq[window+1])): 
					return (1, 1, 0, 0, 0)
				else :
					return (1, 0, 0, 0, 0)
			if (event=="AL"): 
				if (is_null(seq[window]) and is_null(seq[window-1])): 
					return (1, 0, 1, 0, 0)
				else :
					return (1, 0, 0, 0, 0)
			if (event=="DG"):
				window +=1
				if (is_null(seq[window]) and is_not_null(seq[window-1])):
					return (1, 0, 0, 1, 0)
				else :
					return (1, 0, 0, 0, 0)
			if (event=="DL"): 
				if (is_not_null(seq[window]) and is_not_null(seq[window+1])):
					return (1, 0, 0, 0, 1)
				else :
					return (1, 0, 0, 0, 0)
		if (strand=="-"):
			if (event=="AG"): 
				window +=1
				if (is_null(seq[window]) and is_not_null(seq[window-1])):
					return (1, 1, 0, 0, 0)
				else :
					return (1, 0, 0, 0, 0)
			if (event=="AL"): 
				if (is_null(seq[window]) and is_null(seq[window+1])):
					return (1, 0, 1, 0, 0)
				else :
					return (1, 0, 0, 0, 0)
			if (event=="DG"): 
				window -=1
				if (is_null(seq[window]) and is_not_null(seq[window+1])):
					return (1, 0, 0, 1, 0)
				else :
					return (1, 0, 0, 0, 0)
			if (event=="DL"): 
				if (is_not_null(seq[window]) and is_not_null(seq[window-1])):
					return (1, 0, 0, 0, 1)
				else :
					return (1, 0, 0, 0, 0)
	else :
		return (0, 0, 0, 0, 0)

def check_splicing_events(loc, window, gene, ref, EVENT):
	'''
	loc : pandas DataFrame
	window: int
	gene: str
	ref : pandas DataFrame
	EVENT : str
	'''
	strand = get_strand(gene, ref)
	COV, AG, AL, DG, DL = 0,0,0,0,0
	for k in range (len(loc)):
		seq = loc.true_aligned_seq.iloc[k]
		cov, ag, al, dg, dl = check_sequence_for_splicing_events(seq, window, strand, EVENT)
		COV += cov 
		AG += ag
		AL += al 
		DG += dg 
		DL += dl 
	if (COV == 0):
		return (0, 0, 0, 0, 0)
	else : 
		return (COV, AG/COV, AL/COV, DG/COV, DL/COV)

