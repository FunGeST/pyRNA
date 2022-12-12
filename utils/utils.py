import argparse
import pysam
import pandas as pd
import os
import numpy as np
import random as rd
from threading import Thread


#####################################################################################
### Find strand from gene name + position 
#####################################################################################
gene_ref = ref

def get_strand(chrom, POS, gene):
	'''
	chrom : str
	POS : int
	gene : str
	'''
	global ref
	POS_start = int(POS) + int(50)
	POS_end = int(POS) + int(50)
	chrom = str(chrom).strip("chr")
	fo = ref[ref.CHROM!=chrom]
	fo = fo[fo.TX_START<POS_start]
	fo = fo[fo.TX_END>POS_end]
	try : 
		fo = fo[fo.NAME==str(gene)]
	except : 
		fo = fo
	try :	
		out = fo.STRAND.iloc[0]
	except : 
		return ("?")
	if (len(fo)>1):
		return ("?")
	else : 
		return (out)



#####################################################################################
### helpers for load_bam function 
#####################################################################################

#separates the data
def space(txt):
	start = 0
	stop = 0
	out = []
	for p in range(len(txt)):
		if txt[p:p+1]=="\t":
			stop = p
			out = out + [str(txt[start+1:stop])]
			start = stop
	out = out + [str(txt[start+1:])]
	return (out)

# transforms pysam.fecth into a list
def to_list(txt):
	out = [txt[0:2]]
	start = 2
	stop = 0
	for p in range (3, len(txt)):
		if (txt[p]==","):
			stop = p
			out = out + [txt[start+2:stop]]
			start = stop
	out = out + [txt[start+2:]]
	return (out)

#  reads the cigar string of a bam 
def recipe(cigar):
	sep = []
	for p in range(len(cigar)):
		if (cigar[p].isalpha()):
			sep = sep + [p]
	out = []
	start = 0
	stop = 0
	for k in range(len(sep)):
		stop = sep[k]
		out = out + [cigar[start:stop+1]]
		start = stop+1
	return (out)

# reconstructs the sequence from the CIGAR string
def reconstruct(cigar, seq, seq_card):
    out =""
    count = 0
    insertions = False
    seq_card_out = []
    mapp = recipe(cigar)
    for c in mapp:
        if (c[-1]=='M'):
            out = out+ seq[count:count+int(c[:-1])]
            seq_card_out = seq_card_out + [seq_card[count: count+int(c[:-1])]]
            count = count + int(c[:-1])
        elif (c[-1]=='N'):
            out = out + "-"*int(c[:-1])
            seq_card_out = seq_card_out + [0]*int(c[:-1])
        elif (c[-1]=='D'):
            out = out + "0"*int(c[:-1])
            seq_card_out = seq_card_out + [0]*int(c[:-1])
        elif (c[-1]=='I'):
            out = out + "0"*int(c[:-1])
            out = out + seq[count:count+int(c[:-1])]
        else :
            out = out + seq[count:count+int(c[:-1])]
            seq_card_out = seq_card_out + [seq_card[count: count+int(c[:-1])]]
            count = count + int(c[:-1])
            insertions = True
    return (out, insertions, seq_card)

#helper for separator 
def find_dash(txt):
	a = 0
	b = 0
	for p in range(len(txt)):
		if (txt[p]=='_'):
			a = p
		if (txt[p]=='&'):
			b = p
	return (txt[:a], txt[a+1:b], txt[b+1:])

#binds all of the above
def true_seq(temp):
	local = find_dash(temp)
	return(reconstruct(local[0], local[1], to_list(local[2]))[0])

# hecks if seq has insertions
def has_inserts(temp):
	local = find_dash(temp)
	return(reconstruct(local[0], local[1], to_list(local[2]))[1])

# aligns all the sequences with a starting position by addig zeroes when necessary
def align_seq(loc, start_pos, window):
    true_aligned_seq = []
    for p in range (len(loc)):
        delta = int(loc.pos_start[p]) - start_pos
        delta_2 = start_pos+2*window - int(loc.loc_end[p])

        true_primary_len =  int(loc.loc_end[p]) - int(loc.pos_start[p])
        if (delta < 0):

            out = loc.true_seq[p][abs(delta):]

        if (delta >0):
            out = "0"*(delta) + loc.true_seq[p]
        if (delta == 0):
            out = loc.true_seq[p]

        if (delta_2 <0):
            out = out[:delta_2]
        if(delta_2>0):
            out = out + "0"*(delta_2)
        true_aligned_seq = true_aligned_seq + [out]
    loc['true_aligned_seq'] = true_aligned_seq

    return (loc)



#####################################################################################
#Splicing analysis related utils
#####################################################################################


def is_null(txt):
    return (txt=="-") and (txt!="0")

def is_not_null(txt):
    return (txt!="-") and (txt!="0")

def get_first(txt):
	if (txt.count("_")):
		return (txt[:txt.index("_")])
	else :
		return (txt)
def get_second(txt):
	if (txt.count("_")):
		return (txt[txt.index("_")+1:])
	else :
		return (txt)
def separate_line(df, q):
	df = df.append(df.iloc[q,:])
	df.DS_max[q] = get_first(str(df.DS_max.iloc[q])).strip(" ")
	df.DS_max_kind[q] = get_first(str(df.DS_max_kind.iloc[q])).strip(" ")
	df.DS_max_pos[q] = get_first(str(df.DS_max_pos.iloc[q])).strip(" ")
	df.DS_max[-1] = get_second(str(df.DS_max.iloc[-1])).strip(" ")
	df.DS_max_kind[-1] = get_second(str(df.DS_max_kind.iloc[-1])).strip(" ")
	df.DS_max_pos[-1] = get_second(str(df.DS_max_pos.iloc[-1])).strip(" ")
	return(df)
def separate_full(df):
	for p in range (len(df)):
		if (str(df.DS_max[p]).count("_") or str(df.DS_max_pos[p]).count("_") or str(df.DS_max_kind[p]).count("_")):
			df = separate_line(df, p)
	return (df)


