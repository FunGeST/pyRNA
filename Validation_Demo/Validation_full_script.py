import argparse
import pysam
import pandas as pd
import os
import numpy as np
import random as rd
from threading import Thread

lim_runj = 0.01

path_to_bam = "utils/bams/"

parser = argparse.ArgumentParser(description='Variants_list')
parser.add_argument('inputfile', type=str, help='Input file in .vcf format')

args = parser.parse_args()

ref = pd.read_csv("utils/grch37.txt", sep="\t")
ref_projs = pd.read_csv("utils/RNA_project.txt", sep="\t")

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#                       UTILS
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################


#####################################################################################
### read the bam and transform it to a pandas data frame
#####################################################################################
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


#####################################################################################
### reconstruct the real sequence via the read and the CIGAR string
#####################################################################################
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

def find_dash(txt):
    a = 0
    b = 0
    for p in range(len(txt)):
        if (txt[p]=='_'):
            a = p
        if (txt[p]=='&'):
            b = p
    return (txt[:a], txt[a+1:b], txt[b+1:])

def true_seq(temp):
    local = find_dash(temp)
    return(reconstruct(local[0], local[1], to_list(local[2]))[0])

# check if seq has insertions
def has_inserts(temp):
    local = find_dash(temp)
    return(reconstruct(local[0], local[1], to_list(local[2]))[1])

def true_seq_card(temp):
    local = find_dash(temp)
    out = reconstruct(local[0], local[1], to_list(local[2]))[2]
    outta = []
    for c in out:
        try :
            pop = int(c.strip("'").strip(","))
            outta = outta + [pop]
        except :
            continue
    return (outta)

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

def comptage(txt):
    out = []
    for c in txt:
        if (c=="A") or (c=="C") or (c=="G") or (c=="T"):
            out = out + [1]
        else :
            out = out + [0]
    return (out)


#####################################################################################
### Binds all of the above. Full analysis for 1 variant
#####################################################################################
def load_bam(SAMPLE, PROJECT):
    try :
        samfile = pysam.AlignmentFile(path_to_bam+SAMPLE+".bam", "rb")
    except :
        samfile = (0,0,0)
    return(samfile)


def to_dataframe(samfile, CHROM, PROJECT, VAR, SAMPLE, DS, EVENT, WINDOW):
    if (isinstance(samfile, tuple)):
        return (0,0,0)
    else :
        list_reads = []
        for read in samfile.fetch(str(CHROM), VAR-WINDOW, VAR + WINDOW, until_eof=True):
            list_reads = list_reads + [str(read)]
        temp = []
        try :
            for p in range (math.ceil(len(list_reads)/2)-10, math.ceil(len(list_reads)/2)+10):
                u = space(list_reads[p])
                temp = temp + [u]
            loc = pd.DataFrame(temp)
            try :
                loc.columns = ['QNAME', 'FLAG', 'RNAME','pos_start', 'MAPQ', 'Cigar_String', 'G', 'H', 'I', 'seq', 'seq_card', 'L']
            except :
               loc = (0,0,0)
            return (loc)
        except :
            for p in range (len(list_reads)):
                temp = temp + [space(list_reads[p])]
            loc = pd.DataFrame(temp)
            try :
                loc.columns = ['QNAME', 'FLAG', 'RNAME','pos_start', 'MAPQ', 'Cigar_String', 'G', 'H', 'I', 'seq', 'seq_card', 'L']
            except :
               loc = (0,0,0)
            return (loc)

def formatage(loc, CHROM, PROJECT, VAR, SAMPLE, DS, EVENT, WINDOW):
    if (isinstance(loc, tuple)):
        return (0,0,0)
    else :
        for p in range(len(loc)):
            try:
                loc.seq_card[p]=to_list(loc.seq_card.iloc[p][12:-2])
            except :
                continue
        loc = loc.reset_index()
        del loc['index']
        try :
            loc['temp'] = loc['Cigar_String'] + ['_']*len(loc) + loc['seq'] + ['&']*len(loc) + loc['seq_card'].map(str)
            loc['true_seq'] = loc['temp'].map(true_seq)
            loc['has_inserts'] = loc['temp'].map(has_inserts)
            loc['true_seq_card'] = loc['temp'].map(true_seq_card)
            del loc['temp']
            loc['loc_end'] = loc['pos_start'].map(int) + loc['true_seq'].map(len)
            loc = align_seq(loc, VAR+DS-WINDOW-1, WINDOW)
            return (loc)
        except :
            return (0,0,0)


def load_and_do_stuff(CHROM, VAR, SAMPLE, PROJECT, DS, EVENT, WINDOW):
    if ((SAMPLE==SAMPLE)):
        samfile = pysam.AlignmentFile(path_to_bam+SAMPLE+".bam", "rb")
    else :
        return (0,0,0)
    list_reads = []
    for read in samfile.fetch(str(CHROM), VAR-WINDOW, VAR + WINDOW):
        list_reads = list_reads + [space(str(read))]
    loc = pd.DataFrame(list_reads)
    try :
        loc.columns = ['QNAME', 'FLAG', 'RNAME','pos_start', 'MAPQ', 'Cigar_String', 'G', 'H', 'I', 'seq', 'seq_card', 'L']
    except :
        return (0,0, 0)
    for p in range(len(loc)):
        loc.seq_card[p]=to_list(loc.seq_card[p][12:-2])
    loc['temp'] = loc['Cigar_String'] + ['_']*len(loc) + loc['seq'] + ['&']*len(loc) + loc['seq_card'].map(str)
    loc['true_seq'] = loc['temp'].map(true_seq)
    loc['has_inserts'] = loc['temp'].map(has_inserts)
    loc['true_seq_card'] = loc['temp'].map(true_seq_card)
    del loc['temp']
    loc['loc_end'] = loc['pos_start'].map(int) + loc['true_seq'].map(len)
    loc = align_seq(loc, VAR+DS-WINDOW-1, WINDOW)
    return (loc)


def is_in_region(seq, POS, var, window):
    try :
        if ((seq[window-1]=="0") or (seq[window]=="0") or (seq[window+1]=="0")) :
            return False
        else :
            return True
    except :
        False

def is_null(txt):
    return (txt=="-") and (txt!="0")

def is_not_null(txt):
    return (txt!="-") and (txt!="0")

def validate_event(seq, POS, var, window, event, way):
    if (is_in_region(seq, POS, var, window)):
        if (way=="+"):
            if (event=="AG"): 
                window -=1
                if  (is_null(seq[window]) and is_not_null(seq[window+1])): 
                    return (1, 1, seq[:window] + "||" + seq[window:window+2] + "||" +seq[window+2:] + "\n\n Validated  \n\n")
                else :
                    return (1, 0, seq[:window] + "||" + seq[window:window+2] + "||" +seq[window+2:])
            if (event=="DL"): 
                if (is_null(seq[window]) and is_null(seq[window-1])): 
                    return (1, 1, seq[:window-1] + "||" + seq[window-1:window+1] + "||" +seq[window+1:]+ "\n\n Validated \n\n")
                else :
                    return (1, 0,  seq[:window-1] + "||" + seq[window-1:window+1] + "||" +seq[window+1:])
            if (event=="DG"): 
                window +=1
                if (is_null(seq[window]) and is_not_null(seq[window-1])):
                    return (1, 1,  seq[:window-1] + "||" + seq[window-1:window+1] + "||" +seq[window+1:]+ "\n \nValidated \n\n")
                else :
                    return (1, 0,  seq[:window-1] + "||" + seq[window-1:window+1] + "||" +seq[window+1:])
            if (event=="AL"): 
                if (is_not_null(seq[window]) and is_not_null(seq[window+1])):
                    return (1, 1, seq[:window] + "||" + seq[window:window+2] + "||" +seq[window+2:]+ "\n \nValidated \n\n")
                else :
                    return (1, 0, seq[:window] + "||" + seq[window:window+2] + "||" +seq[window+2:])
        if (way=="-"):
            if (event=="AG"): 
                window +=1
                if (is_null(seq[window]) and is_not_null(seq[window-1])):
                    return (1, 1,  seq[:window-1] + "||" + seq[window-1:window+1] + "||" +seq[window+1:]+ "\n\nValidated \n\n")
                else :
                    return (1, 0,  seq[:window-1] + "||" + seq[window-1:window+1] + "||" +seq[window+1:])
            if (event=="DL"): 
                if (is_null(seq[window]) and is_null(seq[window+1])):
                    return (1, 1, seq[:window] + "||" + seq[window:window+2] + "||" +seq[window+2:]+ "\n \n Validated \n\n")
                else :
                    return (1, 0, seq[:window] + "||" + seq[window:window+2] + "||" +seq[window+2:])
            if (event=="DG"): 
                window -=1
                if (is_null(seq[window]) and is_not_null(seq[window+1])):
                    return (1, 1, seq[:window] + "||" + seq[window:window+2] + "||" +seq[window+2:]+ "\n \n Validated \n\n")
                else :
                    return (1, 0, seq[:window] + "||" + seq[window:window+2] + "||" +seq[window+2:]+ "\n \n")
            if (event=="AL"): 
                if (is_not_null(seq[window]) and is_not_null(seq[window-1])):
                    return (1, 1,  seq[:window-1] + "||" + seq[window-1:window+1] + "||" +seq[window+1:]+ "\n \n Validated \n\n")
                else :
                    return (1, 0,  seq[:window-1] + "||" + seq[window-1:window+1] + "||" +seq[window+1:])
    else :
        return (0,0, "-")


def validation(loc, POS, var, window, event, strand):
    cov = 0
    val = 0
    reads = "strand = " + str(strand) +"\n"+"\n"
    if (strand=="?"):
        a_p, b_p, a_m, b_m = 0, 0, 0, 0
        for p in range(len(loc)):
            try :
                if (int(loc.MAPQ.iloc[p])<20):
                    continue
            except :
                pass
            a, b, c = validate_event(loc.true_aligned_seq.iloc[p], POS, var, window, event, "+")
            a_p += a
            b_p += b
            a, b, c_m = validate_event(loc.true_aligned_seq.iloc[p], POS, var, window, event, "-")
            a_m += a
            b_m += b
        if (max([b_m,b_p])==b_m):
            cov = a_m
            val = b_m
            strand = "-"
        else :
            cov = a_p
            val = b_p
            strand = "+"
        return (cov, val, "Unkown strand")
    else :
        for p in range(len(loc)):
            try :
                if (int(loc.MAPQ.iloc[p])<20):
                    continue
            except :
                pass
            if (strand=="+"):
                a, b,c = validate_event(loc.true_aligned_seq.iloc[p], POS, var, window, event, "+")
            if (strand=="-"):
                a, b, c = validate_event(loc.true_aligned_seq.iloc[p], POS, var, window, event, "-")
            cov += a
            val += b
            reads += c + "\n"
        return (cov, val, reads)

def validation_control(loc, POS, var, window, event, strand):
    cov = 0
    val = 0
    reads = ""
    if (strand=="?"):
        a_p, b_p, a_m, b_m = 0, 0, 0, 0
        for p in range(len(loc)):
            try :
                if (int(loc.MAPQ.iloc[p])<20):
                   continue
            except :
                print(" ")
            a, b, c = validate_event(loc.true_aligned_seq.iloc[p], POS, var, window, event, "+")
            a_p += a
            b_p += b
            a, b, c = validate_event(loc.true_aligned_seq.iloc[p], POS, var, window, event, "-")
            a_m += a
            b_m += b
        if (min([b_m,b_p])==b_m):
            cov = a_m
            val = b_m
            strand = "-"
        else :
            cov = a_p
            val = b_p
            strand = "+"
        return (cov, val)
    else :
        for p in range(len(loc)):
            try :
                if (int(loc.MAPQ.iloc[p])<20):
                    continue
            except :
                print(" ")
            try : 
                if (strand=="+"):
                    a, b, c = validate_event(loc.true_aligned_seq.iloc[p], POS, var, window, event, "+")
                if (strand=="-"):
                    a, b, c = validate_event(loc.true_aligned_seq.iloc[p], POS, var, window, event, "-")
                cov += a
                val += b
                reads += c + "\n"
            except : 
                continue
        return (cov, val)

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

def get_strand(gene):
    global ref
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

def control_complet(chrom, var, control, PROJECT_c, DS, event, window, gene):
    POS = int(var) + int(DS)
    loc_c = load_and_do_stuff(chrom, var, control, PROJECT_c, DS, event, window)
    strand = get_strand(gene)
    if (strand=="?"):
        cov_c, val_c, strand= validation_control(loc_c, POS, var, window, event)
    else :
        cov_c, val_c = validation(loc, POS, var, window, event, strand)
    if (cov_c==0):
        raise NameError('Need a new control')
    else :
        return (cov_c, val_c, strand)

def find_sequencage(id):
    loc = ref_projs[ref_projs.RNAseqProject==str(id)].reset_index()
    try :
        return (loc.read.iloc[0])
    except :
        return ("?")

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#               The following parts actually validate predictions
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################


def utilitary_of_validation(df, p, window=5):
    global ref_samples_75
    global ref_samples_100
    chrom = str(df.CHROM.iloc[p])
    var = df.POS_hg19.iloc[p]
    sample = str(p+1)+"_Sliced_"+df['sample'].iloc[p]
    PROJECT = None
    gene = df.SpliceAI_SYMBOL.iloc[p]
    DS = int(float(df.DS_max_pos.iloc[p]))
    POS = int(var) + int(DS)
    strand = get_strand(gene)
    event = df.DS_max_kind.iloc[p]
    loc = load_and_do_stuff(chrom, var, sample, PROJECT, DS, event, window)
    if isinstance(loc, tuple):
        cov, val = -1, 0
        return(cov, val)
    else :
        AL, DL, DG, AG = 0, 0, 0, 0
        if (event=="AL"):
            cov, AL, reads = validation(loc, POS, var, window, "AL", strand)
            cov, DL, reads2 = validation(loc, POS, var, window, "DL", strand)
        elif (event=="DL"):
            cov, DL, reads = validation(loc, POS, var, window, "DL", strand)
            cov, AL, reads2 = validation(loc, POS, var, window, "AL", strand)
        elif (event=="AG"):
            cov, AG, reads = validation(loc, POS, var, window, "AG", strand)
        elif (event=="DL"):
            cov, DG, reads = validation(loc, POS, var, window, "DG", strand)
        r = open(str(sample) + "_" + str(chrom) + "_" + str(var) + "_" + str(event)+'.txt', "w")
        r.write(reads)
        r.close()
        return(cov, AL, DL, DG, AG)

def utilitary_of_validation_POT(df, p, samples_control_100_POT, window=5):
    try :
            chrom = str(int(float(df.CHROM[p])))
    except:
            chrom = str(df.CHROM[p])
    var = df.POS_hg19.iloc[p]
    sample = str(p+1)+"_Sliced_"+df['sample'].iloc[p]
    PROJECT = ""
    gene = df.SpliceAI_SYMBOL.iloc[p]
    DS = int(float(df.DS_max_pos.iloc[p]))
    POS = int(var) + int(DS)
    strand = get_strand(gene)
    event = df.DS_max_kind.iloc[p]

    CONTROL = pd.DataFrame()
    for sample_c in samples_control_100_POT:
        fraction_control = to_dataframe(sample_c, chrom, PROJECT, var, sample, DS, event, window)
        if (isinstance(fraction_control, tuple)):
            continue
        else :
            CONTROL = CONTROL.append(fraction_control)
    CONTROL = formatage(CONTROL, chrom, "ooo", var, sample, DS, event, window)
    AL_c, DL_c, DG_c, AG_c = 0, 0, 0, 0
    if (event=="AL"):
        cov_c, DL_c = validation_control(CONTROL, POS, var, window, "DL", strand)
        cov_c, AL_c = validation_control(CONTROL, POS, var, window, "AL", strand)
            #r = open("Control_"+str(sample) + "_" + str(chrom) + "_" + str(var) + "_" + str(event)+'_AL.txt', "w")
            #r.write(reads)
            #r.close()
    elif (event=="DL"):
        cov_c, DL_c = validation_control(CONTROL, POS, var, window, "DL", strand)
        cov_c, AL_c = validation_control(CONTROL, POS, var, window, "AL", strand)
            #r = open("Control_"+str(sample) + "_" + str(chrom) + "_" + str(var) + "_" + str(event)+'_DL.txt', "w")
            #r.write(reads)
            #r.close()
    elif (event=="DG"):
        cov_c, DG_c = validation_control(CONTROL, POS, var, window, "DG", strand)
    elif (event=="AG"):
        cov_c, AG_c = validation_control(CONTROL, POS, var, window, "AG", strand)

    return(cov_c, AL_c, DL_c, DG_c, AG_c)



def utilitary_of_validation_PON(df, p, samples_control_100_PON, window=5):
    try :
            chrom = str(int(float(df.CHROM[p])))
    except:
            chrom = str(df.CHROM[p])
    var = df.POS_hg19.iloc[p]
    sample = str(p+1)+"_Sliced_"+df['sample'].iloc[p]
    PROJECT = ""
    gene = df.SpliceAI_SYMBOL.iloc[p]
    DS = int(float(df.DS_max_pos.iloc[p]))
    POS = int(var) + int(DS)
    strand = get_strand(gene)
    event = df.DS_max_kind.iloc[p]

    CONTROL = pd.DataFrame()
    for sample_c in samples_control_100_PON:
        fraction_control = to_dataframe(sample_c, chrom, PROJECT, var, sample, DS, event, window)
        if (isinstance(fraction_control, tuple)):
            continue
        else :
            CONTROL = CONTROL.append(fraction_control)
    CONTROL = formatage(CONTROL, chrom, "ooo", var, sample, DS, event, window)
    AL_c, DL_c, DG_c, AG_c = 0, 0, 0, 0
    if (event=="AL"):
        cov_c, DL_c = validation_control(CONTROL, POS, var, window, "DL", strand)
        cov_c, AL_c = validation_control(CONTROL, POS, var, window, "AL", strand)
            #r = open("Control_"+str(sample) + "_" + str(chrom) + "_" + str(var) + "_" + str(event)+'_AL.txt', "w")
            #r.write(reads)
            #r.close()
    elif (event=="DL"):
        cov_c, DL_c = validation_control(CONTROL, POS, var, window, "DL", strand)
        cov_c, AL_c = validation_control(CONTROL, POS, var, window, "AL", strand)
        #r = open("Control_"+str(sample) + "_" + str(chrom) + "_" + str(var) + "_" + str(event)+'_DL.txt', "w")
        #r.write(reads)
        #r.close()
    elif (event=="DG"):
        cov_c, DG_c = validation_control(CONTROL, POS, var, window, "DG", strand)
    elif (event=="AG"):
        cov_c, AG_c = validation_control(CONTROL, POS, var, window, "AG", strand)

    return(cov_c, AL_c, DL_c, DG_c, AG_c)


def uni_thread_validation(df, p, samples_control_100_POT, samples_control_100_PON):
    global meta
    df_loc = df.iloc[p,:]
    try:
        cov, AL, DL, DG, AG= utilitary_of_validation(df, p, window=5)
    except :
        cov, AL, DL, DG, AG = -2, 0, 0, 0, 0#error_code
    df_loc["Coverage"] = cov
    if (df.DS_max_kind.iloc[p]=='AL'):
        df_loc["MA_COUNT"] = AL
        df_loc["all_effects_MA_COUNT"] = AL + DL
    elif (df.DS_max_kind.iloc[p]=='DL'):
        df_loc["MA_COUNT"] = DL
        df_loc["all_effects_MA_COUNT"] = AL + DL
    elif (df.DS_max_kind.iloc[p]=='DG'):
        df_loc["MA_COUNT"] = DG
    elif (df.DS_max_kind.iloc[p]=='AG'):
        df_loc["MA_COUNT"] = AG
    try:
        cov_c, AL_c, DL_c, DG_c, AG_c = utilitary_of_validation_POT(df, p, samples_control_100_POT, window=5)
    except :
        cov_c, AL_c, DL_c, DG_c, AG_c = "NA", "NA", "NA", "NA", "NA"
    df_loc["COV_C_POT"] = cov_c
    if (df.DS_max_kind.iloc[p]=='AL'):
        df_loc["MA_COUNT_POT"] = AL_c
        df_loc["all_effects_MA_COUNT_POT"] = AL_c + DL_c
    elif (df.DS_max_kind.iloc[p]=='DL'):
        df_loc["MA_COUNT_POT"] = DL_c
        df_loc["all_effects_MA_COUNT_POT"] = AL + DL_c
    elif (df.DS_max_kind.iloc[p]=='DG'):
        df_loc["MA_COUNT_POT"] = DG_c
    elif (df.DS_max_kind.iloc[p]=='AG'):
        df_loc["MA_COUNT_POT"] = AG_c
    df_loc['AL_MA_COUNT_c_POT'] = AL_c
    df_loc['DL_MA_COUNT_c_POT'] = DL_c
    df_loc['AG_MA_COUNT_c_POT'] = AG_c
    df_loc['DG_MA_COUNT_c_POT'] = DG_c
    try :
        df_loc["UNJ_AL_c_POT"] = AL_c/cov_c
    except :
        df_loc["UNJ_AL_c_POT"] = "NA"
    try :
        df_loc["UNJ_DL_c_POT"] = DL_c/cov_c
    except :
        df_loc["UNJ_DL_c_POT"] = "NA"
    try :
        df_loc["UNJ_AG_c_POT"] = AG_c/cov_c
    except :
        df_loc["UNJ_AG_c_POT"] = "NA"
    try :
        df_loc["UNJ_DG_c_POT"] = DG_c/cov_c
    except :
        df_loc["UNJ_DG_c_POT"] = "NA"

    #try:
    cov_c, AL_c, DL_c, DG_c, AG_c = utilitary_of_validation_PON(df, p, samples_control_100_PON, window=5)
    #except :
    #    cov_c, AL_c, DL_c, DG_c, AG_c = "NA", "NA", "NA", "NA", "NA"
    df_loc["COV_C_PON"] = cov_c
    if (df.DS_max_kind.iloc[p]=='AL'):
        df_loc["MA_COUNT_PON"] = AL_c
        df_loc["all_effects_MA_COUNT_PON"] = AL_c + DL_c
    elif (df.DS_max_kind.iloc[p]=='DL'):
        df_loc["MA_COUNT_PON"] = DL_c
        df_loc["all_effects_MA_COUNT_PON"] = AL + DL_c
    elif (df.DS_max_kind.iloc[p]=='DG'):
        df_loc["MA_COUNT_PON"] = DG_c
    elif (df.DS_max_kind.iloc[p]=='AG'):
        df_loc["MA_COUNT_PON"] = AG_c
    df_loc['AL_MA_COUNT_c_PON'] = AL_c
    df_loc['DL_MA_COUNT_c_PON'] = DL_c
    df_loc['AG_MA_COUNT_c_PON'] = AG_c
    df_loc['DG_MA_COUNT_c_PON'] = DG_c
    try :
        df_loc["UNJ_AL_c_PON"] = AL_c/cov_c
    except :
        df_loc["UNJ_AL_c_PON"] = "NA"
    try :
        df_loc["UNJ_DL_c_PON"] = DL_c/cov_c
    except :
        df_loc["UNJ_DL_c_PON"] = "NA"
    try :
        df_loc["UNJ_AG_c_PON"] = AG_c/cov_c
    except :
        df_loc["UNJ_AG_c_PON"] = "NA"
    try :
        df_loc["UNJ_DG_c_PON"] = DG_c/cov_c
    except :
        df_loc["UNJ_DG_c_PON"] = "NA"
    meta = meta.append(df_loc)
    return(None)



class multi_threaded_validation(Thread):
    def __init__(self, df, p):
        Thread.__init__(self)
        self.df = df
        self.p = p
        samples_control_POT = [str(p+1)+"_Sliced_"+ref_projs['sample'].iloc[int(index)] for index in ref_samples_100]
        samples_control_100_POT = []
        for index in samples_control_POT :
            samples_control_100_POT = samples_control_100_POT + [load_bam(index, "")]
        self.samples_control_100_POT = samples_control_100_POT
        ref_PON = ['Sliced_TCGA-BC-A10Q-11A','Sliced_TCGA-BC-A10R-11A','Sliced_TCGA-BC-A10T-11A','Sliced_TCGA-BC-A10U-11A', 'Sliced_TCGA-BC-A10W-11A']
        samples_control_100_PON = []
        for name in ref_PON :
            samples_control_100_PON = samples_control_100_PON + [load_bam(str(p+1)+"_"+name, "")]
        self.samples_control_100_PON = samples_control_100_PON

    def run(self):
        uni_thread_validation(self.df, self.p, self.samples_control_100_POT, self.samples_control_100_PON)

def valid_full(df, window=5):
    df['POS_hg19']=df['POS_hg19'].map(int)
    df['DS_max_pos']=df['DS_max_pos'].map(float).map(int)
    threads = []
    df = df.reset_index()
    for ind in df.index:
        threads.append(multi_threaded_validation(df, ind))
    for t in threads:
        t.start()
        t.join()
    return (meta)

def main(name):
    df = pd.read_excel(name)
    df = valid_full(df, window=5)
    return (df)


########################
#POT  CREATOR
########################
ref_projs = pd.read_excel("Splice_mutations_in_TCGA_samples_with_bam_files.xlsx", sep="\t")['sample'].drop_duplicates()
ref_projs = ref_projs.reset_index()
ref_samples_100 = [1, 6, 3, 4, 5]


#################################################
#################################################

# Parsing RNA-seq data and fetching key metrics

#################################################
#################################################
meta = pd.DataFrame()
loc = main(args.inputfile)


#################################################
#################################################

# Computation of scores and cleaning up 

#################################################
#################################################


loc = loc.drop_duplicates()
loc = loc.reset_index()
del loc['index']

loc['Coverage_mutated_sample'] = loc['Coverage']

loc['utils'] = loc['DS_max_kind'].map(str) + ["_"]*len(loc) + loc['MA_COUNT'].map(str) + ["&"]*len(loc) + loc['all_effects_MA_COUNT'].map(str)
loc['utils_c'] = loc['DS_max_kind'].map(str) + ["_"]*len(loc) + loc['MA_COUNT_PON'].map(str) + ["&"]*len(loc) + loc['AL_MA_COUNT_c_PON'].map(str) + ["?"]*len(loc) + loc['DL_MA_COUNT_c_PON'].map(str)

def RT_intron(txt):
    if (txt[:txt.index("_")]=="DL"):
        return(int(float(txt[txt.index("_")+1:txt.index('&')])))
    elif (txt[:txt.index("_")]=="AL"):
        return(int(float(txt[txt.index('&')+1:])) - int(float(txt[txt.index("_")+1:txt.index('&')])))
    else :
        return ("-")
def ES(txt):
    if (txt[:txt.index("_")]=="AL"):
        return(int(float(txt[txt.index("_")+1:txt.index('&')])))
    elif (txt[:txt.index("_")]=="DL"):
        return(int(float(txt[txt.index('&')+1:])) - int(float(txt[txt.index("_")+1:txt.index('&')])))
    else :
        return ("-")
def NS(txt):
    if (txt[txt.index("_")-1:txt.index("_")]=="G"):
        return (int(float(txt[txt.index("_")+1:txt.index('&')])))
    else :
        return ("-")
def RT_intron_c(txt):
    try :
        if (txt[:txt.index("_")]=="DL"):
            return(int(float(txt[txt.index("_")+1:txt.index('&')])))
        elif (txt[:txt.index("_")]=="AL"):
            return(int(float(txt[txt.index("&")+1:txt.index('?')])))
        else :
            return ("-")
    except :
        return ("-")
def ES_c(txt):
    try :
        if (txt[:txt.index("_")]=="AL"):
            return(int(float(txt[txt.index("_")+1:txt.index('&')])))
        elif (txt[:txt.index("_")]=="DL"):
            return(int(float(txt[txt.index("?")+1:])))
        else :
            return ("-")
    except :
        return ("-")
def NS_c(txt):
    if (txt[txt.index("_")-1:txt.index("_")]=="G"):
        return (int(float(txt[txt.index("_")+1:txt.index('&')])))
    else :
        return ("-")





#FOR PON ONLY

loc['RT_intron_PON'] = loc['utils'].map(RT_intron)
loc['RT_intron_PON_c'] = loc['utils_c'].map(RT_intron_c)
loc['ES_PON'] = loc['utils'].map(ES)
loc['ES_PON_c'] = loc['utils_c'].map(ES_c)
loc['NS_PON'] = loc['utils'].map(NS)
loc['NS_PON_c'] = loc['utils_c'].map(NS_c)

#####
loc['utils'] = loc['DS_max_kind'].map(str) + ["_"]*len(loc) + loc['MA_COUNT'].map(str) + ["&"]*len(loc) + loc['all_effects_MA_COUNT'].map(str)
loc['utils_c'] = loc['DS_max_kind'].map(str) + ["_"]*len(loc) + loc['AL_MA_COUNT_c_POT'].map(str) + ["&"]*len(loc) + loc['DL_MA_COUNT_c_POT'].map(str) + ["?"]*len(loc) + loc['AG_MA_COUNT_c_POT'].map(str) + ["!"]*len(loc) + loc['DG_MA_COUNT_c_POT'].map(str)

def RT_intron(txt):
    try :
        if (txt[:txt.index("_")]=="DL"):
            return(int(float(txt[txt.index("_")+1:txt.index('&')])))
        elif (txt[:txt.index("_")]=="AL"):
            return(int(float(txt[txt.index('&')+1:])) - int(float(txt[txt.index("_")+1:txt.index('&')])))
        else :
            return ("-")
    except:
        return ("-")
def ES(txt):
    try:
        if (txt[:txt.index("_")]=="AL"):
            return(int(float(txt[txt.index("_")+1:txt.index('&')])))
        elif (txt[:txt.index("_")]=="DL"):
            return(int(float(txt[txt.index('&')+1:])) - int(float(txt[txt.index("_")+1:txt.index('&')])))
        else :
            return ("-")
    except :
        return ("-")
def NS(txt):
    try :
        if (txt[txt.index("_")-1:txt.index("_")]=="G"):
            return (int(float(txt[txt.index("_")+1:txt.index('&')])))
        else :
            return ("-")
    except:
        return ("-")
def RT_intron_c(txt):
    try:
        if (txt[:txt.index("_")]=="DL"):
            return(int(float(txt[txt.index("&")+1:txt.index('?')])))
        elif (txt[:txt.index("_")]=="AL"):
            return(int(float(txt[txt.index("&")+1:txt.index('?')])))
        else :
            return ("-")
    except:
        return("-")
def ES_c(txt):
    try :
        if (txt[:txt.index("_")]=="AL"):
            return(int(float(txt[txt.index("_")+1:txt.index('&')])))
        elif (txt[:txt.index("_")]=="DL"):
            returnreturn(int(float(txt[txt.index("_")+1:txt.index('&')])))
        else :
            return ("-")
    except :
        return ("-")
def NS_c(txt):
    try :
        if (txt[:txt.index("_")]=="DG"):
            return (int(float(txt[txt.index("!")+1:])))
        elif (txt[:txt.index("_")]=="AG"):
            return (int(float(txt[txt.index("?")+1:txt.index('!')])))
        else :
            return ("-")
    except :
        return ("-")

#FOR POT ONLY

loc['RT_intron_POT'] = loc['utils'].map(RT_intron)
loc['RT_intron_POT_c'] = loc['utils_c'].map(RT_intron_c)
loc['ES_POT'] = loc['utils'].map(ES)
loc['ES_POT_c'] = loc['utils_c'].map(ES_c)
loc['NS_POT'] = loc['utils'].map(NS)
loc['NS_POT_c'] = loc['utils_c'].map(NS_c)


#################################################################################################
#################################################################################################

#################################################################################################
#################################################################################################

#################################################################################################
#################################################################################################
out = loc
def RUNJ_PON(df):
    out = []
    for p in range(len(df)):
        if (float(df.Coverage_mutated_sample[p])==0 or float(df.COV_C_PON[p])==0):
            val=0
        else :
            if (df.DS_max_kind[p]=="AG") or (df.DS_max_kind[p]=="DG"):
                val = float(df.NS_PON[p]) / float(df.Coverage_mutated_sample[p]) - float(df.NS_PON_c[p]) / float(df.COV_C_PON[p])
            elif (df.DS_max_kind[p]=="AL"):
                val = float(df.ES_PON[p]) / float(df.Coverage_mutated_sample[p]) - float(df.ES_PON_c[p]) / float(df.COV_C_PON[p])
            elif (df.DS_max_kind[p]=="DL"):
                val = float(df.RT_intron_PON[p]) / float(df.Coverage_mutated_sample[p]) - float(df.RT_intron_PON_c[p]) / float(df.COV_C_PON[p])
            else :
                val = -1
            if (val !="-"):
                if (float(val)<0):
                    val = -2
        out = out + [val]
    df["RUNJ_PON"] = out
    return (df)

out = RUNJ_PON(out)

def RUNJ_all_effects_PON(df):
    out = []
    for p in range(len(df)):
        if (float(df.Coverage_mutated_sample[p])==0 or float(df.COV_C_PON[p])==0):
            val=0
        elif (True):
            if (df.DS_max_kind[p]=="AG") or (df.DS_max_kind[p]=="DG"):
                val = float(df.NS_PON[p]) / float(df.Coverage_mutated_sample[p]) - float(df.NS_PON_c[p]) / float(df.COV_C_PON[p])
            elif (df.DS_max_kind[p]=="AL") or (df.DS_max_kind[p]=="DL"):
                try :
                    val = float((float(df.ES_PON[p]) + float(df.RT_intron_PON[p]))) / float(df.Coverage_mutated_sample[p]) - float((float(df.ES_PON_c[p]) +float(df.RT_intron_PON_c[p]))) / float(df.COV_C_PON[p])
                except:
                    val = 0
            else :
                val =-1
        else :
            print("Cas impossible")
        out = out + [val]
    df["RUNJ_all_effects_PON"] = out
    return (df)

out = RUNJ_all_effects_PON(out)

def RUNJ_POT(df):
    out = []
    for p in range(len(df)):
        if (float(df.Coverage_mutated_sample[p])==0 or float(df.COV_C_POT[p])==0):
            val=0
        else :
            try :
                if (df.DS_max_kind[p]=="DG") or (df.DS_max_kind[p]=="AG"):
                    val = float(df.NS_POT[p]) / float(df.Coverage_mutated_sample[p]) - float(df.NS_POT_c[p]) / float(df.COV_C_POT[p])
                elif (df.DS_max_kind[p]=="AL"):
                    val = float(df.ES_POT[p]) / float(df.Coverage_mutated_sample[p]) - float(df.ES_POT_c[p]) / float(df.COV_C_POT[p])
                elif (df.DS_max_kind[p]=="DL"):
                    val = float(df.RT_intron_POT[p]) / float(df.Coverage_mutated_sample[p]) - float(df.RT_intron_POT_c[p]) / float(df.COV_C_POT[p])
                else :
                    val =-1
            except :
                val = -1
        out = out + [val]
    df["RUNJ_POT"] = out
    return (df)

out = RUNJ_POT(out)

def RUNJ_all_effects_POT(df):
    out = []
    for p in range(len(df)):
        if (float(df.Coverage_mutated_sample[p])==0 or float(df.COV_C_POT[p])==0):
            val=0
        elif (True):
            try :
                if (df.DS_max_kind[p]=="DG") or (df.DS_max_kind[p]=="AG"):
                    val = float(df.NS_POT[p]) / float(df.Coverage_mutated_sample[p]) - float(df.NS_POT_c[p]) / float(df.COV_C_POT[p])
                elif (df.DS_max_kind[p]=="AL") or (df.DS_max_kind[p]=="DL"):
                    try :
                        val = float((float(df.ES_POT[p]) + float(df.RT_intron_POT[p]))) / float(df.Coverage_mutated_sample[p]) - float((float(df.ES_POT_c[p]) +float(df.RT_intron_POT_c[p]))) / float(df.COV_C_POT[p])
                    except:
                        val = 0
                else :
                    val =-1
            except:
                val = -1
        else :
            print("Cas impossible")
        out = out + [val]
    df["RUNJ_all_effects_POT"] = out
    return (df)

out = RUNJ_all_effects_POT(out)


def validated_RUNJ_PON(df):
    global lim_runj
    out = []
    for p in range (len(df)):
        if (float(df.Coverage_mutated_sample[p])==float(-1)):
            out = out + ["NC"]
        elif (float(df.Coverage_mutated_sample[p])==float(-2)):
            out = out + ["NC"]
        elif (float(df.Coverage_mutated_sample[p])<float(20)):
            out = out + ["NC"]
        elif (float(df.RUNJ_PON[p])>float(lim_runj)):
            out = out +  [True]
        else:
            out = out +  [False]
    df['RUNJ_validation_PON'] = out
    return (df)

out = validated_RUNJ_PON(out)
def validated_all_effects_RUNJ_PON(df):
    global lim_runj
    out = []
    for p in range (len(df)):
        if (float(df.Coverage_mutated_sample[p])==float(-1)):
            out = out + ["NC"]
        elif (df.RUNJ_all_effects_PON[p]=="-"):
            out = out + ["NC"]
        elif (float(df.Coverage_mutated_sample[p])<float(20)):
            out = out + ["NC"]
        elif (float(df.RUNJ_all_effects_PON[p])>float(lim_runj)):
            out = out +  [True]
        else:
            out = out +  [False]
    df['all_effects_RUNJ_validation_PON'] = out
    return (df)
out = validated_all_effects_RUNJ_PON(out)

def validated_RUNJ_POT(df):
    global lim_runj
    out = []
    for p in range (len(df)):
        if (float(df.Coverage_mutated_sample[p])==float(-1)):
            out = out + ["NC"]
        elif (df.RUNJ_POT[p]=="-"):
            out = out + ["NC"]
        elif (float(df.Coverage_mutated_sample[p])<float(20)):
            out = out + ["NC"]
        elif (float(df.RUNJ_POT[p])>float(lim_runj)):
            out = out +  [True]
        else:
            out = out +  [False]
    df['RUNJ_validation_POT'] = out
    return (df)

out = validated_RUNJ_POT(out)
def validated_all_effects_RUNJ_POT(df):
    global lim_runj
    out = []
    for p in range (len(df)):
        if (float(df.Coverage_mutated_sample[p])==float(-1)):
            out = out + ["NC"]
        elif (float(df.Coverage_mutated_sample[p])==float(-2)):
            out = out + ["NC"]
        elif (float(df.Coverage_mutated_sample[p])<float(20)):
            out = out + ["NC"]
        elif (float(df.RUNJ_all_effects_POT[p])>float(lim_runj)):
            out = out +  [True]
        else:
            out = out +  [False]
    df['all_effects_RUNJ_validation_POT'] = out
    return (df)
out = validated_all_effects_RUNJ_POT(out)

out['Validation_POT'] = out['RUNJ_validation_POT'].map(str) + ['_']*len(out) + out['all_effects_RUNJ_validation_POT'].map(str)

def custom_or(txt):
    a = txt[:txt.index("_")]
    b = txt[txt.index("_")+1:]
    if (a=="True" or b == "True"):
        return ("True")
    elif (a=="False" or b=="False"):
        return ("False")
    else:
        return ("NC")


out['Validation_POT'] = out['Validation_POT'].map(custom_or)
out['Validation_PON'] = out['RUNJ_validation_PON'].map(str) + ['_']*len(out) + out['all_effects_RUNJ_validation_PON'].map(str)
out['Validation_PON'] = out['Validation_PON'].map(custom_or)


out['Validation'] = out['Validation_POT'].map(str) + ['_']*len(out) + out['Validation_PON'].map(str)
out['Validation'] = out['Validation'].map(custom_or)

#out = out.drop_duplicates(subset="unique_id")

def harmonisation_1(txt):
    if (str(txt)=="-1.0"):
        return ("NA")
    else :
        return (txt)

out['Coverage_mutated_sample'] = out['Coverage_mutated_sample'].map(harmonisation_1)
out['RUNJ_PON'] = out['RUNJ_PON'].map(harmonisation_1)
out['RUNJ_all_effects_PON'] = out['RUNJ_all_effects_PON'].map(harmonisation_1)
out['RUNJ_POT'] = out['RUNJ_POT'].map(harmonisation_1)
out['RUNJ_all_effects_POT'] = out['RUNJ_all_effects_POT'].map(harmonisation_1)

loc = out
loc = loc.reset_index()
del loc['index']
###### RUNJ_scored_finalisation

loc = loc.rename(columns={"RT_intron_PON":"RT_intron_mutated_sample"})
loc = loc.rename(columns={"RT_intron_PON_c":"RT_intron_PON"})
del loc['RT_intron_POT']
loc = loc.rename(columns={"RT_intron_POT_c":"RT_intron_POT"})
loc = loc.rename(columns={"ES_PON":"ES_mutated_sample"})
loc = loc.rename(columns={"ES_PON_c":"ES_PON"})
del loc['ES_POT']
loc = loc.rename(columns={"ES_POT_c":"ES_POT"})
loc = loc.rename(columns={"NS_PON":"NS_mutated_sample"})
loc = loc.rename(columns={"NS_PON_c":"NS_PON"})
del loc['NS_POT']
loc = loc.rename(columns={"NS_POT_c":"NS_POT"})
loc = loc.rename(columns={"Coverage_mutated_sample":"Coverage_mutated_sample"})

def all_number(cell):
    try :
        out = int(cell)
    except:
        out = int(0)
    return (out)


loc['RT_intron_mutated_sample'] = loc['RT_intron_mutated_sample'].map(all_number)
loc['RT_intron_PON'] = loc['RT_intron_PON'].map(all_number)
loc['RT_intron_POT'] = loc['RT_intron_POT'].map(all_number)

loc['ES_mutated_sample'] = loc['ES_mutated_sample'].map(all_number)
loc['ES_PON'] = loc['ES_PON'].map(all_number)
loc['ES_POT'] = loc['ES_POT'].map(all_number)

loc['NS_mutated_sample'] = loc['NS_mutated_sample'].map(all_number)
loc['NS_PON'] = loc['NS_PON'].map(all_number)
loc['NS_POT'] = loc['NS_POT'].map(all_number)
loc['Coverage_mutated_sample'] = loc['Coverage_mutated_sample'].map(all_number)
loc["MA"] = loc['NS_mutated_sample'] + loc['ES_mutated_sample'] + loc['RT_intron_mutated_sample']
loc['MN'] = loc['Coverage_mutated_sample'] - loc["MA"]
loc["CA_PON"] = loc['NS_PON'] + loc['NS_PON'] + loc['RT_intron_PON']
loc['CN_PON'] = loc['COV_C_PON'].map(all_number)- loc["CA_PON"]
loc["CA_POT"] = loc['NS_POT'] + loc['NS_POT'] + loc['RT_intron_POT']
loc['CN_POT'] = loc['COV_C_POT'].map(all_number) - loc["CA_POT"]



def Validation_multi(df, lim): 
    out =[]
    for p in range(len(df)):
        if (float(df.RUNJ_PON[p])>lim) or (float(df.RUNJ_all_effects_PON[p])>lim) or (float(df.RUNJ_POT[p])>lim) or (float(df.RUNJ_all_effects_POT[p])>lim):
            out = out + ['True']
            continue
        elif (float(df.RUNJ_PON[p])==-1.0) or (float(df.RUNJ_all_effects_PON[p])==-1.0):
            out = out + ['NC']
        elif ((df.Coverage_mutated_sample[p]<20)  or (df.COV_C_PON[p]<20)):
            out = out + ['NC']
            continue
        elif (df.RUNJ_PON[p]=='NC' or df.MA[p]=='NC' or df.RUNJ_all_effects_PON[p]=='NC' or df.RUNJ_PON[p]!=df.RUNJ_PON[p] or df.RUNJ_all_effects_PON[p]!=df.RUNJ_all_effects_PON[p] ):
            out = out + ['NC']
            continue
        elif(float(df.MA[p])>=2 and (float(df.RUNJ_PON[p])>lim or float(df.RUNJ_all_effects_PON[p])>=lim) and (df.RUNJ_POT[p]==-1 or df.RUNJ_all_effects_POT[p]==-1)):
            out = out + ['True']
            continue
        elif (float(df.MA[p])>=2 and (float(df.RUNJ_PON[p])>lim or float(df.RUNJ_all_effects_PON[p])>=lim) and (float(df.RUNJ_POT[p])>=lim or float(df.RUNJ_all_effects_POT[p])>=lim)):
            out = out + ['True']
            continue
        else:
            out = out + ['False']
            continue
    name = 'Validation_'+str(lim)
    df[name] = out


def harmonize(txt):
    if (txt=="NA"):
        return (-1)
    elif (str(txt)=='-2.0'):
        return(-1)
    else :
        return (txt)

loc['MA'] = loc['MA'].map(harmonize).map(float)
loc['RUNJ_PON'] = loc['RUNJ_PON'].map(harmonize).map(float)
loc['RUNJ_all_effects_PON'] = loc['RUNJ_all_effects_PON'].map(harmonize).map(float)
loc['RUNJ_POT'] = loc['RUNJ_POT'].map(harmonize).map(float)
loc['RUNJ_all_effects_POT'] = loc['RUNJ_all_effects_POT'].map(harmonize).map(float)


Validation_multi(loc, 0.05)
Validation_multi(loc, 0.01)


del loc['Validation']
del loc['Validation_PON']
del loc['Validation_POT']
del loc['RUNJ_validation_PON']
del loc['RUNJ_validation_POT']


loc['Gene'] = loc['SpliceAI_SYMBOL']

loc['SpliceAI - prediction'] = [True]*len(loc)

loc['SpliceAI - maximum DS'] = loc['DS_max']
loc['SpliceAI - event with maximum DS'] = loc['DS_max_kind']
loc['SpliceAI - position of splicing alteration relative to mutation'] = loc['DS_max_pos']


loc['RNA-seq validation of SpliceAI prediction'] = loc['Validation_0.01']


keys = ['sample', 'CHROM', 'POS_hg19', "POS_hg38", 'REF', 'ALT', 'Gene',
'Consequence (VEP)', 'Essential or cryptic', 'SpliceAI - prediction', 
'SpliceAI - maximum DS', 'SpliceAI - event with maximum DS', 'SpliceAI - position of splicing alteration relative to mutation',
'RNA-seq validation of SpliceAI prediction']

loc = loc[keys]

loc.to_excel("Validated_test.xlsx", index=False)









