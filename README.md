# pyRNA

pyRNA library allows to mine RNA-seq data in Python. Technically, it converts pysam files into pandas dataframe, allowing easy downstream analyses of RNA sequences, for example to detect sequence variants. pyRNA takes as input an RNA-seq bam file, a gene name and a genomic position of interest (e.g. location of a mutation predicted to impact splicing). It can return data frames representing sequencing reads overlapping the position of interest, and compute the number and proportion of abnormal reads considering the expected structure of the transcript. PyRNA was used to validate mutations predicted to alter splicing using paired RNA-seq data, as detailed in this article: https://doi.org/10.1101/2022.10.14.512264. 



## Usage
> python pyRNA.py --pathM PATH_TO_BAM -S BAM_NAME -A ACTION -C CHROM -P POSITION -L WINDOW -E EVENT -G GENE -R REFERENCE

`--pathM` directory containing the bam file

`-S` name of the bam file

`-A` action required ('reads_to_dataframe', 'write_reads' or 'splicing_analysis')

`-C` chromosome

`-P` position to look at in the bam file

`-L` size of the region around the mutation that will be considered. For example, if '-L 100', all reads overlapping the genomic window from P-100 to P+100 will be included.

`-E` expected splicing alteration (for 'splicing_analysis'). One of: 'AG' (Acceptor Gain), 'AL' (Acceptor Loss), 'DG' (Donor Gain) or 'DL' (Donor Loss)

`-G` gene (required for 'splicing_analysis'). Gene predicted to be altered by the splicing alteration.

`-R` reference genome. One of 'grch37' or 'grch38'. Defines the gtf file that will be used for reference transcript structures.


## Outputs

Depending on the selected action (-A), 3 different outputs will be produced:

`-A reads_to_dataframe` will generate a dataframe comprising an extraction of the bam file restricted to reads overlapping the specified region. The output file will include the reads but also other fields from the bam file like CIGAR string, MAPQ score etc.

`-A write_reads` will generate a dataframe comprising only the sequences of reads overlapping the specified region. This file can be easily explored visually to spot abnormal reads.

`-A splicing_analysis` will compute and print the number of abnormal reads supporting the predicted splicing alteration (specified in -E), the total number of reads at the locus and the proportion of abnormal reads.


## Examples
The **example** folder contains a toy bam file to test the tool. Below are a few example command lines.

* Extract a portion of a bam file corresponding to a region spanning 100 bases before and after chr11:119090189:
> python pyRNA.py --pathM ./example -S test.bam -A reads_to_dataframe -C chr11 -P 119090189 -L 100 -R grch38

* Extract reads overlapping a region spanning 100 bases before and after chr11:119090189 and output them in a data frame:
> python pyRNA.py --pathM ./example -S test.bam -A write_reads -C chr11 -P 119090189 -L 100 -R grch38

* Compute the number and proportion of reads supporting a splice acceptor loss at position chr11:119090189 in HMBS gene:
> python pyRNA.py --pathM ./example -S test.bam -A splicing_analysis -E AL -G HMBS -C chr11 -P 119090189 -L 100 -R grch38




# Validation demo

The **Validation_Demo** folder containing a demo script to validate a few mutations predicted to impact splicing on matched RNA-seq data. The folder also contains sliced RNA-seq bam files containing only reads around the positions of interest that are necessary to run the script. Use the following command to run the demo:
> python Validation_full_script.py test.xlsx



