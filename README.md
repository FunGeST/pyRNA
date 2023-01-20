# pyRNA

pyRNA library allows to mine RNA-seq data in Python. Technically, it converts pysam files into pandas dataframe, allowing easy downstream analyses of RNA sequences, for example to detect sequence variants. PyRNA was used to validate mutations predicted to alter splicing using paired RNA-seq data, as detailed in this article: https://doi.org/10.1101/2022.10.14.512264. 

It can either : 

1) Output a dataframe of a locus from a bam (A='reads_to_dataframe')
2) Output only the read of that said locus (A='write_reads')
3) Search for a splicing anomaly at any given position. (A='splicing_analysis') It will compute and print the following scores : 

The usage of the novel junction created by the abnormal splicing event for AG (Acceptor Gain) / AL (Acceptor Loss) / DG (Donor Gain) / DL (Donor Loss). 

Any score is computed as such : 

DS = (#Abnormal reads corroborrating the event) / (#Total of reads)

Principle : 

python pyRNA.py --pathM PATH_BAM -S SAMPLE_NAME -A ACTION -C CHROM -P POSITION -L WINDOW -E EVENT -G GENE -R REFERENCE


Examples  :

python pyRNA.py --pathM ./example -S test.bam -A write_reads -C chr11 -P 119090189 -L 100 -R grch38

python pyRNA.py --pathM ./example -S test.bam -A reads_to_dataframe -C chr11 -P 119090189 -L 100 -R grch38

python pyRNA.py --pathM ./example -S test.bam -A splicing_analysis -E AL -G CBL -C chr11 -P 119090189 -L 100 -R grch38



# Validation demo

You can run a validation on a few examples with sliced paired RNA-seq data. 

To run it : python Validation_full_script.py test.xlsx

NB : I had to tweak a bit my pipeline to make it work with custom sliced control, but I have commented what to change to use it with your own data. 


