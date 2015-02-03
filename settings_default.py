"""
Default settings for mixed-hcv
Make a copy of this file, modify as needed and save as 'settings.py'
"""

# location of bowtie2 executable - if it is PATH then leave as is
path_to_bowtie = 'bowtie2'

# location of HCV genotype reference sequence set - FASTA
refpath = 'data/HCV_REF_2012_genome'

# number of threads to run bowtie2
bowtie_threads = 6

# filter out reads with match length below this cutoff - to remove primers
min_match_len = 50

# minimum map quality
min_mapq = 10
