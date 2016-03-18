#!/usr/bin/env python
"""
Kive-ready script for running the random primer HCV version of the
mixed-hcv pipeline.  Derived from batch-mpi.py.  Replaces optional 
and keyword command-line arguments with parameters to be specified 
at top of script.  Since batch execution is handled by the Kive
fleet, this script takes only a single pair of (gzip compressed) 
FASTQ files as its inputs.
"""
import argparse
import helper

BOWTIE2_VERSION = "2.2.1"

BOWTIE_INDEX = 'gb-ref_hg38.fa'  # gb-ref+hg38_v2
BOWTIE_THREADS = 4
MIN_MATCH_LEN = 100
MIN_MAP_QUALITY = 0
MIN_ALIGN_SCORE = 0


def main():
    parser = argparse.ArgumentParser(
        description='Map contents of FASTQ R1 and R2 data sets to references using bowtie2.')

    parser.add_argument('fastq1', help='<input> FASTQ containing forward reads')
    parser.add_argument('fastq2', help='<input> FASTQ containing reverse reads')
    parser.add_argument('sam', type=argparse.FileType('w'), 
                        help='<output> SAM output of bowtie2')
    parser.add_argument('counts', type=argparse.FileType('w'),
                        help='<output> CSV, counts of reads mapped to reference types')

    args = parser.parse_args()
    
    helper.mixed_hcv(
        # inputs
        fastq1=args.fastq1, 
        fastq2=args.fastq2,
        # outputs
        outpath=args.counts,  # CSV output
        cache=args.sam,  # unstructured output
        mapq_outfile=None,
        # settings
        refpath=BOWTIE_INDEX,
        bowtie2_version=BOWTIE2_VERSION,
        min_match_len=MIN_MATCH_LEN,
        min_mapq=MIN_MAP_QUALITY,
        min_score=MIN_ALIGN_SCORE,
        n_threads=BOWTIE_THREADS
    )
    

if __name__ == '__main__':
    main()