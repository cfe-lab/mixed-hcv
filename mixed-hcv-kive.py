#!/usr/bin/env python

"""
Calculate the frequency of HCV genotypes and subtypes from a MiSeq
sample (paired FASTQ files) by mapping short reads to a library
of reference genomes.

Refactored for Kive, take positional arguments.
"""

import sys
import os
import argparse
import subprocess
import re
import bowtie2

bowtie2_version = '2.2.3'
refpath='gb-ref.fa'  # CodeResource will get written to pwd
cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token


def is_first_read(flag):
    """
    Interpret bitwise flag from SAM field.
    Returns True or False indicating whether the read is the first read in a pair.
    """
    IS_FIRST_SEGMENT = 0x40
    return (int(flag) & IS_FIRST_SEGMENT) != 0

def count_file_lines(path):
    """ Run the wc command to count lines in a file, as shown here:
    https://gist.github.com/zed/0ac760859e614cd03652
    """
    wc_output = subprocess.check_output(['wc', '-l', path])
    return int(wc_output.split()[0])


def do_map(fastq1, fastq2, bowtie_threads, min_match_len=100, min_mapq=0, min_score=0):
    """
    Process SAM output as it is streamed from bowtie2 to assign
    short reads to HCV genotypes/subtypes.
    """
    counts = {}
    rejects = {'mapq': 0, 'hybrid': 0, 'low score': 0, 'short': 0}
    total_count = 0

    # stream STDOUT from bowtie2
    bowtie2_iter = bowtie2.align_paired(bowtie2_version, refpath, fastq1, fastq2, bowtie_threads)

    for line in bowtie2_iter:
        if line.startswith('@'):
            # skip header line
            continue

        items = line.split('\t')
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = items[:11]

        subtype = rname.split('-')[1]  # 'HCV-1a' -> '1a'
        genotype = subtype[0]

        # ignore second reads
        if not is_first_read(flag):
            continue

        # discard reads with low map quality
        q = int(mapq)
        if q < min_mapq:
            rejects['mapq'] += 1
            continue

        # ignore reads whose mate mapped to a different genotype
        if rnext != '=' and genotype != rnext.split('-')[1][0]:
            rejects['hybrid'] += 1
            continue

        # filter out reads based on match length - e.g., primers
        tokens = cigar_re.findall(cigar)
        match_len = sum([int(token.strip('M')) for token in tokens if token.endswith('M')])
        if match_len < min_match_len:
            rejects['short'] += 1
            continue

        # filter out reads with low alignment score
        optionals = dict([item.split(':i:') for item in items[11:] if ':i:' in item])
        if 'AS' not in optionals or int(optionals['AS']) < min_score:
            rejects['low score'] += 1
            continue

        # update counts
        if subtype not in counts:
            counts.update({subtype: 0})
        counts[subtype] += 1
        total_count += 1


    # output results
    keys = counts.keys()
    keys.sort()

    return counts, rejects


def mixed_hcv(fastq1, fastq2, outpath, n_threads=4):
    """
    Calls do_map and handles writing to output file
    :param refpath: Path to bowtie2 index files
    :param fastq1: Path to R1 FASTQ
    :param fastq2: Path to R2 FASTQ
    :param outpath: Where to write CSV output
    :param n_threads: Number of bowtie2 threads to run
    :return:
    """
    # generate index files
    bowtie2.build(bowtie2_version, refpath)

    # new output file
    outpath.write('subtype,count,total\n')

    counts, discards = do_map(fastq1, fastq2, n_threads)
    n_discard = sum(discards.values())
    total_count = sum(counts.values()) + n_discard

    for subtype, count in counts.iteritems():
        outpath.write('%s,%d,%d\n' % (subtype, count, total_count))

    # record number of reads that failed to map
    outpath.write(',%d,%d\n' % (n_discard, total_count))


def main():
    """
    For shell execution, Kive style (positional arguments)
    :return:
    """
    parser = argparse.ArgumentParser(description='Map contents of paired Illumina MiSeq FASTQ R1 and R2 data sets to '
                                                 'HCV genome references using bowtie2.')

    parser.add_argument('fastq1', help='<input> FASTQ containing forward reads')
    parser.add_argument('fastq2', help='<input> FASTQ containing reverse reads')
    parser.add_argument('outpath',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing counts of reads mapped to HCV genotypes and subtypes.')

    args = parser.parse_args()
    mixed_hcv(fastq1=args.fastq1, fastq2=args.fastq2, outpath=args.outpath)

if __name__ == '__main__':
    main()
