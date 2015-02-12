"""
Calculate the frequency of HCV genotypes and subtypes from a MiSeq
sample (paired FASTQ files) by mapping short reads to a library
of reference genomes.

For execution on cluster.
"""

import sys
import os
import subprocess
import re
from glob import glob
from datetime import datetime
from mpi4py import MPI
from settings import *

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

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


def do_map(fastq1, fastq2, refpath, bowtie_threads, min_match_len, min_mapq):
    """
    Process SAM output as it is streamed from bowtie2 to assign
    short reads to HCV genotypes/subtypes.
    """
    # get size of FASTQ
    nrecords = count_file_lines(fastq1) / 2

    # stream output from bowtie2
    bowtie_args = ['bowtie2',
                   '--quiet',
                   '-x', refpath,
                   '-1', fastq1,
                   '-2', fastq2,
                   '--no-unal',
                   '--local',
                   '-p', str(bowtie_threads)]

    p = subprocess.Popen(bowtie_args, stdout=subprocess.PIPE)

    # collect SAM output by refname
    counts = {}
    n_short = 0
    n_hybrid = 0
    n_mapq = 0
    progress = 0
    total_count = 0

    with p.stdout:
        for line in p.stdout:
            if line.startswith('@'):
                # skip header line
                continue

            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.split('\t')[:11]
            progress += 1
            if progress % 5000 == 0:
                print('[%s] (%d/%d) mapped %d' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                    progress, nrecords, 2*total_count))

            subtype = rname.split('-')[1]  # 'HCV-1a' -> '1a'

            # ignore second reads
            if not is_first_read(flag):
                continue

            # discard reads with low map quality
            q = int(mapq)
            if q < min_mapq:
                n_mapq += 1
                continue

            # ignore reads whose mate mapped to a different reference
            if rnext != '=':
                #print rname, rnext
                n_hybrid += 1
                continue

            # filter out primers based on match length
            tokens = cigar_re.findall(cigar)
            match_len = sum([int(token.strip('M')) for token in tokens if token.endswith('M')])
            if match_len < min_match_len:
                n_short += 1
                continue

            # update counts
            if subtype not in counts:
                counts.update({subtype: 0})
            counts[subtype] += 1
            total_count += 1

        if p.returncode:
            raise subprocess.CalledProcessError(p.returncode, bowtie_args)

    # output results
    keys = counts.keys()
    keys.sort()

    return counts, n_short, n_mapq, n_hybrid


################
##  MAIN LOOP
###############


runs = [#'141001_M01841_0078_000000000-AB843', '141020_M01841_0081_000000000-A8CD8',
        '141029_M01841_0086_000000000-A8EJB', '141209_M01841_0092_000000000-AC7T1',
        '150112_M01841_0095_000000000-AA4L1', '150114_M01841_0096_000000000-AC27B',
        '150116_M01841_0097_000000000-AC23N', '150119_M01841_0098_000000000-AA8AV',
        '150121_M01841_0099_000000000-AA2YY', '150123_M01841_0100_000000000-AA4DA',
        '150127_M01841_0101_000000000-AA4GE', '150203_M01841_0103_000000000-AA8B3',
        '140910_M01841_0076_000000000-A8EER']


src = '/data/miseq/'

handle = open('merck-output.csv', 'w')

for run in runs[::-1]:
    files = glob(src+run+'/*_R1_001.fastq')
    for i, f1 in enumerate(files):
        if i % nprocs != my_rank:
            continue

        f2 = f1.replace('_R1_', '_R2_')
        filename = os.path.basename(f1)
        if filename.startswith('Undetermined'):
            continue

        print('[%s] node %d start processing %s\n' % (my_rank,
                                                      datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                      filename))

        sample, snum = filename.split('_')[:2]

        counts, n_short, n_mapq, n_hybrid = do_map(f1, f2, refpath, bowtie_threads, min_match_len, min_mapq)
        n_discard = n_short + n_mapq + n_hybrid
        total_count = sum(counts.values()) + n_discard

        for subtype, count in counts.iteritems():
            handle.write('%s,%s,%s,%s,%d,%d\n' % (runname, sample, snum, subtype, count, total_count))

        # record number of reads that failed to map
        handle.write('%s,%s,%s,,%d,%d\n' % (runname, sample, snum, n_discard, total_count))
        handle.flush()  # clear write buffer

        print('[%s] node %d end processing %s\n' % (my_rank,
                                                    datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                    filename))

handle.close()

