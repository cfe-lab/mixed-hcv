#!/usr/bin/env python

"""
Calculate the frequency of HCV genotypes and subtypes from a MiSeq
sample (paired FASTQ files) by mapping short reads to a library
of reference genomes.
"""

import os
import subprocess
import re
import bowtie2
import csv
import shutil
from datetime import datetime
import errno

cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token

CODON_DICT = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
              'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'TC-':'S', 'TCN':'S',
              'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
              'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
              'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CT-':'L', 'CTN':'L',
              'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CC-':'P', 'CCN':'P',
              'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
              'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'CG-':'R', 'CGN':'R',
              'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
              'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'AC-':'T', 'ACN':'T',
              'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
              'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
              'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'GT-':'V', 'GTN':'V',
              'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GC-':'A', 'GCN':'A',
              'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
              'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'GG-':'G', 'GGN':'G',
              '---':'-', 'XXX':'X'}




class SamFlag:
    IS_PAIRED =                0x001
    IS_MAPPED_IN_PROPER_PAIR = 0x002
    IS_UNMAPPED =              0x004
    IS_MATE_UNMAPPED =         0x008
    IS_REVERSE =               0x010
    IS_MATE_REVERSE =          0x020
    IS_FIRST =                 0x040
    IS_SECOND =                0x080
    IS_SECONDARY_ALIGNMENT =   0x100
    IS_FAILED =                0x200
    IS_DUPLICATE =             0x400
    IS_CHIMERIC_ALIGNMENT =    0x800

def is_first_read(flag):
    """
    Interpret bitwise flag from SAM field.
    Returns True or False indicating whether the read is the first read in a pair.
    """
    IS_FIRST_SEGMENT = 0x40
    return (int(flag) & IS_FIRST_SEGMENT) != 0

def is_primary(flag):
    """
    :return bool:  Returns whether the current record is primary alignment
    """
    if flag is None:
        raise ValueError("No flag associated with this record")
    return not SamFlag.IS_UNMAPPED & flag and not SamFlag.IS_SECONDARY_ALIGNMENT & flag


def is_chimeric(flag):
    """
    :return bool:  Returns whether the current record is chimeric alignment
    """
    return SamFlag.IS_CHIMERIC_ALIGNMENT & flag

def translate(seq, shift=0):
    """
    Translate nucleotide sequence to AA.  Allows gaps, mixtures, Ns.
    :param str seq: nucleotide sequence
    :param int shift: frameshift.  Will leftpad with gaps before translation.
    :return str:  aa sequence.  Unknown AA are represented with ?
    """
    # Does not handle mixtures
    seq = "-"*shift + seq
    seqlen = len(seq)
    aa_seq = ""
    for codon_site in xrange(0, seqlen, 3):
        if codon_site+3 <= seqlen:
            codon = seq[codon_site:codon_site+3]
            aa_seq += CODON_DICT.get(codon, "?")
    return aa_seq


def count_file_lines(path):
    """ Run the wc command to count lines in a file, as shown here:
    https://gist.github.com/zed/0ac760859e614cd03652
    """
    if path.endswith(".gz"):
        return count_zipped_file_lines(path)
    wc_output = subprocess.check_output(['wc', '-l', path])
    return int(wc_output.split()[0])


def count_zipped_file_lines(path):
    """ Run the wc command to count lines in a gzipped file, as shown here:
    https://gist.github.com/zed/0ac760859e614cd03652
    http://stackoverflow.com/questions/4846891/python-piping-output-between-two-subprocesses
    http://superuser.com/questions/135329/count-lines-in-a-compressed-file
    """
    unzipped_pipe = subprocess.Popen(['gunzip', '-c', path], stdout=subprocess.PIPE)
    wc_process = subprocess.Popen(['wc', '-l'], stdin=unzipped_pipe.stdout, stdout=subprocess.PIPE)
    unzipped_pipe.stdout.close() # enable write error in zcat if wc dies
    wc_out, wc_err = wc_process.communicate()

    return int(wc_out.split()[0])


def do_map(fastq1, fastq2, refpath, bowtie_threads, min_match_len, min_mapq, min_score, bowtie2_version,
           cache=None):
    """
    Process SAM output as it is streamed from bowtie2 to assign
    short reads to HCV genotypes/subtypes.

    """
    # get size of FASTQ
    progress = 0
    nrecords = 0

    # collect SAM output by refname
    counts = {}
    rejects = {'unknown': 0, 'mapq': 0, 'hybrid': 0, 'low score': 0, 'short': 0}
    total_count = 0

    # This is for collecting read qualities.
    mapqs = {}

    # stream STDOUT from bowtie2
    flags = ['--quiet', '--local', '--no-head']
    
    bowtie2_iter = bowtie2.align_paired(bowtie2_version, refpath, fastq1, fastq2, bowtie_threads, flags=flags)

    for line, line2 in bowtie2_iter:
        items = line.split('\t')
        qname, flag1, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = items[:11]
        qname2, flag2, rname2, pos2, mapq2, cigar2, _, _, _, seq2, qual2 = line2.split('\t')[:11]
        assert qname == qname2, 'qnames do not match'
        # TODO: make this loop over bowtie2 output two lines at a time

        if not rname.startswith('HCV') or not rname2.startswith('HCV'):
            rejects['unknown'] += 1
            continue

        # ignore reads whose mate mapped to a different genotype
        geno1 = rname.split('-')[1][0] if '-' in rname else ''
        geno2 = rname2.split('-')[1][0] if '-' in rname2 else ''
        if geno1 != geno2 or cigar == '*' or cigar2 == '*':
            rejects['hybrid'] += 1
            continue

        # filter out reads based on match length - e.g., primers
        for cig in [cigar, cigar2]:
            tokens = cigar_re.findall(cig)
            match_len = sum([int(token.strip('M')) for token in tokens if token.endswith('M')])
            if match_len < min_match_len:
                rejects['short'] += 1
                continue

        if cache:
            # write SAM output to file, but only if the read mapped to HCV
            cache.write(line)
            cache.write(line2)

        # filter out reads with low alignment score
        #optionals = dict([item.split(':i:') for item in items[11:] if ':i:' in item])
        #if 'AS' not in optionals or int(optionals['AS']) < min_score:
        #    rejects['low score'] += 1
        #    continue

        # update counts
        if geno1 not in counts:
            counts.update({geno1: 0})
        counts[geno1] += 1
        total_count += 1


    # output results
    keys = counts.keys()
    keys.sort()

    bowtie2_iter.close()

    return counts, rejects, mapqs


def mixed_hcv(fastq1, fastq2, outpath, refpath, bowtie2_version,
              min_match_len=100, min_mapq=0, min_score=0, n_threads=4, 
              mapq_outfile=None, cache=None):
    """
    Calls do_map and handles writing to output file
    ASSUMES:
    - runname is the name of the parent directory of the fastq file.
    - fastq file follows miseq sample naming conventions:
        <enum>_S<sample num>_L001_R2_001.fastq.gz

    :param refpath: Path to bowtie2 index files
    :param fastq1: Path to R1 FASTQ
    :param fastq2: Path to R2 FASTQ
    :param outpath: Where to write CSV output
    :param n_threads: Number of bowtie2 threads to run
    :param mapq_path: where to write mapping quality information
    :return:
    """
    # check if index files exist
    extensions = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    if not all([os.path.exists(refpath+ext) for ext in extensions]):
        # generate index files
        bowtie2.build(bowtie2_version, refpath)

    # new output file
    #runname,sample,snum,subtype,count,total,perc
    outpath.write('subtype,count,total,perc\n')

    counts, discards, mapqs = do_map(
        fastq1, fastq2, refpath=refpath, bowtie_threads=n_threads, bowtie2_version=bowtie2_version,
        min_match_len=min_match_len, min_mapq=min_mapq, min_score=min_score, cache=cache
    )
    n_discard = sum(discards.values())
    total_count = sum(counts.values()) + n_discard

    #runname, sample, snum, subtype, count, total_count, perc)

    for subtype, count in counts.iteritems():
        perc = 0 if not total_count else count*100/float(total_count)
        outpath.write('%s,%d,%d,%.4g\n' % (subtype, count, total_count, perc))

    # record number of reads that failed to map
    disc_perc = 0 if not total_count else n_discard*100/float(total_count)
    outpath.write(',%d,%d,%.4g\n' % (n_discard, total_count, disc_perc))

    # If mapq_outfile is specified, write the mapping quality information out to it.
    if mapq_outfile:
        mapq_csv_writer = csv.writer(mapq_outfile)
        mapq_csv_writer.writerow(("qname", "flag", "mapq"))
        for qname, flag in mapqs:
            mapq_csv_writer.writerow((qname, flag, mapqs[(qname, flag)]))

