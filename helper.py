#!/usr/bin/env python

"""
Calculate the frequency of HCV genotypes and subtypes from a MiSeq
sample (paired FASTQ files) by mapping short reads to a library
of reference genomes.
"""

from errno import ENOENT
import os
import re
import subprocess

import bowtie2

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

    # collect SAM output by refname
    counts = {}
    rejects = {'unknown': 0, 'mapq': 0, 'hybrid': 0, 'low score': 0, 'short': 0}

    # stream STDOUT from bowtie2
    flags = ['--quiet', '--local', '--no-head']
    
    bowtie2_iter = bowtie2.align_paired(bowtie2_version, refpath, fastq1, fastq2, bowtie_threads, flags=flags)

    for line, line2 in bowtie2_iter:
        qname, flag1, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.split('\t')[:11]
        qname2, flag2, rname2, pos2, mapq2, cigar2, _, _, _, seq2, qual2 = line2.split('\t')[:11]
        assert qname == qname2, 'qnames do not match'

        # ignore pairs if one read fails to map to HCV
        if not rname.startswith('HCV') or not rname2.startswith('HCV'):
            rejects['unknown'] += 1
            continue

        if cigar == '*' or cigar2 == '*':
            rejects['unknown'] += 1
            continue

        # ignore reads whose mate mapped to a different subtype
        subtype1 = rname.split('-')[1] if '-' in rname else ''
        subtype2 = rname2.split('-')[1] if '-' in rname2 else ''
        if subtype1 != subtype2:
            # genotypes do not match
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
        if subtype1 not in counts:
            counts.update({subtype1: 0})
        counts[subtype1] += 1

    bowtie2_iter.close()

    return counts, rejects


def _create_link(fastq, link_root):
    """ Create a symbolic link to a FASTQ file with the right extension.

    bowtie2 decides whether to unzip a FASTQ file based on its file name
    extension. Kive will copy raw FASTQ and compressed FASTQ with the .RAW
    extension, so we need to figure out which we've received. Then create
    a symbolic link with the correct extension to pass to bowtie2.
    :param str fastq: path to a FASTQ file that is text or gzipped
    :param str link_root: root of the path where the symbolic link will be
        created. It will have an extension added to it.
    :return: the path to the new symbolic link
    """
    with open(fastq, 'rb') as f:
        # Check the header for the Gzip magic number.
        # See gzip.GzipFile._read_zip_header().
        magic = f.read(2)
    extension = '.fastq.gz' if magic == '\037\213' else '.fastq'
    link_path = link_root + extension
    try:
        os.remove(link_path)
    except OSError as ex:
        if ex.errno == ENOENT:
            pass  # The file didn't exist, nothing to remove.
        else:
            raise
    os.symlink(fastq, link_path)
    return link_path


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

    fastq1_link = _create_link(fastq1, 'r1')
    fastq2_link = _create_link(fastq2, 'r2')

    # new output file
    # runname,sample,snum,subtype,count,total,perc
    outpath.write('subtype,count,total,perc\n')

    counts, discards = do_map(fastq1_link,
                              fastq2_link,
                              refpath=refpath,
                              bowtie_threads=n_threads,
                              bowtie2_version=bowtie2_version,
                              min_match_len=min_match_len,
                              min_mapq=min_mapq,
                              min_score=min_score,
                              cache=cache)
    n_discard = sum(discards.values())
    total_count = sum(counts.values()) + n_discard

    # runname, sample, snum, subtype, count, total_count, perc)

    for subtype, count in counts.iteritems():
        perc = 0 if not total_count else count*100/float(total_count)
        outpath.write('%s,%d,%d,%.4g\n' % (subtype, count, total_count, perc))

    # record number of reads that failed to map
    disc_perc = 0 if not total_count else n_discard*100/float(total_count)
    outpath.write(',%d,%d,%.4g\n' % (n_discard, total_count, disc_perc))
