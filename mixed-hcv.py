"""
Calculate the frequency of HCV genotypes and subtypes from a MiSeq
sample (paired FASTQ files) by mapping short reads to a library
of reference genomes.
"""
import os
import argparse
import subprocess
import re

cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token

# default settings
refpath = 'data/HCV_REF_2012_genome'
bowtie_threads = 2
min_match_len = 50
min_mapq = 10


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


def map(fastq1, fastq2, refpath, bowtie_threads, min_match_len, min_mapq):
    """
    Process SAM output as it is streamed from bowtie2 to assign
    short reads to HCV genotypes/subtypes.
    """
    # get size of FASTQ
    nrecords = count_file_lines(args.fastq1) / 4
    
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
            if progress % 1000 == 0:
                print progress, nrecords, total_count
                
            subtype = rname.split('-')[-1]  # 'HCV-1a' -> '1a'
            if subtype not in counts:
                counts.update({subtype: 0})
            
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
                print rname, rnext
                n_hybrid += 1
                continue
                
            # filter out primers based on match length
            tokens = cigar_re.findall(cigar)
            match_len = sum([int(token.strip('M')) for token in tokens if token.endswith('M')])
            if match_len < min_match_len:
                n_short += 1
                continue
            
            counts[subtype] += 1
            total_count += 1
            
        if p.returncode:
            raise subprocess.CalledProcessError(p.returncode, bowtie_args)
    
    # output results
    keys = counts.keys()
    keys.sort()
    
    for subtype in keys:
        print subtype, counts[subtype]

    print 'too short', n_short
    print 'low mapq', n_mapq
    print 'mates mapped to different subtype', n_hybrid


def main():
    parser = argparse.ArgumentParser(
        description='Map contents of FASTQ R1 and R2 data sets to references using bowtie2.')
    
    parser.add_argument('fastq1', help='<input> FASTQ containing forward reads')
    parser.add_argument('fastq2', help='<input> FASTQ containing reverse reads')
    parser.add_argument('-x', help='path to bowtie2 index (.bt2)', default=refpath)
    parser.add_argument('-p', type=int, help='number of bowtie2 threads', default=bowtie_threads)
    parser.add_argument('-minlen', type=int, help='minimum match length (CIGAR M)', default=min_match_len)
    parser.add_argument('-minq', type=int, help='minimum mapping quality (MAPQ)', default=min_mapq)
    
    args = parser.parse_args()
    
    # check that we have access to bowtie2
    try:
        subprocess.check_output(['bowtie2', '-h'])
    except OSError:
        raise RuntimeError('bowtie2 not found; check if it is installed and in $PATH\n')
        
            # check that the inputs exist
    if not os.path.exists(args.fastq1):
        print('No FASTQ found at %s', args.fastq1)
        sys.exit(1)

    if not os.path.exists(args.fastq2):
        print('No FASTQ found at %s', args.fastq2)
        sys.exit(1)
    
    map(args.fastq1, args.fastq2, )



if __name__ == '__main__':
    main()
