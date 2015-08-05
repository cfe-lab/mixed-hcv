"""
Calculate the frequency of HCV genotypes and subtypes from a MiSeq
sample (paired FASTQ files) by mapping short reads to a library
of reference genomes.
"""
import sys
import os
import argparse
import subprocess
import re
from glob import glob
from datetime import datetime
import helper

from mpi4py import MPI
my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token
BOWTIE2_VERSION = "2.2.1"



def main():
    parser = argparse.ArgumentParser(
        description='Map contents of FASTQ R1 and R2 data sets to references using bowtie2.')
    
    parser.add_argument('-path', help='<input> folder containing gzipped FASTQ files.  Assumes that parent dir of fastq.gz files is the runname')
    parser.add_argument('-pathlist', help='<input> file containing list of folders')
    parser.add_argument('-output', help='<output CSV> file to write results')
    parser.add_argument('-log', help='path to logfile', default='mixed-hcv.batch.log')
    parser.add_argument('-x', help='path to bowtie2 index (.bt2)', default='data/gb-ref')
    parser.add_argument('-p', type=int, help='number of bowtie2 threads', default=6)
    parser.add_argument('-minlen', type=int, help='minimum match length (CIGAR M)', default=100)
    parser.add_argument('-minq', type=int, help='minimum mapping quality (MAPQ)', default=0)
    parser.add_argument('-mins', type=int, help='minimum alignment score', default=0)

    args = parser.parse_args()
    
    paths = []
    if args.pathlist:
        handle = open(args.pathlist, 'rU')
        for line in handle:
            paths.append(line.strip('\n'))
        handle.close()
    elif args.path:
        paths.append(args.path)
    else:
        parser.print_help()
        sys.exit()
    
    # check that we have access to bowtie2
    try:
        subprocess.check_output(['bowtie2', '-h'])
    except OSError:
        raise RuntimeError('bowtie2 not found; check if it is installed and in $PATH\n')
    
    
    complete = {}
    # new output file

    for pindex, path in enumerate(paths):
        #if pindex % nprocs != my_rank:
        #    continue
        
        # check that the inputs exist
        if not os.path.exists(path):
            print('No folder found at %s', path)
            sys.exit(1)

        if not path.endswith('/'):
            path += '/'

        # get a list of all FASTQ R1 files in this folder
        files = glob(path + '*_R1_001.fastq.gz')
        if len(files) == 0:
            print 'ERROR: No FASTQ R1 files found at', path
            sys.exit()


        runname = os.path.basename(os.path.dirname(os.path.abspath(files[0])))

        for i, f1 in enumerate(files):
            # distribute across nodes
            if i % nprocs != my_rank:
                continue

            f2 = f1.replace('_R1_', '_R2_')
            filename = os.path.basename(f1)
            if filename.startswith('Undetermined'):
                continue

            sample, snum = filename.split('_')[:2]

            if (runname, sample, snum) in complete:
                # skip previously processed sample
                continue


            # 57368-2-HLA-B_S2_L001_I1_001.fastq.gz
            fastq_output_csv = args.output + "." + os.path.basename(f1).split("_L001_R1")[0]
            logfilename = args.log+'.'+str(my_rank)
            with open(fastq_output_csv, 'w') as fh_out_csv, open(logfilename, 'a') as fh_log:
                fh_log.write('[%s] start processing %s\n' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), filename))

                helper.mixed_hcv(fastq1=f1, fastq2=f2, outpath=fh_out_csv, refpath=args.x, bowtie2_version=BOWTIE2_VERSION,
                                 min_match_len=args.minlen, min_mapq=args.minq, min_score=args.mins,
                                 n_threads=args.p, is_show_progress=True)

                fh_log.write('[%s] end processing %s\n' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), filename))



if __name__ == '__main__':
    main()
