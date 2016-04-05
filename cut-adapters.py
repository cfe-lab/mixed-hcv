#!/usr/bin/env python

import subprocess
import argparse

def cutadapt(fastq1, fastq2, out1, out2,
             version='1.9.1',
             adapt1='AATGATACGGCGACCACCGAGATCTACACGCCTCCCTCGCGCCATCAG',
             adapt2='CAAGCAGAAGACGGCATACGAGATNNNNNNNNCGGTCTGCCTTGCCAGCCCGCTCAG',
             times='2'):
    """
    :param fastq1: FASTQ.gz R1 input file
    :param fastq2: FASTQ.gz R2 input file
    :param out1: FASTQ.gz R1 output file
    :param out2: FASTQ.gz R2 output file
    :param adapt1: sequence of adapter ligated to 3' end
    :param adapt2: sequence of adapter ligated to 5' end
    :param times: remove up to N adapters from each read
    """

    # check that cutadapt exists on this system
    try:
        stdout = subprocess.check_output(['cutadapt', '--version'])
    except OSError:
        raise RuntimeError('cutadapt not found; check if it is installed and in $PATH\n')

    # check that the installed version is the expected version
    local_version = stdout.split('\n')[0]
    assert version == local_version, 'cutadapt version incompatibility %s != %s' % (version, local_version)

    subprocess.check_call(['cutadapt', '-n', times, '-a', adapt1, '-A', adapt2,
                           '-o', out1, '-p', out2, fastq1, fastq2])


def main():
    parser = argparse.ArgumentParser(description='Python wrapper for cutadapt')

    parser.add_argument('fastq1', help='<input> FASTQ.gz R1 file')
    parser.add_argument('fastq2', help='<input> FASTQ.gz R2 file')
    parser.add_argument('out1', help='<output> FASTQ.gz R1 file')
    parser.add_argument('out2', help='<output> FASTQ.gz R2 file')
    parser.add_argument('-adapt1', help="Adapter sequence ligated to 3' end", default='AATGATACGGCGACCACCGAGATCTACACGCCTCCCTCGCGCCATCAG')
    parser.add_argument('-adapt2', help="Adapter sequence ligated to 5' end", default='CAAGCAGAAGACGGCATACGAGATNNNNNNNNCGGTCTGCCTTGCCAGCCCGCTCAG')

    args = parser.parse_args()
    cutadapt(args.fastq1, args.fastq2, args.out1, args.out2)

if __name__ == '__main__':
    main()
