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

cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token
gpfx = re.compile('^[-]+')  # length of gap prefix

class HCVDeli():
    """
    Because it cuts HCV reads into delicious slices.
    """
    def __init__(self, x, p, minlen, minq, mins, coords_file='data/gb-ref2.coords'):
        self.targets = {'NS5a': [(0, 250)], 'NS5b': [(100, 350), (600, 850)], 'NS3': [(700, 950)]}
        self.min_overlap = 250

        self.refpath = x
        self.bowtie_threads = p
        self.min_match_len = minlen
        self.min_mapq = minq
        self.min_score = mins

        # load reference gene coordinates
        handle = open(coords_file, 'rU')
        self.coords = {}
        for line in handle:
            rname, gene, left, right = line.strip('\n').split(',')
            if rname not in self.coords:
                self.coords.update({rname: {}})
            self.coords[rname].update({gene: (int(left), int(right))})

    def timestamp(self, msg):
        return '[%s] %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), msg)

    def is_first_read(self, flag):
        """
        Interpret bitwise flag from SAM field.
        Returns True or False indicating whether the read is the first read in a pair.
        """
        IS_FIRST_SEGMENT = 0x40
        return (int(flag) & IS_FIRST_SEGMENT) != 0

    def count_file_lines(self, path):
        """ Run the wc command to count lines in a file, as shown here:
        https://gist.github.com/zed/0ac760859e614cd03652
        """
        wc_output = subprocess.check_output(['wc', '-l', path])
        return int(wc_output.split()[0])

    def apply_cigar (self, cigar, seq, qual):
        """
        Use CIGAR string (Compact Idiosyncratic Gapped Alignment Report) in SAM data
        to apply soft clips, insertions, and deletions to the read sequence.
        Any insertions relative to the sample consensus sequence are discarded to
        enforce a strict pairwise alignment, and returned separately in a
        dict object.
        """
        newseq = ''
        newqual = ''
        insertions = {}
        tokens = cigar_re.findall(cigar)
        if len(tokens) == 0:
            return None, None, None
        # Account for removing soft clipped bases on left
        shift = 0
        if tokens[0].endswith('S'):
            shift = int(tokens[0][:-1])
        left = 0
        for token in tokens:
            length = int(token[:-1])
            # Matching sequence: carry it over
            if token[-1] == 'M':
                newseq += seq[left:(left+length)]
                newqual += qual[left:(left+length)]
                left += length
            # Deletion relative to reference: pad with gaps
            elif token[-1] == 'D':
                newseq += '-'*length
                newqual += ' '*length  # Assign fake placeholder score (Q=-1)
            # Insertion relative to reference: skip it (excise it)
            elif token[-1] == 'I':
                insertions.update({left: (seq[left:(left+length)], qual[left:(left+length)])})
                left += length
                continue
            # Soft clipping leaves the sequence in the SAM - so we should skip it
            elif token[-1] == 'S':
                left += length
                continue
            else:
                print "Unable to handle CIGAR token: {} - quitting".format(token)
                sys.exit()

        return shift, newseq, newqual, insertions

    def merge_pairs (self, seq1, seq2, qual1, qual2, q_cutoff=10, minimum_q_delta=5):
        """
        Combine paired-end reads into a single sequence by managing discordant
        base calls on the basis of quality scores.
        """
        mseq = ''
        # force second read to be longest of the two
        if len(seq1) > len(seq2):
            seq1, seq2 = seq2, seq1
            qual1, qual2 = qual2, qual1

        for i, c2 in enumerate(seq2):
            q2 = ord(qual2[i])-33
            if i < len(seq1):
                c1 = seq1[i]
                q1 = ord(qual1[i])-33
                if c1 == '-' and c2 == '-':
                    mseq += '-'
                    continue
                if c1 == c2:  # Reads agree and at least one has sufficient confidence
                    if q1 > q_cutoff or q2 > q_cutoff:
                        mseq += c1
                    else:
                        mseq += 'N'  # neither base is confident
                else:
                    if abs(q2 - q1) >= minimum_q_delta:
                        if q1 > max(q2, q_cutoff):
                            mseq += c1
                        elif q2 > max(q1, q_cutoff):
                            mseq += c2
                        else:
                            mseq += 'N'
                    else:
                        mseq += 'N'  # cannot resolve between discordant bases
            else:
                # past end of read 1
                if c2 == '-' and q2 == 0:
                    mseq += 'n'  # interval between reads
                    continue
                mseq += c2 if (q2 > q_cutoff) else 'N'
        return mseq

    def align(self, fastq1, fastq2):
        """
        Process SAM output as it is streamed from bowtie2 to assign
        short reads to HCV genotypes/subtypes.
        """
        # get size of FASTQ
        nrecords = self.count_file_lines(fastq1) / 2
        print self.timestamp('total reads %d' % nrecords)

        # stream output from bowtie2
        bowtie_args = ['bowtie2',
                       '--quiet',
                       '-x', self.refpath,
                       '-1', fastq1,
                       '-2', fastq2,
                       '--no-unal',
                       '--no-mixed',  # mates must map to same reference
                       '--local',
                       '-p', str(self.bowtie_threads)]

        p = subprocess.Popen(bowtie_args, stdout=subprocess.PIPE)

        progress = 0
        read_cache = {}
        aligned = {}
        # stream bowtie2 output, gather and merge reads
        with p.stdout:
            for line in p.stdout:
                if line.startswith('@'):
                    # skip header line
                    continue

                progress += 1
                if progress % 5000 == 0:
                    print(self.timestamp('mapped %d of %d' % (progress, nrecords)))

                items = line.split('\t')
                qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = items[:11]
                pos = int(pos)-1  # convert to 0-index

                _, seq1, qual1, inserts = self.apply_cigar(cigar, seq, qual)
                seq2 = '-' * pos + seq1
                qual2 = '!' * pos + qual1

                if qname not in read_cache:
                    read_cache.update({qname: (rname, seq2, qual2)})
                    continue

                # pop mate from cache
                _, seq1, qual1 = read_cache.pop(qname)
                mseq = self.merge_pairs(seq1, seq2, qual1, qual2, q_cutoff=15)
                if mseq.count('N') / float(len(mseq.strip('-'))) > 0.2:
                    continue

                key = (mseq, rname, pos)
                if key not in aligned:
                    aligned.update({key: 0})
                aligned[key] += 1

            if p.returncode:
                raise subprocess.CalledProcessError(p.returncode, bowtie_args)

        return aligned


    def slice(self, aligned):
        """
        Slice segments out of aligned reads if they map to correct gene regions.
        Map position should be relative to the reference genome that it mapped to.
        Insertions will have been filtered out by apply_cigar().
        :param aligned:
        :return:
        """
        counts = {}
        slices = {}
        for (mseq, rname, pos), count in aligned.iteritems():
            coords = self.coords[rname]
            subtype = rname.split('-')[1]  # 'HCV-1a' -> '1a'
            genotype = subtype[0]

            # update counts
            if subtype not in counts:
                counts.update({subtype: 0})
            counts[subtype] += 1

            # check where this read mapped
            read_start = int(pos)
            read_end = read_start + len(mseq.strip('-'))

            coords = self.coords[rname]
            if read_end < coords['Core'][0] or read_start > coords['NS5b'][1]:
                # read falls outside of ORF
                continue

            # did it map to one of the target genes?
            for target_gene, target_coords in self.targets.iteritems():
                if target_gene not in slices:
                    slices.update({target_gene: {}})
                left, right = coords[target_gene]

                for tc in target_coords:
                    if tc not in slices[target_gene]:
                        slices[target_gene].update({tc: []})
                    gene_left, gene_right = tc  # unpack tuple

                    # adjust gene coordinates to genome coordinates
                    genome_left = gene_left + left
                    genome_right = gene_right + left
                    if read_start > genome_left or read_end < genome_right:
                        # full coverage of target not possible
                        continue

                    slice = mseq[genome_left:genome_right]
                    slices[target_gene][tc].append((slice, count))

        return slices, counts


    def run(self, f1, f2, handle, log, runname='', complete=[]):
        """
        Analyze a pair of FASTQ files.
        :param f1: FASTQ R1 input
        :param f2: FASTQ R2 input
        :param handle: open file handle to write output
        :param log:
        :param runname:
        :param complete:
        :return:
        """
        filename = os.path.basename(f1)
        if filename.startswith('Undetermined'):
            return

        sample, snum = filename.split('_')[:2]

        if (runname, sample, snum) in complete:
            # skip previously processed sample
            return

        logfile = open(log, 'a')
        logfile.write('[%s] start processing %s\n' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), filename))
        logfile.close()

        aligned = self.align(f1, f2)
        slices, _ = self.slice(aligned)
        # write out to file
        for gene, subsets in slices.iteritems():
            for coords, subset in subsets.iteritems():
                left, right = coords
                intermed = [(count, seq) for seq, count in subset]
                intermed.sort(reverse=True)
                for rank, (count, seq) in enumerate(intermed):
                    handle.write('%s,%s,%s,%s,%d,%d,%d,%d,%s\n' % (runname, sample, snum, gene,
                                                                   left, right, rank+1, count, seq))
                    handle.flush()  # make sure we write intact lines

        logfile = open(log, 'a')
        logfile.write('[%s] end processing %s\n' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), filename))
        logfile.close()


    def batch(self, paths, output, log):
        """
        Set up batch run on multiple sets of FASTQ files in a path or list of paths.
        :param paths:
        :param output:
        :param log:
        :param bowtie_args:
        :return:
        """

        # enable script to restart from interrupted run
        complete = {}
        if os.path.exists(output):
            handle = open(output, 'rU')
            _ = handle.next()  # skip header
            for line in handle:
                runname, sample, snum = line.strip('\n').split(',')[:3]
                complete.update({(runname, sample, snum): 1})
            handle.close()
            # re-open file for appending
            handle = open(output, 'a')
        else:
            # new output file
            handle = open(output, 'w')
            handle.write('runname,sample,snum,gene,start,end,rank,count,seq\n')

        for path in paths:
            # check that the inputs exist
            if not os.path.exists(path):
                print('No folder/file found at %s', path)
                sys.exit(1)

            if not path.endswith('/'):
                path += '/'

            files = glob(path + '*_R1_001.fastq')
            if len(files) == 0:
                print 'ERROR: No FASTQ R1 files found at', path
                sys.exit()

            runname = path.split('/')[4]
            for f1 in files:
                f2 = f1.replace('_R1_', '_R2_')
                self.run(f1, f2, handle, log, runname=runname, complete=complete)

        handle.close()


def main():
    # check that we have access to bowtie2
    try:
        subprocess.check_output(['bowtie2', '-h'])
    except OSError:
        raise RuntimeError('bowtie2 not found; check if it is installed and in $PATH\n')

    parser = argparse.ArgumentParser(
        description='Map contents of FASTQ R1 and R2 data sets to references using bowtie2.')

    # TODO: add subparser
    parser.add_argument('-R1', help='<input> FASTQ R1')
    parser.add_argument('-R2', help='<input> FASTQ R2')
    parser.add_argument('-path', help='<input> folder containing FASTQ files')
    parser.add_argument('-pathlist', help='<input> file containing list of folders')

    parser.add_argument('output', help='<output CSV> file to write results')

    parser.add_argument('-log', help='path to logfile', default='hcvdeli.log')
    parser.add_argument('-x', help='path to bowtie2 index (.bt2)', default='data/gb-ref2')
    parser.add_argument('-p', type=int, help='number of bowtie2 threads', default=6)
    parser.add_argument('-minlen', type=int, help='minimum match length (CIGAR M)', default=100)
    parser.add_argument('-minq', type=int, help='minimum mapping quality (MAPQ)', default=0)
    parser.add_argument('-mins', type=int, help='minimum alignment score', default=0)

    args = parser.parse_args()
    if not ((args.R1 and args.R2) or args.path or args.pathlist):
        parser.error('Must set one of the following: {-R1 and -R2, -path, -pathlist}.')
    deli = HCVDeli(x=args.x, p=args.p, minlen=args.minlen, minq=args.minq, mins=args.mins)

    paths = []
    if args.pathlist:
        handle = open(args.pathlist, 'rU')
        for line in handle:
            paths.append(line.strip('\n'))
        handle.close()
        deli.batch(paths=paths, output=args.output, log=args.log)
    elif args.path:
        paths.append(args.path)
        deli.batch(paths=paths, output=args.output, log=args.log)
    elif args.R1:
        # no continuation of run for single file mode!
        handle = open(args.output, 'w')
        handle.write('runname,sample,snum,gene,start,end,rank,count,seq\n')
        deli.run(f1=args.R1, f2=args.R2, handle=handle, log=args.log)
        handle.close()
    else:
        # this should never happen!
        parser.print_help()
        sys.exit()


if __name__ == '__main__':
    main()
