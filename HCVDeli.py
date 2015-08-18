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
import bowtie2
import settings
import helper
import csv

cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token
gpfx = re.compile('^[-]+')  # length of gap prefix

class HCVDeli():
    """
    Because it cuts HCV reads into delicious slices.
    :param str coords_file:  filepath to comma-separated file with cols:  subtype ref name, gene, 0based start nuc coord wrt subtype gene, 0based end nuc coord+1 wrt subtype gene
    :param str targets_file:  filepath to comma-separated target regions file with cols:  gene, 0based nuc start coord wrt H77 gene, 0based nuc end coord+1 wrt H77 gene
            colnames:  Gene,Startnuc_0based,AfterEndNuc_0based
    :param float min_target_width:  fraction of target region width that must be covered by read in order to be considered a hit.  Value in [0.0, 1.0]
    """
    def __init__(self, x, p, minlen, minq, mins, min_target_width, coords_file='data/gb-ref2.coords', targets_file=None, cache_path=None):

        if not targets_file:
            # H77 nuc gene coordinates, 0 based.  Start pos, end pos + 1
            # These targets aren't in frames (multiples of 3).  Instead, they are based on 250bp sections that determine genotyping accuracy as per
            # http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0122082#pone-0122082-t002
            self.targets = {'NS5a': [(0, 250)], 'NS5b': [(100, 350), (600, 850)], 'NS3': [(700, 950)]}

        else:
            self.targets = dict()
            with open(targets_file, 'rU') as fh_in_tgt:
                csvreader = csv.DictReader(fh_in_tgt)
                for line in csvreader:
                    gene = line["Gene"]
                    start_0based = line["StartNuc_0based"]
                    after_end_0based = line["AfterEndNuc_0based"]
                    if gene not in self.targets:
                        self.targets.update({gene:[]})
                    self.targets[gene].append((int(start_0based), int(after_end_0based)))

        self.min_overlap = 250

        self.p = None  # this will be assigned the bowtie2 process
        self.refpath = x
        self.bowtie_threads = p
        self.min_match_len = minlen
        self.min_mapq = minq
        self.min_score = mins
        self.min_target_width = min_target_width
        self.cache_path = cache_path

        # load reference gene coordinates
        handle = open(coords_file, 'rU')
        self.coords = {}
        for line in handle:
            rname, gene, left, right = line.strip('\n').split(',')
            if rname not in self.coords:
                self.coords.update({rname: {}})
            self.coords[rname].update({gene: (int(left), int(right))})
        handle.close()

    # def __del__(self):
    #     self.p.terminate()  # close down bowtie2, esp. in case of Python exiting

    def timestamp(self, msg):
        return '[%s] %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), msg)



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
        :param str seq1:  mate1 sequence, left-padded wrt reference
        :param str seq2:  mate2 sequence, left-padded wrt reference
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


    def align(self, fastq1, fastq2, is_show_progress=False, cache=None):
        """
        Process SAM output as it is streamed from bowtie2 to assign
        short reads to HCV genotypes/subtypes.
        """
        # get size of FASTQ
        nrecords = 0
        if is_show_progress:
            nrecords = helper.count_file_lines(fastq1) / 2
            print self.timestamp('total reads %d' % nrecords)

        flags = ['--quiet', 
                 '--local',
                 '--no-unal',  # ignore any reads in which both mates don't map to ref
                 '--no-mixed']  # mates must map to same reference])

        # Instantiate a cache object
        cache_file = None
        close_sam = False

        # Check to see if we have a cache, and use it if we do
        if cache is None or not cache.check_sam(fastq1, fastq2, flags):
            bowtie2_iter = bowtie2.align_paired(settings.bowtie2_version, self.refpath, fastq1, fastq2, self.bowtie_threads, flags=flags)

            if cache is not None:
                cache_file = cache.open_sam_cache(fastq1, fastq2, flags, bowtie2_iter)
        else:
            bowtie2_iter = cache.get_sam(fastq1, fastq2, flags)
            close_sam = True


        # We don't care about reads where only single mate maps
        # We don't care about reads where mates map to different reference
        bowtie2_iter = bowtie2.align_paired(version=settings.bowtie2_version,
                                            refpath=self.refpath,
                                            fastq1=fastq1,
                                            fastq2=fastq2, nthreads=self.bowtie_threads,
                                            flags=flags)



        progress = 0
        read_cache = {}
        aligned = {}
        total_indel_ignore = 0  # TODO:  hack to keep track of alignments to ref sequences with indels wrt H77

        # stream bowtie2 output, gather and merge reads
        for line in bowtie2_iter:
            if line.startswith('@'):
                # skip header line
                continue

            if cache_file is not None:
                cache_file.write(line)

            # FIXME: DEBUGGING
            #if progress > 100000:
            #    break

            items = line.split('\t')
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = items[:11]
            pos = int(pos)-1  # convert to 0-index
            flag = int(flag)

            # Don't want chimeric or secondary alignments
            if not helper.is_primary(flag) or helper.is_chimeric(flag):
                continue

            if mapq < self.min_mapq:
                continue

            if is_show_progress:  # Increment progress after checking for chimeric/secondary alignments otherwise, can have more alignments than reads.
                progress += 1
                if progress % 5000 == 0:
                    print(self.timestamp('mapped %d of %d' % (progress, nrecords)))


            # We only slice HCV
            if rname == '*' or rname.startswith('hg38'):
                continue

            _, seq1, qual1, inserts = self.apply_cigar(cigar, seq, qual)
            seq2 = '-' * pos + seq1
            qual2 = '!' * pos + qual1

            # TODO:  remove this.  hack to check for 1a ref sequences with indels wrt h77 in NS5a gene
            endpos = pos + len(seq1) - 1
            for indel_ref, indel_ref_start, indel_ref_end in [("HCV-1a-EU781824", 6204, 7544),
                                                              ("HCV-1a-EU781823", 6181, 7539),
                                                              ("HCV-1a-EU256096", 6178, 7530),
                                                              ("HCV-1a-EU256104", 6175, 7530),
                                                              ("HCV-1a-KC283194", 5960, 7304),
                                                              ("HCV-1a-EU482836", 6189, 7529)]:
                if rname == indel_ref and pos >= indel_ref_start and endpos <= indel_ref_end:
                    total_indel_ignore += 1
                    

            # Once we loop, we lose track of this mate.  If the other mate is unmapped, we never output this mate.
            # I.e.  we exclude single mates
            if qname not in read_cache:
                read_cache.update({qname: (rname, seq2, qual2, pos)})
                continue

            # pop mate from cache
            rname1, seq1, qual1, pos1 = read_cache.pop(qname)
            mseq = self.merge_pairs(seq1, seq2, qual1, qual2, q_cutoff=15)
            if mseq.count('N') / float(len(mseq.strip('-'))) > 0.2:
                continue

            # use left-most coordinate between mate1 and mate2 for merged mates
            read_start_1based_wrt_subtype_fullgenome = min(pos, pos1)
            key = (mseq, rname, read_start_1based_wrt_subtype_fullgenome)
            if key not in aligned:
                aligned.update({key: 0})
            aligned[key] += 1

        if cache_file is not None:
            cache_file.close()
        if close_sam:
            bowtie2_iter.close()

        print "ERROR: Total hit 1a ref with indel wrt H77  in NS5A - " + str(total_indel_ignore)

        return aligned


    def slice(self, aligned):
        """
        Slice segments out of aligned reads if they map to correct gene regions.
        Map position should be relative to the reference genome that it mapped to.
        Insertions will have been filtered out by apply_cigar().
        Do not allow reads that do not entirely cover target region.
        :param aligned:
        :return:
        """
        slices = {}
        for (mseq, rname, pos), count in aligned.iteritems():
            # check where this read mapped
            read_start = int(pos)  # 1based coord wrt subtype full genome  corresponding to read start
            read_end = read_start + len(mseq.strip('-'))  # 1based coord wrt subtype full genome corresponding to read end + 1
            if rname not in self.coords:
                print "Read %s didn't have a corresponding map" % rname
                continue
            coords = self.coords[rname]
            if read_end < coords['Core'][0] or read_start > coords['NS5b'][1]:
                # read falls outside of ORF
                continue

            # did it map to one of the target genes?
            for target_gene, target_coords in self.targets.iteritems():
                if target_gene not in slices:
                    slices.update({target_gene: {}})
                left, right = coords[target_gene]  # 0based nuc coordinates wrt subtype full genome corresponding to gene start, gene end+1

                for tc in target_coords:  # 0based nuc coordinates wrt H77 gene corresponding to target start, target end + 1
                    if tc not in slices[target_gene]:
                        slices[target_gene].update({tc: []})
                    # 0based nuc coordinates wrt H77 gene corresponding to target start, target end + 1
                    gene_left, gene_right = tc  # unpack tuple

                    # adjust gene coordinates to genome coordinates.  Assume no indels wrt H77
                    genome_left = gene_left + left  # 0based nuc coordinates wrt subtype full genome corresponding to H77 target start
                    genome_right = gene_right + left # 0based nuc coordinates wrt subtype full genome corresponding to H77 target end + 1

                    # If the read begins after the target region or ends before the target region, read_slice_size will be negative.
                    read_slice_size = min(genome_right, read_end-1) - max(genome_left, read_start-1)
                    target_slice_size = genome_right - genome_left
                    if read_slice_size / float(target_slice_size) < self.min_target_width:
                        # merged read does not cover the minimum required width of the target region
                        continue

                    slice = mseq[genome_left:genome_right]  # mseq should be left-padded wrt subtype full genome

                    # After slicing merged sequences, we may end up with duplicate slices for the same target gene & coords
                    #   against the same reference
                    slices[target_gene][tc].append((rname, slice, count))

        return slices


    def run(self, f1, f2, handle, log, runname='', complete={}, is_show_progress=False, cache=None, min_target_width_cover=1.0):
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

        # Instantiate a cache object if we can

        aligned = self.align(f1, f2, is_show_progress=is_show_progress, cache=cache)
        slices = self.slice(aligned)
        # write out to file
        for gene, subsets in slices.iteritems():
            for coords, subset in subsets.iteritems():
                # 0based nuc coordinates wrt H77 gene corresponding to target start, target end + 1
                target_left_0based_wrt_h77gene, target_right_0based_wrt_h77_gene = coords

                # After slicing merged sequences, we may end up with duplicate slices for the same target gene & coords
                #   against the same reference.  Merge these entries.
                subtype_slice_to_count = {}
                for rname, nucseq, count in subset:
                    subtype = rname.split('-')[1]
                    if not subtype_slice_to_count.get((subtype, nucseq)):
                        subtype_slice_to_count[(subtype, nucseq)] = 0

                    subtype_slice_to_count[(subtype, nucseq)] += count


                # Sort by count
                intermed = sorted(subtype_slice_to_count.items(), key=lambda x:x[1], reverse=True)

                for rank, ((subtype, nucseq), count) in enumerate(intermed):

                    # Assume that all deletions in the nucleotide sequence have been removed
                    # and that there are no deletions prior to the nucleotide sequence
                    # so that the left coordinate is true.
                    frameshift = target_left_0based_wrt_h77gene % 3
                    aaseq = helper.translate(seq=nucseq, shift=frameshift)
                    handle.write('%s,%s,%s,%s,%d,%d,%d,%d,%s,%s,%s\n' % (runname, sample, snum, gene,
                                                                   target_left_0based_wrt_h77gene, target_right_0based_wrt_h77_gene,
                                                                   rank+1, count, subtype, nucseq, aaseq))
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

        # Are we in mpi mode?
        try:
            from mpi4py import MPI
            my_rank = MPI.COMM_WORLD.Get_rank()
            nprocs = MPI.COMM_WORLD.Get_size()
            log += '.' + str(my_rank)

        except ImportError:
            print("Disabling MPI")
            my_rank = 0
            nprocs = 1



        # enable script to restart from interrupted run
        complete = {}

        # if os.path.exists(output) and os.path.getsize(output):
        #     handle = open(output, 'rU')
        #     _ = handle.next()  # skip header
        #     for line in handle:
        #         runname, sample, snum = line.strip('\n').split(',')[:3]
        #         complete.update({(runname, sample, snum): 1})
        #     handle.close()
        #     # re-open file for appending
        #     handle = open(output, 'a')
        # else:
        #     # new output file
        #     handle = open(output, 'w')
        #     handle.write('runname,sample,snum,gene,start,end,rank,count,subtype,nucseq,aaseq\n')

        for path in paths:
            # check that the inputs exist
            if not os.path.exists(path):
                print('No folder/file found at %s', path)
                sys.exit(1)

            if not path.endswith('/'):
                path += '/'

            files = glob(path + '*_R1_001.fastq*')
            if len(files) == 0:
                print 'ERROR: No FASTQ R1 files found at', path
                sys.exit()

            #runname = path.split('/')[4]


            for i, f1 in enumerate(files):
                if i % nprocs != my_rank or my_rank >= len(files):
                    continue
                f2 = f1.replace('_R1_', '_R2_')

                # Hopefully this will give us more consistency with our run names :S
                try:
                    filename = os.path.abspath(f1) 
                    match = re.search('/([0-9]{6}_M[0-9]{5}_[0-9]{4}_[0-9]{9}-[A-Za-z0-9]{5})', filename)
                    runname = match.group(1)
                except:
                    runname = os.path.basename(os.path.dirname(os.path.abspath(f1)))

                sample, snum = os.path.basename(f1).split('_')[:2]
                out_filename = output + "." + sample  + "_" + snum + ".csv"
                cache = None

                # If we're allowed to cache, check to see if the result is in the cache 
                if self.cache_path is not None:
                    cache = helper.Cache(runname, self.min_mapq, self.refpath, "deli", self.cache_path)

                    if cache.check_result(out_filename):
                        cache.decache_result(out_filename)
                        continue

                with open(out_filename, 'w') as handle:
                    handle.write('runname,sample,snum,gene,start,end,rank,count,subtype,nucseq,aaseq\n')
                    self.run(f1, f2, handle, log, runname=runname, complete=complete, is_show_progress=True, cache=cache)
        
                # Cache it
                if self.cache_path is not None:
                    cache.cache_result(out_filename)


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
    parser.add_argument('-min_target_width', type=float, help='minimum fraction of target width covered to be considered a hit', default=1.0)
    parser.add_argument('-targetfile', help='<input> file containing list of folders.  If unspecified, uses predetermined 250bp gene regions best for genotyping')
    parser.add_argument("--cache", help="the cache folder that holds all results and sam files", default=None)

    args = parser.parse_args()
    if not ((args.R1 and args.R2) or args.path or args.pathlist):
        parser.error('Must set one of the following: {-R1 and -R2, -path, -pathlist}.')
    deli = HCVDeli(x=args.x, p=args.p, minlen=args.minlen, minq=args.minq, mins=args.mins, min_target_width=args.min_target_width, targets_file=args.targetfile,
        cache_path=args.cache)

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
        handle.write('runname,sample,snum,gene,start,end,rank,count,subtype,nucseq,aaseq\n')
        deli.run(f1=args.R1, f2=args.R2, handle=handle, log=args.log, is_show_progress=True, cache=args.cache)
        handle.close()
    else:
        # this should never happen!
        parser.print_help()
        sys.exit()



if __name__ == '__main__':
    main()

