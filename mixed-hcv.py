__author__ = 'art'

import os
import atexit

import Tkinter as tk
import ttk
import tkFileDialog

from Bio import SeqIO
from settings import *

import subprocess
import re
import time



class App:
    def __init__(self, master):
        self.__version__ = '0.1'
        self.master = master

        frame = tk.Frame(master)  # simple container to hold other widgets
        frame.pack()  # resize to contents and make visible

        self.button_load = tk.Button(
            frame, text="Select FASTQ", command=self.open_files
        )
        self.button_load.grid(row=0, column=0, sticky='W')

        self.button_quit = tk.Button(
            frame, text='Quit', fg='red', command=frame.quit
        )
        self.button_quit.grid(row=0, column=1, sticky='E')

        self.progress_bar = ttk.Progressbar(
            frame, orient='horizontal', length=500, mode='determinate'
        )
        self.progress_bar.grid(row=1, columnspan=2)


        self.console = tk.Text(
            frame, height=30, width=80, bg='black', fg='white', cursor='xterm'
        )
        self.console.grid(row=2, columnspan=2, sticky='W')

        #frame2 = tk.Frame(master, bd=10, relief=tk.GROOVE)
        #frame2.pack()
        #self.canvas = tk.Canvas(
        #    frame2, height=300, width=560, bd=1
        #)
        #self.canvas.grid(row=3, columnspan=2)
        #self.canvas.create_rectangle(0, 0, 560, 300)  # from top left


        self.console.insert(tk.END, 'mixed-hcv v%s\n' % self.__version__)
        self.console.insert(tk.END, '===========%s\n' % ('='*len(self.__version__)))
        self.console.insert(tk.END, 'Choose a FASTQ R1 file by clicking on "Select FASTQ" button.  The program\n'
                                    'automatically locates the corresponding R2 file.  If the files are valid then\n'
                                    'mixed-hcv will map the reads to the HCV genotype reference set.\n\n')


        self.fastq1 = ''
        self.fastq2 = ''

        self.cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token


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


    def do_map(self, fastq1, fastq2, refpath, bowtie_threads, min_match_len, min_mapq):
        """
        Process SAM output as it is streamed from bowtie2 to assign
        short reads to HCV genotypes/subtypes.
        """

        # get size of FASTQ
        nrecords = self.count_file_lines(fastq1) / 2
        self.progress_bar['value'] = 0
        self.progress_bar['maximum'] = nrecords

        # stream output from bowtie2
        bowtie_args = ['bowtie2',
                       '--quiet',
                       '-x', refpath,
                       '-1', fastq1,
                       '-2', fastq2,
                       #'--no-unal',
                       '--local', '--sensitive-local',
                       '-p', str(bowtie_threads)]

        p = subprocess.Popen(bowtie_args, stdout=subprocess.PIPE)
        #atexit.register(p.terminate)

        # collect SAM output by refname
        genotype_counts = dict([(str(i+1), 0) for i in range(7)])
        counts = {}
        n_short = 0
        n_hybrid = 0
        n_mapq = 0
        progress = 0
        total_count = 0

        with p.stdout:
            # FIXME: This locks up Tkinter event handler - split into separate thread
            for line in p.stdout:
                if line.startswith('@'):
                    # skip header line
                    continue

                qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.split('\t')[:11]
                progress += 1
                if progress % progress_update_interval == 0:
                    self.progress_bar['value'] = progress
                    #self.console.insert(tk.END, '| 1: %d | 2: %d | 3: %d | 4: %d | 5: %d | 6: %d '
                    #                            '| 7: %d |\n' % (genotype_counts['1'],
                    #                                             genotype_counts['2'], genotype_counts['3'],
                    #                                             genotype_counts['4'], genotype_counts['5'],
                    #                                             genotype_counts['6'], genotype_counts['7']))
                    self.master.update_idletasks()
                    #self.console.see(tk.END)

                # read failed to map to anything
                if rname == '*':
                    continue

                subtype = rname.split('-')[1]  # 'HCV-1a-[accno]' -> '1a'
                genotype = subtype[0]

                # ignore second reads
                if not self.is_first_read(flag):
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
                tokens = self.cigar_re.findall(cigar)
                match_len = sum([int(token.strip('M')) for token in tokens if token.endswith('M')])
                if match_len < min_match_len:
                    n_short += 1
                    continue

                # update counts
                if subtype not in counts:
                    counts.update({subtype: 0})
                counts[subtype] += 1
                genotype_counts[genotype] += 1
                total_count += 1

            if p.returncode:
                raise subprocess.CalledProcessError(p.returncode, bowtie_args)

        # output results
        keys = counts.keys()
        keys.sort()

        return counts, nrecords, total_count, n_short, n_mapq, n_hybrid

    def is_valid_fastq(self, path):
        if path == '':
            # user called cancel
            return False

        handle = SeqIO.parse(path, format='fastq')
        try:
            _ = handle.next()
        except:
            self.console.insert(tk.END, 'File %s does not appear to be FASTQ format\n' % path)
            self.console.see(tk.END)
            return False
        handle.close()
        return True

    def open_files(self):
        """
        User specifies a FASTQ file containing first reads.
        Search this location in filesystem for the complementary
        second read set.
        """
        self.fastq1 = tkFileDialog.askopenfilename()
        if not self.is_valid_fastq(self.fastq1):
            return

        self.console.insert(tk.END, 'Opened %s\n' % self.fastq1)
        self.console.see(tk.END)

        self.fastq2 = self.fastq1.replace('_R1_', '_R2_')
        if not os.path.exists(self.fastq2):
            self.console.insert(tk.END, 'Failed to locate R2 FASTQ file\n')
            self.console.see(tk.END)
            return

        if not self.is_valid_fastq(self.fastq2):
            return

        self.console.insert(tk.END, 'Opened %s\n' % self.fastq2)
        self.console.see(tk.END)
        self.master.update_idletasks()

        # okay to process
        t0 = time.time()
        try:
            counts, total_count, total_mapped, n_short, n_mapq, n_hybrid = self.do_map(self.fastq1, self.fastq2,
                                                                                       refpath, bowtie_threads,
                                                                                       min_match_len, min_mapq)
        except Exception, e:
            self.console.insert(tk.END, 'ERROR: %s\n' % e)
            self.console.see(tk.END)
        elapsed = time.time() - t0

        # spool output to console
        self.console.insert(tk.END, '===========================\n')
        self.console.insert(tk.END, 'Total number of pairs: %d\n' % (int(total_count)/2,))
        self.console.insert(tk.END, 'Total pairs mapped: %d\n' % total_mapped)
        self.console.insert(tk.END, '-------------------------\n')

        keys = counts.keys()
        keys.sort()
        for key in keys:
            self.console.insert(tk.END, '%s: %d (%1.2f%%)\n' % (key, counts[key],
                                                                100.*counts[key]/float(total_mapped)))

        self.console.insert(tk.END, '-------------------------\n')
        self.console.insert(tk.END, 'Reads too short: %d\n' % n_short)
        self.console.insert(tk.END, 'Low map quality: %d\n' % n_mapq)
        self.console.insert(tk.END, 'Discordant pair mapping: %d\n' % n_hybrid)
        self.console.insert(tk.END, '===========================\nDONE - running time %d seconds\n\n' % elapsed)
        self.console.see(tk.END)



root = tk.Tk()  # parent widget
root.wm_title('mixed_hcv')

app = App(root)

root.mainloop()  # enter Tk event loop
root.destroy()  # only required in certain dev environments
