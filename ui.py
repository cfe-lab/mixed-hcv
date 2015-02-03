__author__ = 'art'

import mapper
from Tkinter import *
import tkFileDialog
import os
import sys
from Bio import SeqIO
from settings import *

class App:
    def __init__(self, master):
        self.__version__ = '0.1'

        frame = Frame(master)  # simple container to hold other widgets
        frame.pack()  # resize to contents and make visible

        self.button_load = Button(
            frame, text="Select FASTQ", command=self.open_files
        )
        self.button_load.grid(row=0, column=0, sticky='W')

        self.button_quit = Button(
            frame, text='Quit', fg='red', command=frame.quit
        )
        self.button_quit.grid(row=0, column=1, sticky='E')

        self.console = Text(
            frame, height=26, width=80, bg='black', fg='white', cursor='xterm'
        )
        self.console.grid(row=1, columnspan=2)

        self.console.insert(END, 'mixed-hcv Tkinter interface v%s\n' % self.__version__)
        self.console.insert(END, '================================\n')
        self.console.insert(END, 'Select a FASTQ R1 file using the "Open file" dialog.  The program will\n'
                                 'automatically locate corresponding R2 file.  If the files are valid then\n'
                                 'mixed-hcv will map the reads to the HCV genotype reference set.\n\n')


        self.fastq1 = ''
        self.fastq2 = ''

    def is_valid_fastq(self, path):
        handle = SeqIO.parse(path, format='fastq')
        try:
            _ = handle.next()
        except:
            self.console.insert(END, 'File %s does not appear to be FASTQ format\n' % path)
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

        self.console.insert(END, 'Opened %s' % self.fastq1)

        self.fastq2 = self.fastq1.replace('_R1_', '_R2_')
        if not os.path.exists(self.fastq2):
            self.console.insert(END, 'Failed to locate R2 FASTQ file\n')
            return

        if not self.is_valid_fastq(self.fastq2):
            return

        # okay to process
        self.process_files()


    def process_files(self):
        """
        Call mapper module
        :return:
        """





root = Tk()  # parent widget
root.wm_title('mixed_hcv')

app = App(root)

root.mainloop()  # enter Tk event loop
root.destroy()  # only required in certain dev environments