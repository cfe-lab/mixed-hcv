"""
Generate gene coordinates for gb reference set
"""

from seqUtils import convert_fasta
import HyPhy
import hyphyAlign as align
import os

hyphy = HyPhy._THyPhy (os.getcwd(), 1) # instance of HyPhy

# configure for nucleotide alignment
align.change_settings(hyphy, 
                      alphabet = align.nucAlphabet, 
                      scoreMatrix = align.nucScoreMatrix, 
                      gapOpen = 20,
                      gapOpen2 = 20,
                      gapExtend = 10,
                      gapExtend2 = 10,
                      noTerminalPenalty = 1)

	
handle = open('../data/gb-ref2.fa', 'rU')
fasta = convert_fasta(handle)
handle.close()

handle = open('../data/h77-genes.fa', 'rU')
h77 = convert_fasta(handle)
handle.close()

outfile = open('../data/gb-ref2.coords', 'w')

for h, s in fasta:
    print h
    for gene, ref in h77:
        # locate the H77 gene in this reference genome
        aquery, aref, ascore = align.pair_align(hyphy, ref, s)
        left, right = align.get_boundaries(aref)
        aquery2 = aquery[left:right].replace('-', '')
        coord1 = s.index(aquery2)
        coord2 = coord1 + len(aquery2)
        
        outfile.write('%s,%s,%d,%d\n' % (h, gene, coord1, coord2))

outfile.close()

