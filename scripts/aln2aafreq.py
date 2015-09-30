"""
Convert *.aligned.csv file into amino acid frequencies by:
(1) generating reference genome-specific amino acid frequency
  tables
(2) aligning each consensus against H77
(3) merging tables into a single frequency table in H77
  coordinate space
Do this for NS3, NS5a and NS5b
"""
import sys
import os
from seqUtils import translate_nuc
from itertools import groupby
from csv import DictReader
from glob import glob

# settings
files = glob('../working/coverage/15072[04]*/*/*.aligned.csv')
targets = ['NS3', 'NS5a', 'NS5b']
aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
          'R', 'S', 'T', 'V', 'W', 'Y', '*']


# load gene coordinates (0-index as per Python's str.index())
gene_coords = {}
with open('../data/gb-ref2.coords', 'rU') as f:
    for line in f:
        refname, gene, left, right = line.strip('\n').split(',')
        if gene not in targets:
            continue
        if refname not in gene_coords:
            gene_coords.update({refname: {}})
        gene_coords[refname].update({gene: (int(left), int(right))})



for f in files:
    reader = DictReader(open(f, 'rU'))
    print f
    filename = os.path.basename(f)

    outfile = open(f.replace('.aligned.csv', '.amino.csv'), 'w')
    outfile.write('rname,gene,aa.pos,' + ','.join(aminos) + '\n')

    for rname, group in groupby(reader, lambda row: row['refname']):
        counts = dict([(target, {}) for target in targets])

        # gather amino acid counts
        for row in group:
            start = int(row['offset'])
            seq = row['seq']
            end = start + len(seq)
            for gene in targets:
                left, right = gene_coords[rname][gene]

                if start > right or end < left:
                    # no overlap with gene
                    continue

                seq2 = '-'*start + seq
                clip = seq2[left:right]  # should be in frame
                aaseq = translate_nuc(clip.upper(), 0)
                for pos, aa in enumerate(aaseq):
                    if aa in ['-', '?']:
                        continue
                    if pos not in counts[gene]:
                        counts[gene].update({pos: dict(zip(aminos, [0 for a in aminos]))})
                    counts[gene][pos][aa] += int(row['count'])

        for gene, group in counts.iteritems():
            for pos, group2 in group.iteritems():
                astr = ','.join(map(str, [group2[a] for a in aminos]))
                outfile.write('%s,%s,%d,%s\n' % (rname, gene, pos, astr))

    outfile.close()
    #break  # run on one file only


