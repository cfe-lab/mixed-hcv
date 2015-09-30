# batch script for applying sam2aln.py to the cached SAM files
# corresponding to the Merck HCV runs

from glob import glob
import os

nthreads = 6 # turn off multithreading
samfiles = glob('/Volumes/MACDATAFILE/mixed-hcv/cache/*/sam/*/*.sam')
dest = '/Users/art/git/mixed-hcv/working/coverage/'
targets = ['150720_M01841_0148_000000000-AE93J', '150724_M01841_0150_000000000-AE8E0']

for f in samfiles:
    # contains bowtie2 arguments
    _, _, _, _, _, runname, _, dbname, filename = f.split('/')
    if runname not in targets:
        continue
    if dbname != 'gb-ref+hg38_v2':
        continue

    print runname, dbname, filename

    tokens = filename.split('_')
    _, sample, snum = tokens[:3]
    bowtie2_args = tokens[-1].replace('.sam', '')

    if not os.path.exists(dest+runname):
        os.mkdir(dest+runname)

    # subfolder for reference database used
    dest2 = os.path.join(dest, runname, dbname)
    if not os.path.exists(dest2):
        os.mkdir(dest2)

    prefix = '%s_%s_%s' % (sample, snum, bowtie2_args)
    aligned_csv = os.path.join(dest2, prefix+'.aligned.csv')
    insert_csv = aligned_csv.replace('.aligned.csv', '.insert.csv')
    failed_csv = aligned_csv.replace('.aligned.csv', '.failed.csv')
    if os.path.exists(aligned_csv):
        print 'skipping'
        continue

    os.system('python sam2aln.py %s %s %s %s -p %d' % (f, aligned_csv, insert_csv, failed_csv, nthreads))
    #break  # debugging
