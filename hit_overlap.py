# For each run, finds the hits that overlap between Full-Genome (batch-mpi.py) pipelines run with different parameters.
# Calls R to make venn diagrams
import os
import subprocess
import helper

RUNS=['150720_M01841_0148_000000000-AE93J',
      '150729_M01841_0152_000000000-AE8B3',
      '150802_M01841_0154_000000000-AE8FV',
      '150806_M01841_0156_000000000-AE8F9',
      '150724_M01841_0150_000000000-AE8E0',
      '150731_M01841_0153_000000000-AE96E',
      '150804_M01841_0155_000000000-AE95T']

#REFS=['gb-ref+hg38_v2', 'gb-ref2', 'gb-ref']
REFS=['gb-ref+hg38_v2']

MAPQS=[0, 10]

CONCAT_HIT_CSV = "/media/macdatafile/mixed-hcv/staging/hits.csv"


with open(CONCAT_HIT_CSV, 'w') as fh_out:
    fh_out.write("run,ref,sample,read,subtype,mapq\n")
    for run in RUNS:
        for ref in REFS:
                # /media/macdatafile/mixed-hcv/cache/150720_M01841_0148_000000000-AE93J/sam/gb-ref+hg38_v2/F_56585A-HCV_S1_L001_R1_001.fastq.gz_R56585A-HCV_S1_L001_R2_001.fastq.gz_--quiet--local--no-unal--no-mixed.sam
            samdir = "/media/macdatafile/mixed-hcv/cache/" + run + "/sam/" + ref
            if not os.path.exists(samdir):
                print "WARN:  missing sam cache dir " + samdir
                continue

            sams=os.listdir(samdir)

            for sam in sams:
                print "run=" + run + " ref=" + ref + " sam=" + sam
                sample = sam.split("_L001_R")[0]
                # Only first reads are counted in full genome hits
                sampipe = subprocess.Popen(["samtools", "view",
                                       "-S",  # input is sam
                                       "-f", str(helper.SamFlag.IS_FIRST),  # only spit out alignments that have this flag set
                                       "-T", "/media/macdatafile/mixed-hcv/ref/" + ref + ".fa",  # since sam files don't have headers, we need to feed it the reference fasta
                                       samdir + os.sep + sam],
                                      stdout=subprocess.PIPE)
                for line in sampipe.stdout:
                    items = line.split('\t')
                    qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = items[:11]
                    if rname.startswith("hg38"):
                        rname_nice="Human"
                    # We only care about differentiating 1a and 1b subtypes.  For the rest, we only need genotype granularity
                    elif rname.startswith("HCV-1a"):
                        rname_nice="1a"
                    elif rname.startswith("HCV-1b"):
                        rname_nice="1b"
                    elif rname.startswith("HCV"):
                        subtype = rname.split("-")[1][0]
                        genotype = subtype[0]
                        rname_nice = genotype
                    else:
                        rname_nice="Unaligned"


                    fh_out.write(",".join([run, ref, sample, qname, rname_nice, mapq]) + "\n")


subprocess.check_call(["Rscript",
                      "launch_knitr_report.R",
                      "hit_overlap.R",
                      "/media/macdatafile/mixed-hcv/mixture_reports/hit_overlap.html",
                      CONCAT_HIT_CSV
])