# mixed-hcv
Detects mixed genotype HCV infections by mapping short read data to reference genomes

**Full-Genome Hit Pipeline:**

Aligns reads against reference fasta and counts hits.  Outputs sam file including unaligned reads, secondary alignments, and chimeric alignments.  Any first-mate that aligns to HCV\* or hg38\* reference sequence counts as a hit.  Ignores low mapping quality, low alignment score (AS field) alignments.

  

- Execute run.sh to run entire pipeline.

  - Batch-mpi.py:  Aligns reads, counts hits.  Caches sam, per-sample hit count CSV to hard-drive.  Will not overwrite existing cached files.  You must execute this with MPI.
  
  - collectcat.py:  Collates all per-sample hit count CSV files into a single per-run hit count CSV.


**Target-Region Hit Pipeline:**

Aligns reads against reference fasta, counts hits against target-regions, and outputs sequences that hit target regions.  Outputs sam file excluding unaligned reads, and reads where both mates do not map to the same reference.  If the merged sequence from both mates covers the target region over the minimum width threshold, it is considered a hit.  Ignores secondary, chimeric, low mapping quality alignments.  Target regions are defined in /mixed-hcv/data/gb-ref2.coords by default, but can be overridden via commandline.  The target region coordinates CSV file specifies the 0-based target region coordinates with respect to the reference sequence.

- Execute run_HCVDeli.sh to run entire pipeline.

  - HCVDeli.py:  Aligns reads, merges mates, counts hits.  Caches sam, per-sample target-region sequence count CSV to hard-drive.  Will not overwrite existing cached files. You may execute this script with or without MPI.
  
  - collectcat.py:  Collates all per-sample hit target-region sequence count CSV files into per-run sequence count CSV.
  


**Reference Fastas:**

- gb-ref+hg38_v2: The HCV+Human Genome+Human Mitochondria reference fasta leads to fewest false hits, but at 3GB, it is too large to upload to github.  Stored in /macdatafile/mixed-hcv/ref/gb-ref+hg38_v2.fa

- gb-ref: Only HCV sequences from genbank.  Some of the sequences don't have accessions.  Stored in /macdatafile/mixed-hcv/ref/gb-ref.fa

- gb-ref2:  Only HCV sequences from genbank.  All sequences will have accessions.  Stored in /macdatafile/mixed-hcv/ref/gb-ref2.fa

  
