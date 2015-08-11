#!/bin/bash

runs=$@
quality=0

for run in $runs; do

mpirun --machinefile ~/gitrepo/mixed-hcv/out/machines.txt --output-filename ~/gitrepo/mixed-hcv/log/deli/HCVDeli.log python HCVDeli.py -pathlist ~/gitrepo/mixed-hcv/reads/allreads.txt -log ~/gitrepo/mixed-hcv/log/deli/HCVDeli.mpi.log -x ~/gitrepo/mixed-hcv/data/gb-ref+hg38_v2 ~/gitrepo/mixed-hcv/out/deli/HCVDeli.mpi.out.csv

    mpirun -hostfile ./out/machines.txt -output-filename ./log/$run/$run.deli.stdout.log -v python HCVDeli.py  -path ./reads/$run  -log ./log/$run/$run.deli.log -x data/gb-ref+hg38_v2 -minq $quality  ./out/$run/$run.deli.csv

    bash catcsv.sh ./out/$run/$run.deli.csv.* >  /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/deli/$run.csv


    Rscript launch_knitr_report.R  subtype_slice_perc.R  /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/deli/$run.csv  /media/macdatafile/mixed-hcv/expected_mixture/$run.expected_mixture.csv  $run

done;
