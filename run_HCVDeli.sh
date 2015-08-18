#!/bin/bash

runs=$@
quality=0
minwidth=0.5
ref=gb-ref+hg38_v2

for run in $runs; do

    mkdir -p ./log/$run
    mkdir -p ./out/$run

    mpirun -hostfile ./out/machines.txt -output-filename ./log/$run/q$quality.$ref.minwidth$minwidth.$run.deli.stdout.log -v python HCVDeli.py  -path ./reads/$run  -log ./log/$run/q$quality.$ref.minwidth$minwidth.$run.deli.log -x data/$ref -minq $quality  -min_target_width $minwidth ./out/$run/q$quality.$ref.minwidth$minwidth.$run.deli.csv

    bash catcsv.sh ./out/$run/q$quality.$ref.minwidth$minwidth.$run.deli.csv.* >  /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/deli/q$quality.$ref.minwidth$minwidth.$run.csv


    Rscript launch_knitr_report.R  subtype_slice_perc.R  /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/deli/q$quality.$ref.minwidth$minwidth.$run.csv  /media/macdatafile/mixed-hcv/expected_mixture/$run.expected_mixture.csv  q$quality.$ref.minwidth$minwidth.$run

done;
