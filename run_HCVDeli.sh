#!/bin/bash

runs=$@
quality=0
minwidth=0.5
ref=gb-ref+hg38_v2


for run in $runs; do

    mkdir -p ./log/$run
    mkdir -p ./out/$run


    mpirun -hostfile ./out/machines.txt --cache /media/macdatafile/mixed-hcv/cache/ -output-filename ./log/$run/$run__deli__HCVhuman__q$quality__mw$minwidth.stdout.log -v python HCVDeli.py  -path ./reads/$run  -log ./log/$run/$run__deli__HCVhuman__q$quality__mw$minwidth.log -x data/$ref -minq $quality  -min_target_width $minwidth ./out/$run/$run__deli__HCVhuman__q$quality__mw$minwidth.csv

    bash catcsv.sh ./out/$run/$run__deli__HCVhuman__q$quality__mw$minwidth.csv.* >  /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/deli/$run__deli__HCVhuman__q$quality__mw$minwidth.csv


    Rscript launch_knitr_report.R  subtype_slice_perc.R  /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/deli/$run__deli__HCVhuman__q$quality__mw$minwidth.csv  /media/macdatafile/mixed-hcv/expected_mixture/$run.expected_mixture.csv  $run__deli__HCVhuman__q$quality__mw$minwidth.csv

done;
