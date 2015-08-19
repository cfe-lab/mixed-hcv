#!/bin/bash

runs=$@
qualities=(0 10)
minwidths=(0.25 0.5 0.75 1.0)
ref=gb-ref+hg38_v2
refnice=HCVhuman

for run in $runs; do

    mkdir -p ./log/$run
    mkdir -p ./out/$run

    for quality in $qualities; do
        for minwidth in $minwidths; do
        	echo "run=$run quality=$quality minwidth=$minwidth ref=$ref refnice=$refnice"
            mpirun -hostfile ./out/machines.txt  -output-filename ./log/$run/${run}__deli__${refnice}__q${quality}__mw${minwidth}.stdout.log -v python HCVDeli.py  -path ./reads/$run  --cache /media/macdatafile/mixed-hcv/cache/ -log ./log/$run/${run}__deli__${refnice}__q${quality}__mw${minwidth}.log -x data/$ref -minq $quality  -min_target_width $minwidth ./out/$run/${run}__deli__${refnice}__q${quality}__mw${minwidth}.csv

            bash catcsv.sh ./out/$run/${run}__deli__${refnice}__q${quality}__mw${minwidth}.csv.* >  /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/deli/${run}__deli__${refnice}__q${quality}__mw${minwidth}.csv


            Rscript launch_knitr_report.R  subtype_slice_perc.R  /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/deli/${run}__deli__${refnice}__q${quality}__mw${minwidth}.csv  /media/macdatafile/mixed-hcv/expected_mixture/$run.expected_mixture.csv  ${run}__deli__${refnice}__q${quality}__mw${minwidth}.csv
        done;
    done;
done;
