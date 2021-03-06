#!/bin/bash

runs=$@
qualities=(0 7 10 20)
minwidths=(0.125)
ref=gb-ref+hg38_v2
refnice=HCV_Human

for run in ${runs[@]}; do

    mkdir -p ./log/$run
    mkdir -p ./out/$run

    for quality in ${qualities[@]}; do
        for minwidth in ${minwidths[@]}; do
        	echo "run=$run quality=$quality minwidth=$minwidth ref=$ref refnice=$refnice"
            mpirun -hostfile ./out/machines.txt  -output-filename ./log/$run/${run}__deli__${refnice}__q${quality}__mw${minwidth}.stdout.log -v python HCVDeli.py  -path ./reads/$run  --cache /media/macdatafile/mixed-hcv/cache/ -log ./log/$run/${run}__deli__${refnice}__q${quality}__mw${minwidth}.log -x data/$ref -minq $quality  -min_target_width $minwidth ./out/$run/${run}__deli__${refnice}__q${quality}__mw${minwidth}.csv


            python collectcat.py -runname ${run} -output /media/macdatafile/mixed-hcv/staging/ -runtype deli -x ${ref} -minq ${quality} -min_target_width ${minwidth} --cache /media/macdatafile/mixed-hcv/cache/ ./out/$run/${run}__deli__${refnice}__q${quality}__mw${minwidth}.csv.*


            Rscript launch_knitr_report.R  subtype_slice_perc.R  /media/macdatafile/mixed-hcv/mixture_reports/${run}__deli__${refnice}__q${quality}__mw${minwidth}.subtype_slice_perc.html  /media/macdatafile/mixed-hcv/staging/$run/${run}__deli__${refnice}__q${quality}__mw${minwidth}.csv  /media/macdatafile/mixed-hcv/expected_mixture/$run.expected_mixture.csv  ${run}__deli__${refnice}__q${quality}__mw${minwidth}
            
            rm ./out/$run/${run}__deli__${refnice}__q${quality}__mw${minwidth}.csv*
        done;
    done;
done;
