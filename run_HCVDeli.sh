#!/bin/bash

runs=$@
qualities=(0 10)
minwidths=(0.25 0.5 0.75 1.0)
ref=gb-ref+hg38_v2
refnice=$refnice

for run in $runs; do

    mkdir -p ./log/$run
    mkdir -p ./out/$run

    for quality in $qualities; do
        for minwidth in minwidths; do
            mpirun -hostfile ./out/machines.txt  -output-filename ./log/$run/$run__deli__$refnice__q$quality__mw$minwidth.stdout.log -v python HCVDeli.py  -path ./reads/$run  --cache /media/macdatafile/mixed-hcv/cache/ -log ./log/$run/$run__deli__$refnice__q$quality__mw$minwidth.log -x data/$ref -minq $quality  -min_target_width $minwidth ./out/$run/$run__deli__$refnice__q$quality__mw$minwidth.csv

            bash catcsv.sh ./out/$run/$run__deli__$refnice__q$quality__mw$minwidth.csv.* >  /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/deli/$run__deli__$refnice__q$quality__mw$minwidth.csv


            Rscript launch_knitr_report.R  subtype_slice_perc.R  /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/deli/$run__deli__$refnice__q$quality__mw$minwidth.csv  /media/macdatafile/mixed-hcv/expected_mixture/$run.expected_mixture.csv  $run__deli__$refnice__q$quality__mw$minwidth.csv
        done;
    done;
done;
