#1/bin/bash

runs=$@
quality=0

for run in $runs; do

    mpirun -hostfile hostfile -v python batch-mpi.py  -path /media/RAW_DATA/MiSeq/runs/$run/Data/Intensities/BaseCalls/ -runname $run -output $run.csv -log $run.log -x data/gb-ref+hg38_v2 -minq $quality #--mapqfile /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/q7/$run.mapq.csv 
    bash catcsv.sh $run.csv.* >  /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/q$quality.$run.csv
    zip gb-ref+hum.$run.ind.zip $run.csv.*
    Rscript launch_knitr_report.R subtype_perc.R /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/q$quality.$run.csv  /media/macdatafile/mixed-hcv/expected_mixture/$run.expected_mixture.csv  $run


    mpirun -hostfile hostfile -v python batch-mpi.py  -path /media/RAW_DATA/MiSeq/runs/$run/Data/Intensities/BaseCalls/ -runname $run -output $run.csv -log $run.log -x data/gb-ref -minq $quality #--mapqfile /media/macdatafile/mixed-hcv/gb-ref/q7/$run.mapq.csv

    bash catcsv.sh $run.csv.* >  /media/macdatafile/mixed-hcv/gb-ref/q$quality.$run.csv
    zip gb-ref.$run.ind.zip $run.csv.*
    Rscript launch_knitr_report.R subtype_perc.R /media/macdatafile/mixed-hcv/gb-ref/q$quality.$run.csv  /media/macdatafile/mixed-hcv/expected_mixture/$run.expected_mixture.csv  $run
done;

