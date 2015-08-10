#1/bin/bash

runs=$@
quality=10

for run in $runs; do

    mpirun -hostfile hostfile -v python batch-mpi.py  -path /media/RAW_DATA/MiSeq/runs/$run/Data/Intensities/BaseCalls/ -runname $run -output $run.csv -log $run.log -x data/gb-ref+hg38_v2 -minq $quality #--mapqfile /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/q7/$run.mapq.csv 
    bash catcsv.sh $run.csv.* >  /media/macdatafile/mixed-hcv/gb-ref+hg38_v2/q7/$run.csv

    mpirun -hostfile hostfile -v python batch-mpi.py  -path /media/RAW_DATA/MiSeq/runs/$run/Data/Intensities/BaseCalls/ -runname $run -output $run.csv -log $run.log -x data/gb-ref -minq $quality #--mapqfile /media/macdatafile/mixed-hcv/gb-ref/q7/$run.mapq.csv

    bash catcsv.sh $run.csv.* >  /media/macdatafile/mixed-hcv/gb-ref/q7/$run.csv
done;

