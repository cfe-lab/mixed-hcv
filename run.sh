#1/bin/bash

runs=$@
qualities=(0 10)
refs=(gb-ref+hg38_v2 gb-ref2 gb-ref)

for run in ${runs[@]}; do
    for ref in ${refs[@]}; do
        for quality in ${qualities[@]}; do

    echo mpirun -hostfile hostfile -v python batch-mpi.py --cache ./cache/ -path /media/RAW_DATA/MiSeq/runs/$run/Data/Intensities/BaseCalls/ -runname $run -output $run.csv -log $run.log -x data/$ref -minq $quality  
    mpirun -hostfile hostfile -v python batch-mpi.py --cache ./cache/  -path /media/RAW_DATA/MiSeq/runs/$run/Data/Intensities/BaseCalls/ -runname $run -output $run.csv -log $run.log -x data/$ref -minq $quality  
    echo ""
    echo collectcat.py -runname ${run} -runtype full -output /media/macdatafile/mixed-hcv/staging/ -x ${ref} -minq ${quality} -min_target_width "0" --cache /media/macdatafile/mixed-hcv/cache/ $run.csv.* 
    python collectcat.py -runname ${run} -runtype full -output /media/macdatafile/mixed-hcv/staging/ -x ${ref} -minq ${quality} -min_target_width "0" --cache ./cache/ $run.csv.* 
echo ""
echo ""
        done
    done
done
