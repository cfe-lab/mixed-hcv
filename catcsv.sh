#!/bin/bash
iter=0

for input in $@ 
do 
    if [ iter = 0 ] ; then
        head -n1 $input;
    fi
    tail -n+2 $input;
done;




#{ head -n1 M00325.csv.D75150-HCV_S20; for f in M00325.csv.D751*; do tail -n+2 "$f"; done; } > new.csv
