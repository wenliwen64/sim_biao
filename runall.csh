#!/bin/sh

random=$1

while [ "$random" -lt $2 ];
do
    qsub -hard -l h_vmem=3G -l projectio=1,scratchfree=500 -o ./log/${ifile}.out -e ./log/${ifile}.err sim.csh $random  
    echo "The $random task has been submitted!"
    sleep 0.1
    let "random+=1"
done
