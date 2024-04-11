#!/bin/bash

cd ~/Documents/BoneRugosity/SIMSONIC/SimSonic2D_SourceCode/
export OMP_NUM_THREADS=16
SECONDS=0
for k in {19..19}
do
        l=$(printf '%02d' $(( 5*k - 4)))
        m=$(printf '%02d' $(( 5*k - 3)))
        n=$(printf '%02d' $(( 5*k - 2)))
        o=$(printf '%02d' $(( 5*k - 1)))
        p=$(printf '%02d' $(( 5*k - 0)))

	echo $1/'tx_'$l'/'
        ./simsonic2D $1/'tx_'$l'/' &
        ./simsonic2D $1/'tx_'$m'/' &
        ./simsonic2D $1/'tx_'$n'/' &
        ./simsonic2D $1/'tx_'$o'/' &
        ./simsonic2D $1/"tx_"$p"/" 
        wait
done

eval "echo Elapsed time: $(date -ud "@$SECONDS" +'$((%s/3600/24)) days %H hr %M min %S sec')"
