#!/bin/bash

cd ~/Documents/BoneRugosity/SIMSONIC/SimSonic2D_SourceCode/
export OMP_NUM_THREADS=16
SECONDS=0

./simsonic2D $1/tx_01/
eval "echo Elapsed time: $(date -ud "@$SECONDS" +'$((%s/3600/24)) days %H hr %M min %S sec')"
