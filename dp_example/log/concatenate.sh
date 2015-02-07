#!/bin/bash

ONELINES_DIR="./onelines"
if [ ! -d $ONELINES_DIR ]
then
    mkdir -p $ONELINES_DIR
fi

N=100
NA=1
T=`expr $N \* $N / 2`
NUM_BATCHES=10
EPS="1e+20"
MAX_IT=120

SRF_MIN=0.2
SRF_MAX=0.2
SRF_INC=0.1

MF_MIN=0.2
MF_MAX=2.0
MF_INC=0.2

START=1
END=50

PREFIX="e_emsf_"$N"_"$NA"_"$T"_"$NUM_BATCHES"_"$EPS"_"$MAX_IT

SR_MIN=$(echo "$SRF_MIN * $N" | bc -l)
SR_MIN=`printf %.0f%s ${SR_MIN%}`

SR_MAX=$(echo "$SRF_MAX * $N" | bc -l)
SR_MAX=`printf %.0f%s ${SR_MAX%}`

SR_INC=$(echo "$SRF_INC * $N" | bc -l)
SR_INC=`printf %.0f%s ${SR_INC%}`

SR=$SR_MIN
while [ $SR -le $SR_MAX ]
do
    SR=`printf %04d%s ${SR%}`

    M_MIN=$(echo "$MF_MIN * $SR" | bc -l)
    M_MIN=`printf %.0f%s ${M_MIN%}`

    M_MAX=$(echo "$MF_MAX * $SR" | bc -l)
    M_MAX=`printf %.0f%s ${M_MAX%}`

    M_INC=$(echo "$MF_INC * $SR" | bc -l)
    M_INC=`printf %.0f%s ${M_INC%}`

    PREFIX_SR=$PREFIX"_"$SR
    PLOT_FILENAME="plot_"$PREFIX_SR".R"
    echo "D <- NULL" > $PLOT_FILENAME
    echo >> $PLOT_FILENAME

    M=$M_MIN
    while [ $M -le $M_MAX ]
    do
        M=`printf %04d%s ${M%}`

        PREFIX_SR_M=$PREFIX"_"$SR"_"$M
        PREFIX_SR_M_LOG=$PREFIX_SR_M".log"

        if [ -f $PREFIX_SR_M_LOG ]
        then
            rm $PREFIX_SR_M_LOG
        fi

        ALL_EXIST=1
        for RUN in $(eval echo {$START..$END})
        do
            PREFIX_SR_M_RUN_LOG=$PREFIX_SR_M"_"`printf %02d%s ${RUN%}`".log"
            if [ -f $PREFIX_SR_M_RUN_LOG ]
            then
                (cat $PREFIX_SR_M_RUN_LOG; echo) >> $PREFIX_SR_M_LOG
            else
                ALL_EXIST=0
            fi
        done

        if [ $ALL_EXIST -eq 1 ]
        then
            echo "Generated    $PREFIX_SR_M_LOG"

            echo "tmp <- read.table(\"$PREFIX_SR_M_LOG\")" >> $PLOT_FILENAME
            echo "tmp <- apply(tmp, 2, mean)" >> $PLOT_FILENAME
            echo "D <- cbind(D, tmp)" >> $PLOT_FILENAME
            echo >> $PLOT_FILENAME
        fi            
        
        M=`expr $M + $M_INC`
    done

    echo "matplot(D, t=\"o\", log=\"y\", main=\"sr=$SR; maxit=$MAX_IT; D, K aleat.\")" >> $PLOT_FILENAME
    echo >> $PLOT_FILENAME

    SR=`expr $SR + $SR_INC`
    echo
done
