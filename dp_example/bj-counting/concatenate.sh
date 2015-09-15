#!/bin/bash

NUM_BATCHES=30000
NUM_EPISODES=1000000
MIN_BATCHES=0
NUM_POINTS=100
N=203
SR=20
NA=2

START=1
END=30
GROUP="v_cnt_bj_"$SR"_"$NA"_"$NUM_BATCHES"_"$NUM_EPISODES"_"$MIN_BATCHES"_"$NUM_POINTS
PLOT_FILE="plot_"$GROUP".R"

echo "D <- NULL" > $PLOT_FILE
echo "legends=list()" >> $PLOT_FILE
echo >> $PLOT_FILE
for EPSILON in "0.01" "0.05" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1"
do
    GROUP_FILE=$GROUP"_"$EPSILON".log"
    if [ -f $GROUP_FILE ]
    then
        rm $GROUP_FILE
    fi

    ALL_EXIST=1
    for RUN in $(seq $START $END)
    do
        RUN_FILE=$GROUP"_"$EPSILON"_"$(printf %02d%s ${RUN%})".log"
        if [ -f $RUN_FILE ]
        then
            (cat $RUN_FILE; echo) >> $GROUP_FILE
        else
            ALL_EXIST=0
            echo "Está faltando arquivo de resultado: $RUN_FILE"
        fi
    done

    echo "tmp <- read.table(\"$GROUP_FILE\")" >> $PLOT_FILE
    echo "tmp <- apply(tmp, 2, mean)" >> $PLOT_FILE
    echo "D <- cbind(D, tmp)" >> $PLOT_FILE
    echo "legends[[length(legends)+1]]=paste(\""$EPSILON"\")" >> $PLOT_FILE
    echo >> $PLOT_FILE
done
echo "matplot(D, t=\"o\")" >> $PLOT_FILE
echo "grid()" >> $PLOT_FILE
echo "legend(\"topright\", legend=legends, col=cols, lwd=2, bty=\"n\")" >> $PLOT_FILE
