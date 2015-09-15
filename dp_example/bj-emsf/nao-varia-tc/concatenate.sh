#!/bin/bash

N=203
SR=20
NA=2
NUM_BATCHES=5000
NUM_EPISODES=1000000
MIN_BATCHES=0
NUM_POINTS=100
MAX_IT=$NUM_BATCHES
GAMMA_PISF=1
MAX_IT_PISF=300

START=1
END=1

EXPERIMENT="v_aaai_bj_"$N"_"$SR"_"$NA"_"$NUM_BATCHES"_"$NUM_EPISODES"_"$MIN_BATCHES"_"$NUM_POINTS"_"$MAX_IT"_"$GAMMA_PISF"_"$MAX_IT_PISF
PLOT_FILE="plot_"$EXPERIMENT".R"

echo "D <- NULL" > $PLOT_FILE
echo "legends=list()" >> $PLOT_FILE
echo >> $PLOT_FILE
for M in "10" "100"
do
    for ALPHA in "0.1" "0.4" "1"
    do
        for EPSILON in "0.1" "0.4" "1"
        do
            GROUP=$EXPERIMENT"_"$M"_"$EPSILON"_"$ALPHA
            GROUP_FILE=$GROUP".log"

            if [ -f $GROUP_FILE ]
            then
                rm $GROUP_FILE
            fi

            ALL_EXIST=1
            for RUN in $(seq $START $END)
            do
                RUN_FILE=$GROUP"_"$(printf %02d%s ${RUN%})".log"
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
    done
done

echo "matplot(D, t=\"o\")" >> $PLOT_FILE
echo "grid()" >> $PLOT_FILE
