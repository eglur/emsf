#!/bin/bash

NUM_BATCHES=30000
NUM_EPISODES=1000000
MIN_BATCHES=0
NUM_POINTS=100
N=203
SR=20
NA=2

START=1
END=15

PREFIX="v_cnt_bj_"$SR"_"$NA"_"$NUM_BATCHES"_"$NUM_EPISODES"_"$MIN_BATCHES"_"$NUM_POINTS

PLOT_FILENAME="plot_"$PREFIX".R"
echo "D <- NULL" > $PLOT_FILENAME
echo >> $PLOT_FILENAME

PREFIX_LOG=$PREFIX".log"
if [ -f $PREFIX_LOG ]
then
    rm $PREFIX_LOG
fi

ALL_EXIST=1
for RUN in $(eval echo {$START..$END})
do
    PREFIX_RUN_LOG=$PREFIX"_"`printf %02d%s ${RUN%}`".log"
    if [ -f $PREFIX_RUN_LOG ]
    then
        (cat $PREFIX_RUN_LOG; echo) >> $PREFIX_LOG
    else
        ALL_EXIST=0
    fi
done

if [ $ALL_EXIST -eq 1 ]
then
    echo "tmp <- read.table(\"$PREFIX_LOG\")" >> $PLOT_FILENAME
    echo "tmp <- apply(tmp, 2, mean)" >> $PLOT_FILENAME
    echo "D <- cbind(D, tmp)" >> $PLOT_FILENAME

    echo "SUCCESS on generating    $PLOT_FILENAME"
fi            

echo >> $PLOT_FILENAME
echo "matplot(D, t=\"o\")" >> $PLOT_FILENAME
echo "grid()" >> $PLOT_FILENAME
