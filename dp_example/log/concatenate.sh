#!/bin/bash

ONELINES_DIR="./onelines"
if [ ! -d $ONELINES_DIR ]
then
    mkdir -p $ONELINES_DIR
fi

    COMMAND="./blackjack $i 10000 1000000 0 100"

NUM_BATCHES=100000
NUM_EPISODES=1000000
MIN_BATCHES=0
NUM_POINTS=100
N=203
SR=20
NA=2

START=1
END=50

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
    echo $PREFIX_RUN_LOG
    if [ -f $PREFIX_RUN_LOG ]
    then
        (cat $PREFIX_RUN_LOG; echo) >> $PREFIX_LOG
    else
        ALL_EXIST=0
    fi
done

if [ $ALL_EXIST -eq 1 ]
then
    echo "Generated    $PREFIX_LOG"

    echo "tmp <- read.table(\"$PREFIX_LOG\")" >> $PLOT_FILENAME
    echo "tmp <- apply(tmp, 2, mean)" >> $PLOT_FILENAME
    echo "D <- cbind(D, tmp)" >> $PLOT_FILENAME
    echo >> $PLOT_FILENAME
fi            

echo "matplot(D, t=\"o\", main=\"Counting + PI utilizando $NUM_BATCHES lotes (avaliando em $NUM_EPISODES jogos)\")" >> $PLOT_FILENAME
echo "grid()" >> $PLOT_FILENAME
echo >> $PLOT_FILENAME
echo
