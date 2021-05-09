#!/bin/bash


die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 2 ] || die "Usage: concatenate.sh START END"


FILENAME_ROOT="t_cnt"
if [ -f $FILENAME_ROOT".log" ]
then
    rm $FILENAME_ROOT".log"
fi

for i in $(eval echo {$1..$2})
do
    FILENAME=$FILENAME_ROOT"_"`printf %02d%s ${i%}`".log"
    (cat $FILENAME; echo) >> $FILENAME_ROOT".log"
done


ROOT="t_emsf"
for l in h i j k
do
    LETTER=$ROOT"_"$l
    LETTERLOG=$LETTER".log"
    if [ -f $LETTERLOG ]
    then
b        rm $LETTERLOG
    fi

    for i in $(eval echo {$1..$2})
    do
        NUMBERLOG=$LETTER"_"`printf %02d%s ${i%}`".log"
        (cat $NUMBERLOG; echo) >> $LETTERLOG
    done
done
