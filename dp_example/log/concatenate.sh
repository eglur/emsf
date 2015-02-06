#!/bin/bash


die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 3 ] || die "Usage: concatenate.sh PREFIX START END"

PREFIX=$1
echo $PREFIX

for SR in 20 30 40 50
do
    PREFIX_SR=$PREFIX"_"$SR
    PREFIX_SR_LOG=$PREFIX_SR".log"
    if [ -f $PREFIX_SR_LOG ]
    then
        rm $PREFIX_SR_LOG
    fi

    for i in $(eval echo {$2..$3})
    do
        PREFIX_SR_RUN_LOG=$PREFIX_SR"_"`printf %02d%s ${i%}`".log"
        (cat $PREFIX_SR_RUN_LOG; echo) >> $PREFIX_SR_LOG
    done
done
