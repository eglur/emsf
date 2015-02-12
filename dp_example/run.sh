#!/bin/bash

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 2 ] || die "Usage: run.sh START END"

CORES=14
START=$1
END=$2

LOCAL_START=$START
LOCAL_END=`expr $START + $CORES - 1`

if [ $LOCAL_END -gt $END ]
then
    LOCAL_END=$END
fi

while [ $LOCAL_START -le $END ]
do
    while [ "$(pidof blackjack)" ]
    do
	sleep 1
    done

    for M in 10 50 100
    do
	while [ "$(pidof blackjack)" ]
	do
	    sleep 1
	done

	for i in $(eval echo {$LOCAL_START..$LOCAL_END})
	do
	    COMMAND="./blackjack $i 30000 1000000 0 100 30 1.0 300 $M"
	    echo $COMMAND
	    $COMMAND &
	done
    done

    LOCAL_START=`expr $LOCAL_END + 1`
    LOCAL_END=`expr $LOCAL_START + $CORES - 1`

    if [ $LOCAL_END -gt $END ]
    then
	LOCAL_END=$END
    fi
done
