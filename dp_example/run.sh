#!/bin/bash

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 2 ] || die "Usage: run.sh START END"

for M in 10 50 100
do
    while [ "$(pidof blackjack)" ]
    do
	sleep 1
    done

    for i in $(eval echo {$1..$2})
    do
	COMMAND="./blackjack $i 30000 1000000 0 100 30 1.0 300 $M"
	echo $COMMAND
	$COMMAND &
    done
done
