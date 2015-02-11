#!/bin/bash

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 2 ] || die "Usage: run.sh START END"

for i in $(eval echo {$1..$2})
do
    COMMAND="./blackjack $i 30000 1000000 0 100"
    echo $COMMAND
    $COMMAND &
done
