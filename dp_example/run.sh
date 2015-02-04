#!/bin/bash

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 2 ] || die "Usage: run.sh START END"

for i in $(eval echo {$1..$2})
do
    COMMAND="./emsf 100 10 1 100000 10 1e20 30 $i"
    echo $COMMAND
    $COMMAND &
done

