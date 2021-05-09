#!/bin/bash

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 3 ] || die "Usage: run.sh START END n"

for i in $(eval echo {$1..$2})
do
    COMMAND="./emsf $i $3"
    echo $COMMAND
    $COMMAND &
done

