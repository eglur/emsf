#!/bin/bash

for M in "10" "50" "100"
do
    for TCC in "100" "50"
    do
        for ALPHA in "1.0" "0.5" "0.1"
        do
            echo "bj.emsf.comp.experiment(m=$M, alpha=$ALPHA, tcc=$TCC, num.episodes=5e3, epsilon=0.15, tc=100, num.avg=50, ne=1e6)"
        done
    done
done
