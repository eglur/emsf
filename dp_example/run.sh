#!/bin/bash

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 2 ] || die "Usage: run.sh START END"

NCORES=15
START=$1
END=$2

PIDS=()
for i in $(seq $START $END)
do
    for EPSILON in "0.01" "0.05" "0.10" "0.20" "0.30" "0.40" "0.50" "0.60" "0.70" "0.80" "0.90" "1.0"
    do
        COMMAND="./blackjack $i 30000 1000000 0 100 $EPSILON"
        echo $COMMAND
        $COMMAND &

        # Adiciona pid do emsf que acabou de ser executado à lista de emsfs em execução
        PIDS+=($!)

        # Somente continua se o tamanho da lista de pids em execução ainda não chegou ao limite (i.e., quantidade de cores)
        # Enquanto isso, checa se os pids da lista ainda estão em execução e a atualiza
        while [ ${#PIDS[@]} == $NCORES ]; do
            sleep 1
            # Para cada pid da lista
            for i in $(seq 0 $((${#PIDS[@]} - 1))); do
                if ! $(ps -p ${PIDS[$i]} &> /dev/null); then
                    # Remove pid da lista se processo já terminou sua execução
                    PIDS=(${PIDS[@]:0:$i} ${PIDS[@]:(($i +1))})
                    # Sair, pois não é sensato iterar por índice em uma lista que acabou de ser alterada
                    break
                fi
            done
        done
    done
done
