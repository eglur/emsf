#!/bin/bash

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 2 ] || die "Usage: run.sh START END"

NUM_BATCHES=5000
NUM_EPISODES=1000000
MIN_BATCHES=0
NUM_POINTS=100
MAX_IT=$NUM_BATCHES
TC=100
GAMMA_PISF=0.999
MAX_IT_PISF=300

case "$HOSTNAME" in
    prjcuda02)
        NCORES=15
        ;;
    turing)
        NCORES=11
        ;;
    skye)
        NCORES=3
        ;;
    *)
        echo "Numero de cores desconhecido."
        exit 1
esac

START=$1
END=$2

PIDS=()
for RUN in $(seq $START $END)
do
    for M in "5" "10"
    do
        for ALPHA in "0.2" "0.4" "0.6" "0.8" "1.0"
        do
            for EPSILON in "0.15"
            do
                COMMAND="./blackjack $RUN $NUM_BATCHES $NUM_EPISODES $MIN_BATCHES $NUM_POINTS $MAX_IT $TC $GAMMA_PISF $MAX_IT_PISF $M $EPSILON $ALPHA"
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
    done
done
