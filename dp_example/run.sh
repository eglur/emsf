#!/bin/bash

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 2 ] || die "Usage: run.sh START END"

#############################
# Parâmetros do experimento #
#############################
NUM_EPISODES=1000000    # Quantidade de episódios

NUM_BATCHES=30000       # Quantidade final de batches utilizados
MIN_BATCHES=0           # Quantidade inicial de batches utilizados
NUM_POINTS=100          # Número total de pontos; o programa obtém internamente o incremento

MAX_IT=30               # Número máximo de iterações do EM-SF

GAMMA_PISF=1            # Valor do gamma utilizado no PISF
MAX_IT_PISF=300         # Número máximo de iterações do PISF

START=$1                # Número do experimento inicial (utilizado para diferenciar os nomes
END=$2                  # Número do experimento final (utilizado para diferenciar os nomes dos arquivos de logs)

CORES=10                # Quantidade de CPU cores da máquina que serão utilizados

# LOCAL_START E LOCAL_END são utilizados para limitar a utilização de
# cores de acordo com os parâmetros setados.
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

    for M in 5 15 20 25
    do
	while [ "$(pidof blackjack)" ]
	do
	    sleep 1
	done

	for i in $(eval echo {$LOCAL_START..$LOCAL_END})
	do
	    COMMAND="./blackjack $i $NUM_BATCHES $NUM_EPISODES $MIN_BATCHES $NUM_POINTS $MAX_IT $GAMMA_PISF $MAX_IT_PISF $M"
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
