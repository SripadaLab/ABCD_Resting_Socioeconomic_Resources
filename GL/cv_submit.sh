#!/bin/bash

START=$1
STOP=$2

SCRIPTDIR=/home/mangstad/work/ser
OUTPUTDIR=/home/mangstad/work/ser
LOG_DIR=${OUTPUTDIR}/logs/

for SEED in `seq ${START} ${STOP}`
do
    export SEED
    echo sbatch --job-name bs_${SEED} \
	--output=${LOG_DIR}/slurm_%j.log \
	${SCRIPTDIR}/cv.master
    sbatch --job-name bs_${SEED} \
        --output=${LOG_DIR}/slurm_%j.log \
	${SCRIPTDIR}/cv.master
done
