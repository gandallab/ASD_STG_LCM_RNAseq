#!/bin/bash
#$ -l h_data=1G,h_rt=00:20:00
#$ -wd /u/project/gandalm/pampas/ASD_STG_LCM_RNAseq/junc
#$ -t 1-3

source /u/home/p/pampas/.bashrc
INPUT_FILE=/u/project/gandalm/pampas/ASD_STG_LCM_RNAseq/junc/cell_bamList
INBAM=`head -n ${SGE_TASK_ID} ${INPUT_FILE} | tail -n 1`
tp=3
len=`expr ${#INBAM} - $tp`
OUTJUC=${INBAM:0:$len}junc

/u/home/p/pampas/local/bin/leafcutter/scripts/bam2junc.sh ${INBAM} ${OUTJUC}

touch ${INBAM}.done