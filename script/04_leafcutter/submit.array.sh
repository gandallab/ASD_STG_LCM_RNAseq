#!/bin/bash
#$ -l h_data=1G,h_rt=00:20:00
#$ -wd /u/project/gandalm/pampas/ASD_STG_LCM_RNAseq/junc
#$ -t 1-3

./bam2junc.jobarray.sh $SGE_TASK_ID