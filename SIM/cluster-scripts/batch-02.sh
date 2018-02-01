#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH --mem 4000
#SBATCH --output=./../cluster-out/02-%a.out
#SBATCH --error=./../cluster-err/02-%a.err
#SBATCH --array=1-546

## add SAS
module add sas/9.4


## run SAS command
sas -work /dev/shm -noterminal ./../programs/02-perform-simulations.sas -log "./../cluster-logs/02-$SLURM_ARRAY_TASK_ID.log" -sysparm "$SLURM_ARRAY_TASK_ID"