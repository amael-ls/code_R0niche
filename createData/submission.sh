#!/bin/bash

#SBATCH --array=1-14
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=12288
#SBATCH --ntasks=1

#SBATCH --job-name=matlabData
#SBATCH --mail-user=amael.lesquin@usherbrooke.ca
#SBATCH --mail-type=ALL

echo `pwd`

echo $SLURM_ARRAY_TASK_ID

cd $SLURM_SUBMIT_DIR
Rscript ~/projects/def-dgravel/amael/article1/progToSendToReview/puteTG.R
