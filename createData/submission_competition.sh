#!/bin/bash

#SBATCH --array=1-5000%1000
#SBATCH -t 00:30:00
#SBATCH --mem-per-cpu=512
#SBATCH --ntasks=1
#SBATCH --job-name=competition
#SBATCH --mail-user=amael.lesquin@usherbrooke.ca
#SBATCH --mail-type=ALL

echo `pwd`

echo $SLURM_ARRAY_TASK_ID

cd $SLURM_SUBMIT_DIR
Rscript ~/scratch/progToSendToReview/createData/calculateCompetition.R
