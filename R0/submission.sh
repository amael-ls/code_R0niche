#!/bin/bash

#SBATCH --account=def-dgravel
#SBATCH --array=1-3
#SBATCH --time 00:05:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --job-name=matlab_integral
#SBATCH --mem=8192
#SBATCH --mail-user=amael.lesquin@usherbrooke.ca
#SBATCH --mail-type=ALL

echo `pwd`
echo $SLURM_ARRAY_TASK_ID

# Create a local work directory
# trick from https://www.rc.fas.harvard.edu/resources/documentation/software/matlab-pct-simultaneous-job-problem/
mkdir -p /scratch/$USER/$SLURM_JOB_ID

# Change to the good directory and run
cd $SLURM_SUBMIT_DIR

srun matlab -nodisplay -r "main"

# Cleanup local work directory
rm -rf /scratch/$USER/$SLURM_JOB_ID
