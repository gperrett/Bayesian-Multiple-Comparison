#!/bin/bash
#
#SBATCH --job-name=multiBART
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --array=1-1000 # this creates an array!
#SBATCH --time=130:00
#SBATCH --mem=8GB
#SBATCH --output=run_%A_%a.out
#SBATCH --error=run_%A_%a.err
module purge
module load r/gcc/4.5.0

Rscript dgp1_study1.R
Rscript dgp2_study1.R
Rscript dgp1_study2.R
Rscript dgp2_study2.R
