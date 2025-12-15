#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=32g
#SBATCH -t 07-00:00:00

module load matlab/2025a
matlab -nodesktop -nosplash -singleCompThread -r run_batch -logfile spp.out
matlab -nodesktop -nosplash -singleCompThread -r post -logfile post.out