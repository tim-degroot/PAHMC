#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH --tasks-per-node 36

module load 2022
module load Python/3.10.4-GCCcore-11.3.0
module load Anaconda3/2022.05



python src/main.py $1 36



