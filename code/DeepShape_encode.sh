#!/bin/bash
#
#SBATCH --job-name=DESSO
#SBATCH --output=res.txt
#
#SBATCH --time=99:00:00
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000
#SBATCH --gres gpu:10
#SBATCH --gres-flags=enforce-binding

srun python DeepShape_encode.py