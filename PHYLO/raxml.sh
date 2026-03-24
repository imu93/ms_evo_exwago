#!/usr/bin/bash

#SBATCH --job-name=raxml-ng  # Job name
#SBATCH --export=ALL
#SBATCH --ntasks=1                   # CPU request
#SBATCH --cpus-per-task=60
#SBATCH --mem=280gb                     # Job memory request
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=ss3.%j.log   # Standard output and error log

source /home/s2234986/.miniforge3/etc/profile.d/conda.sh
conda activate phylo_env
conda info | grep "active environment"


aln=gbe_all_agos.aln
model=LG+I+G4
raxml-ng --all --msa $aln --model $model --prefix T1  --seed 2 --threads 60 --workers 1 --bs-metric fbp --bs-trees 100
