#!/bin/bash -l
#SBATCH --job-name=run_GWAS
#SBATCH --mem-per-cpu=64000
#SBATCH --account=mignot
#SBATCH --time=12:00:00

module load python/3.6.4
python run_GWAS.py

