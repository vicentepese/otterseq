#!/bin/bash -l
#SBATCH --job-name=run_GWAS
#SBATCH --mem-per-cpu=64000
#SBATCH --account=mignot
#SBATCH --time=12:00:00

python run_GWAS.py

