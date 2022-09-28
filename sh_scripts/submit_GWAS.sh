#!/bin/bash

#SBATCH --job-name=GWAS
#SBATCH --mem-per-cpu=16000
#SBATCH --output=GWAS.out
#SBATCH --error=GWAS.err
#SBATCH --account=mignot
#SBATCH --time=08:00:00

python3 run_GWAS.py