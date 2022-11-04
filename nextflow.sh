#!/bin/bash
#SBATCH --output=./slurm-out/nextflowmain-%j.out
#SBATCH --job-name=nf_main
#SBATCH --mem=32000
#SBATCH -c 1
mkdir -p slurm-out

rm -rf work/; rm .nextflow.log*; rm -rf results/
find ./slurm-out/ -type f -not -name "*-${SLURM_JOBID}.out" -delete
nextflow run Duncan_Lab_Pipeline.nf --genome 'GRCm39' --rmrRNA