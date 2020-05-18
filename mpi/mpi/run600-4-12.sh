#!/bin/bash

#SBATCH --job-name=cpd-omp-mpi
#SBATCH --output=600-4-12.out
#SBATCH --error=600-4-12.err
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=4

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun time ./matFact-mpi 4 tests/inst600-10000-10-40-400.in