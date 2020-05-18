#!/bin/bash

#SBATCH --job-name=cpd-omp-mpi
#SBATCH --output=600-1-4.out
#SBATCH --error=600-1-4.err
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun time ./matFact-mpi 1 tests/inst600-10000-10-40-400.in