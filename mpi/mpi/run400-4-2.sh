#!/bin/bash

#SBATCH --job-name=cpd-omp-mpi
#SBATCH --output=400-4-2.out
#SBATCH --error=400-4-2.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun time ./matFact-mpi 4 tests/inst400-50000-30-200-500.in