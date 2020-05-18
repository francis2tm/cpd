#!/bin/bash

#SBATCH --job-name=cpd-omp-mpi
#SBATCH --output=1M-4-12.out
#SBATCH --error=1M-4-12.err
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=4

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun time ./matFact-mpi 4 tests/instML1M.in