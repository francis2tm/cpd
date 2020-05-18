#!/bin/bash

#SBATCH --job-name=cpd-omp-mpi
#SBATCH --output=30-4-2.out
#SBATCH --error=30-4-2.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun time ./matFact-mpi 4 tests/inst30-40-10-2-10.in