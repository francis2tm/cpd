#!/bin/bash

#SBATCH --job-name=cpd-omp-mpi
#SBATCH --output=100k-2-4.out
#SBATCH --error=100k-2-4.err
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun time ./matFact-mpi 2 tests/instML100k.in