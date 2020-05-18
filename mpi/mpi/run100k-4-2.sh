#!/bin/bash

#SBATCH --job-name=cpd-omp-mpi
#SBATCH --output=100k-4-2.out
#SBATCH --error=100k-4-2.err
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun time ./matFact-mpi 4 tests/instML100k.in