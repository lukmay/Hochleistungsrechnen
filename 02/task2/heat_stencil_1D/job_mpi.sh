#!/bin/bash

#SBATCH --partition=lva
#SBATCH --job-name heat_stencil_mpi
#SBATCH --output heat_stencil_mpi.log
#SBATCH --ntasks 4
#SBATCH --cpus-per-task 1
#SBATCH --exclusive

module load gcc/12.2.0-gcc-8.5.0-p4pe45v
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55

mpiexec -n $SLURM_NTASKS ./heat_stencil_1D_mpi 400