#!/bin/bash

#SBATCH --partition=lva
#SBATCH --job-name pi_mpi
#SBATCH --output pi_mpi_70.log
#SBATCH --ntasks 70
#SBATCH --cpus-per-task 1
#SBATCH --exclusive

module load gcc/12.2.0-gcc-8.5.0-p4pe45v
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55

mpiexec -n $SLURM_NTASKS ./pi_mpi 100
mpiexec -n $SLURM_NTASKS ./pi_mpi 500
mpiexec -n $SLURM_NTASKS ./pi_mpi 1000
mpiexec -n $SLURM_NTASKS ./pi_mpi 5000
mpiexec -n $SLURM_NTASKS ./pi_mpi 10000
mpiexec -n $SLURM_NTASKS ./pi_mpi 50000
mpiexec -n $SLURM_NTASKS ./pi_mpi 100000
mpiexec -n $SLURM_NTASKS ./pi_mpi 500000
mpiexec -n $SLURM_NTASKS ./pi_mpi 1000000
mpiexec -n $SLURM_NTASKS ./pi_mpi 5000000
mpiexec -n $SLURM_NTASKS ./pi_mpi 10000000
mpiexec -n $SLURM_NTASKS ./pi_mpi 50000000
mpiexec -n $SLURM_NTASKS ./pi_mpi 100000000
mpiexec -n $SLURM_NTASKS ./pi_mpi 500000000
mpiexec -n $SLURM_NTASKS ./pi_mpi 1000000000