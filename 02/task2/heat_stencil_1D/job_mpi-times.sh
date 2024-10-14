#!/bin/bash

#SBATCH --partition=lva
#SBATCH --job-name=heat_stencil_mpi
#SBATCH --output=heat_stencil_mpi.log
#SBATCH --ntasks=96
#SBATCH --cpus-per-task=1
#SBATCH --exclusive

module load gcc/12.2.0-gcc-8.5.0-p4pe45v
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55

# Define the problem sizes and number of ranks
problem_sizes=(1024 2048 4096 6144)
ranks=(2 4 8 16 32 64)

# Loop through both the problem sizes and number of ranks
for size in "${problem_sizes[@]}"
do
    for rank in "${ranks[@]}"
    do
        echo "Running MPI Program with problem size: $size and ranks: $rank" >> heat_stencil_mpi.log
        echo "Problem Size: $size, Ranks: $rank" >> heat_stencil_mpi.log
        { time mpiexec -n $rank ./heat_stencil_1D_mpi $size; } 2>> heat_stencil_mpi.log
        echo "---------------------------------------------" >> heat_stencil_mpi.log
    done
done
