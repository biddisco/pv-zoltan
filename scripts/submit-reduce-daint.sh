#!/bin/bash

#SBATCH --job-name=reduce
#SBATCH --time=01:00:00
#SBATCH --nodes=256
#SBATCH --constraint=startx
#SBATCH --partition=viz

export DISPLAY=:0.0
export MPICH_MAX_THREAD_SAFETY=multiple

module load boost/1.56.0
module load cudatoolkit
#module load viz

# boost
export LD_LIBRARY_PATH=/apps/daint/boost/1.56.0/gnu_482/lib:\$LD_LIBRARY_PATH
# hdf5
export LD_LIBRARY_PATH=/users/biddisco/apps/daint/hdf5_1_8_cmake/lib:\$LD_LIBRARY_PATH
# nvidia GL
export LD_LIBRARY_PATH=/opt/cray/nvidia/default/lib64:\$LD_LIBRARY_PATH
# paraview
export LD_LIBRARY_PATH=/scratch/daint/biddisco/egpgv/:\$LD_LIBRARY_PATH

# plugins
export LD_LIBRARY_PATH=/scratch/daint/biddisco/egpgv:\$LD_LIBRARY_PATH
export PV_PLUGIN_PATH=/scratch/daint/biddisco/egpgv

export OMP_NUM_THREADS=1

aprun -n 256 -N 1 /scratch/daint/biddisco/egpgv/pvbatch /scratch/daint/biddisco/egpgv/reduce_points.py -r $1
